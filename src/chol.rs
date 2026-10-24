//! Supernodal numeric Cholesky of `A D A'` and the triangular solves of the reference's barrier, bit for bit (port of
//! `reverse-engineer/tools/barrier_chol.py` and `tools/barrier_solve.py`; specs `barrier-chol.md`, `barrier-solve.md`).
//!
//! The reference's own code (assembly `0x2f83e0`, dense `LDL'` `0x2f1940`, the solves `0x305f10`) uses SSE2 without FMA; its DGEMM, DGEMV and
//! DTRSM are the AVX2 kernels of the statically linked MKL, which use FMA. Both are reproduced with explicit operation order: `mul_add` where
//! the reference fuses, plain `a * b + c` (two roundings) elsewhere. Nothing here may be rewritten in a form that changes the order of a sum
//! (no iterator sums, no auto-vectorised reductions: Rust never contracts `a * b + c` into an fma by itself).
//!
//! Gaps (never approximated): regularisation and 2x2 pivots of `0x2f1940` (never triggered in the corpus; this port has the normal pivot test
//! only), dense columns, free variables, the iterative refinement of `0x30fc90`.
use crate::barrier::{
    BarrierError, DirKind, Direction, NormalSolver, Rhs, StartKind, State, StdForm,
};
use crate::dump;
use crate::order::Symbolic;
use serde_json::json;

#[inline(always)]
fn fma(a: f64, b: f64, c: f64) -> f64 {
    a.mul_add(b, c)
}

const NONE: usize = usize::MAX;

// ---------------------------------------------------------------------------------------------------------------------------------
// Supernode plan: fronts, update lists, workspace size.

/// One supernode of the symbolic factorization: columns `col..col + w`, `c` rows in its list (the first `w` are the diagonal block),
/// and the earlier supernodes (indices, ascending) that update it.
#[derive(Debug, Clone)]
pub struct Supernode {
    pub col: usize,
    pub w: usize,
    pub c: usize,
    pub rows: Vec<usize>,
    pub upd: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct Plan {
    pub sns: Vec<Supernode>,
    /// `SYM+0x80` = max(480, n, widest supernode squared) (`0x2e8680`): the workspace size that bounds the row blocks of a child update.
    pub sym80: usize,
    pub n: usize,
}

impl Plan {
    pub fn new(sym: &Symbolic) -> Plan {
        let n = sym.perm.len();
        let mut sns: Vec<Supernode> = Vec::new();
        for j in 0..n {
            if sym.width[j] > 0 {
                // Numeric supernodes (SYM width2, `0x302277..0x3023e4`): a symbolic supernode of width W is cut into chunks of at most 192 columns (240 when
                // the chunk's remaining row count is >= 0x780); a later chunk keeps the tail of the row list. The other cap (`0x302277`, r9d: a function of
                // the flops and SYM+0x48/0x60, at least 8) is >= 240 on every captured model and not decoded.
                let (c0, roff, wt) = (
                    sym.count[j] as usize,
                    sym.rowoff[j] as usize,
                    sym.width[j] as usize,
                );
                let mut off = 0;
                while off < wt {
                    let cap = if c0 - off >= 0x780 { 240 } else { 192 };
                    let w = (wt - off).min(cap);
                    sns.push(Supernode {
                        col: j + off,
                        w,
                        c: c0 - off,
                        rows: sym.idx[roff + off..roff + c0]
                            .iter()
                            .map(|&r| r as usize)
                            .collect(),
                        upd: Vec::new(),
                    });
                    off += w;
                }
            }
        }
        // sn_of[column] = supernode index.
        let mut sn_of = vec![NONE; n];
        for (k, s) in sns.iter().enumerate() {
            for t in 0..s.w {
                sn_of[s.col + t] = k;
            }
        }
        // A supernode s is updated by every earlier supernode d (ascending) with an off-diagonal row in [s.col, s.col + w_s).
        let mut targets: Vec<Vec<usize>> = vec![Vec::new(); sns.len()];
        for (k, d) in sns.iter().enumerate() {
            let mut t: Vec<usize> = d.rows[d.w..].iter().map(|&r| sn_of[r]).collect();
            t.sort_unstable();
            t.dedup();
            targets[k] = t;
        }
        for (k, t) in targets.iter().enumerate() {
            for &s in t {
                if s > k {
                    sns[s].upd.push(k);
                }
            }
        }
        let maxw = sns.iter().map(|s| s.w).max().unwrap_or(0);
        let sym80 = 480.max(n).max(maxw * maxw);
        Plan { sns, sym80, n }
    }
}

// ---------------------------------------------------------------------------------------------------------------------------------
// Assembly of A D A' (0x2f83e0, generic 0x2f77f0).

/// Row-wise view of the column storage: `rows[i]` = `(column, position in the column storage)` in column order.
pub fn model_rows(f: &StdForm) -> Vec<Vec<(usize, usize)>> {
    let mut rows = model_rows_all(f);
    if !f.dense.is_empty() {
        // The dense columns are nodes of the factor, not terms of `A D A'`.
        for r in rows.iter_mut() {
            r.retain(|&(j, _)| !f.dense.contains(&j));
        }
    }
    rows
}

/// [`model_rows`] with every column (the products `A w` read the dense columns too).
fn model_rows_all(f: &StdForm) -> Vec<Vec<(usize, usize)>> {
    let mut rows = vec![Vec::new(); f.m];
    for j in 0..f.n2 {
        for q in f.col_start[j]..f.col_start[j + 1] {
            rows[f.row_idx[q]].push((j, q));
        }
    }
    rows
}

#[inline(always)]
fn add_m(mm: &mut [f64], r: &mut [f64], w: usize, a: usize, mq: usize, v: f64) {
    if mq < w {
        let k = mq * w + a;
        mm[k] = v + mm[k];
    } else {
        let k = (mq - w) * w + a;
        r[k] = v + r[k];
    }
}

/// `0x2f83e0` for one supernode (columns `col..col + w`, `c` rows in its list, `mp` = global row -> local index of the list, `NONE` for other rows):
/// `A D A'` plus the slack `D` of the non-'=' rows into `M` (`w x w` row-major, lower part) and `R` (`(c - w) x w` row-major).
/// Columns of `A` are visited in the order (first row inside the supernode, column index); the entries of a column inside the supernode are taken
/// three at a time (column a gets `(A_a d) A_q`, b `(A_b d) A_q`, cc `(A_cc d) A_q`), the rest one at a time with `(d A_j) A_q`.
pub fn assemble(
    f: &StdForm,
    rows: &[Vec<(usize, usize)>],
    d: &[f64],
    mp: &[usize],
    col: usize,
    w: usize,
    c: usize,
) -> (Vec<f64>, Vec<f64>) {
    let (idx, val, st) = (&f.row_idx, &f.val, &f.col_start);
    let (n1, n2) = (f.m, f.n2);
    let mut mm = vec![0.0; w * w];
    let mut r = vec![0.0; (c - w) * w];
    for i in col..col + w {
        for &(cj, pos) in &rows[i] {
            let s = st[cj];
            let e = st[cj + 1];
            if !(pos == s || idx[pos - 1] < col) {
                continue;
            }
            let mut l24 = pos;
            loop {
                l24 += 1;
                if e <= l24 || mp[idx[l24]] >= w {
                    break;
                }
            }
            let dd = d[cj];
            let mut j = pos;
            let sj = f.scale.as_ref().map_or(0.0, |v| v[cj]);
            if sj != 0.0 {
                // Scaled column (`Q+0x648[cj] != 0`): every contribution is the number (s s) d (triples) or s (d s) (the rest), added without a product with A.
                let (sv3, sv1) = ((sj * sj) * dd, sj * (dd * sj));
                while j + 2 < l24 {
                    let (a, b, cc) = (mp[idx[j]], mp[idx[j + 1]], mp[idx[j + 2]]);
                    mm[a * w + a] = sv3 + mm[a * w + a];
                    mm[b * w + a] = sv3 + mm[b * w + a];
                    mm[b * w + b] = sv3 + mm[b * w + b];
                    for q in j + 2..l24 {
                        let mq = mp[idx[q]];
                        add_m(&mut mm, &mut r, w, a, mq, sv3);
                        add_m(&mut mm, &mut r, w, b, mq, sv3);
                        add_m(&mut mm, &mut r, w, cc, mq, sv3);
                    }
                    for q in l24..e {
                        let mq = mp[idx[q]];
                        if mq < n1 {
                            add_m(&mut mm, &mut r, w, a, mq, sv3);
                            add_m(&mut mm, &mut r, w, b, mq, sv3);
                            add_m(&mut mm, &mut r, w, cc, mq, sv3);
                        }
                    }
                    j += 3;
                }
                while j < l24 {
                    let a = mp[idx[j]];
                    mm[a * w + a] = sv1 + mm[a * w + a];
                    for q in j + 1..l24 {
                        add_m(&mut mm, &mut r, w, a, mp[idx[q]], sv1);
                    }
                    for q in l24..e {
                        let mq = mp[idx[q]];
                        if mq < n1 {
                            add_m(&mut mm, &mut r, w, a, mq, sv1);
                        }
                    }
                    j += 1;
                }
                continue;
            }
            while j + 2 < l24 {
                let (a, b, cc) = (mp[idx[j]], mp[idx[j + 1]], mp[idx[j + 2]]);
                let (a0, a1, a2) = (val[j], val[j + 1], val[j + 2]);
                let (t0, t1, t2) = (a0 * dd, a1 * dd, a2 * dd);
                mm[a * w + a] = a0 * t0 + mm[a * w + a];
                mm[b * w + a] = a1 * t0 + mm[b * w + a];
                mm[b * w + b] = a1 * t1 + mm[b * w + b];
                for q in j + 2..l24 {
                    let (mq, v) = (mp[idx[q]], val[q]);
                    add_m(&mut mm, &mut r, w, a, mq, t0 * v);
                    add_m(&mut mm, &mut r, w, b, mq, t1 * v);
                    add_m(&mut mm, &mut r, w, cc, mq, t2 * v);
                }
                for q in l24..e {
                    let (mq, v) = (mp[idx[q]], val[q]);
                    if mq < n1 {
                        add_m(&mut mm, &mut r, w, a, mq, t0 * v);
                        add_m(&mut mm, &mut r, w, b, mq, t1 * v);
                        add_m(&mut mm, &mut r, w, cc, mq, t2 * v);
                    }
                }
                j += 3;
            }
            while j < l24 {
                let (a, v0) = (mp[idx[j]], val[j]);
                let t = dd * v0;
                mm[a * w + a] = v0 * t + mm[a * w + a];
                for q in j + 1..l24 {
                    add_m(&mut mm, &mut r, w, a, mp[idx[q]], val[q] * t);
                }
                for q in l24..e {
                    let mq = mp[idx[q]];
                    if mq < n1 {
                        add_m(&mut mm, &mut r, w, a, mq, val[q] * t);
                    }
                }
                j += 1;
            }
        }
        if f.sense[i] != b'=' {
            let k = (i - col) * w + (i - col);
            mm[k] = mm[k] + d[n2 + i];
        }
    }
    (mm, r)
}

/// `0x2f77f0` (used for `w < 3`, `c * w > 30000`): for every supernode column `i` and every column `cj` of row `i` (ascending), `t = D[cj] A[i][cj]`,
/// then `acc[row_q] = A[q] t + acc[row_q]` over the entries of `cj` from `i`'s position on; column `i` of the front is `acc` (the diagonal gets the slack `D`).
pub fn assemble_generic(
    f: &StdForm,
    rows: &[Vec<(usize, usize)>],
    d: &[f64],
    mp: &[usize],
    col: usize,
    w: usize,
    c: usize,
) -> (Vec<f64>, Vec<f64>) {
    let (idx, val, st) = (&f.row_idx, &f.val, &f.col_start);
    let n2 = f.n2;
    let mut mm = vec![0.0; w * w];
    let mut r = vec![0.0; (c - w) * w];
    let mut acc = vec![0.0; c];
    let mut seen = vec![false; c];
    let mut touched: Vec<usize> = Vec::new();
    for i in col..col + w {
        let a = i - col;
        for &(cj, pos) in &rows[i] {
            let e = st[cj + 1];
            let sj = f.scale.as_ref().map_or(0.0, |v| v[cj]);
            // Scaled column (`0x2f81fd`): every contribution is the number (s s) D[cj].
            let t = if sj != 0.0 {
                (sj * sj) * d[cj]
            } else {
                d[cj] * val[pos]
            };
            for q in pos..e {
                let rr = idx[q];
                let m = mp[rr];
                if m < c && rr >= i {
                    if !seen[m] {
                        seen[m] = true;
                        touched.push(m);
                    }
                    acc[m] = if sj != 0.0 {
                        t + acc[m]
                    } else {
                        val[q] * t + acc[m]
                    };
                }
            }
        }
        for &m in &touched {
            if m < w {
                mm[m * w + a] = acc[m];
            } else {
                r[(m - w) * w + a] = acc[m];
            }
            acc[m] = 0.0;
            seen[m] = false;
        }
        touched.clear();
        if f.sense[i] != b'=' {
            mm[a * w + a] = mm[a * w + a] + d[n2 + i];
        }
    }
    (mm, r)
}

// ---------------------------------------------------------------------------------------------------------------------------------
// MKL kernels.

/// One column of DGEMV('T') (`mkl_blas_avx2_xdgemv_t`): `y = fma(sum, alpha, y)`, `sum = (acc0 + acc1) + (acc2 + acc3)`, four lane accumulators.
/// 8-blocks: `t = fma(x_l, a_l, x_{l+4} a_{l+4})` (the upper half is rounded first), `acc_l = t + acc_l`; a 4-block: variant 4 `acc_l = fma(x, a, acc_l)`,
/// variants 2 and 1 `acc_l = x a + acc_l`; the last `L % 4` entries by masked fma on lanes `0..r`.
fn dgemv_t_col(a: &[f64], x: &[f64], variant: u8, alpha: f64, y: f64) -> f64 {
    let l = a.len();
    let mut acc = [0.0f64; 4];
    let main = l - (l % 8);
    let tail4 = (l % 8) & !3;
    let r = l % 4;
    let mut c = 0;
    while c < main {
        for k in 0..4 {
            let t = fma(x[c + k], a[c + k], x[c + 4 + k] * a[c + 4 + k]);
            acc[k] = t + acc[k];
        }
        c += 8;
    }
    let mut p = main;
    if tail4 != 0 {
        for k in 0..4 {
            acc[k] = if variant == 4 {
                fma(x[p + k], a[p + k], acc[k])
            } else {
                x[p + k] * a[p + k] + acc[k]
            };
        }
        p += 4;
    }
    for k in 0..r {
        acc[k] = fma(x[p + k], a[p + k], acc[k]);
    }
    let s = (acc[0] + acc[1]) + (acc[2] + acc[3]);
    fma(s, alpha, y)
}

/// DGEMV('T', L, n, alpha, A, lda, x, 1, beta = 1, y, incy): column `j` of `A` against `x`; four columns at a time, then a pair, then one.
fn dgemv_t(
    a: &[f64],
    lda: usize,
    l: usize,
    n: usize,
    x: &[f64],
    alpha: f64,
    y: &mut [f64],
    incy: usize,
) {
    let n4 = n - (n % 4);
    for j in 0..n {
        let var = if j < n4 {
            4
        } else if (n % 4) & 2 != 0 && j < n4 + 2 {
            2
        } else {
            1
        };
        y[j * incy] = dgemv_t_col(&a[j * lda..j * lda + l], x, var, alpha, y[j * incy]);
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    S,
    A,
    N,
}

/// One element of DGEMM('T','N') beta = 1 in `dgemm_kernel_nocopy_TN_b1` (alpha = 1): 'S' one sequential fma chain; 'A' even and odd fma chains
/// (`k - k % 2` entries) added, an odd last entry fma'd into the sum; 'N' as 'A' with rounded products. Then `c0 + s`.
fn dgemm_tn_elem(kind: Kind, a: &[f64], b: &[f64], c0: f64) -> f64 {
    let k = a.len();
    if kind == Kind::S {
        let mut s = 0.0;
        for i in 0..k {
            s = fma(a[i], b[i], s);
        }
        return c0 + 1.0 * s;
    }
    let (mut e, mut o) = (0.0, 0.0);
    for i in 0..k - (k % 2) {
        if i % 2 == 0 {
            e = if kind == Kind::A {
                fma(a[i], b[i], e)
            } else {
                e + a[i] * b[i]
            };
        } else {
            o = if kind == Kind::A {
                fma(a[i], b[i], o)
            } else {
                o + a[i] * b[i]
            };
        }
    }
    let mut s = e + o;
    if k % 2 == 1 {
        s = if kind == Kind::A {
            fma(a[k - 1], b[k - 1], s)
        } else {
            s + a[k - 1] * b[k - 1]
        };
    }
    c0 + 1.0 * s
}

/// Tile kernel of element (row `i`, column `j`) of an `m x n` product: rows in blocks of 8 (a last block of <= 4 rows is a 4-row block, 5..7 rows a
/// masked 8-row block), columns in blocks of 6 with a remainder of `n % 6`; remainder 1 or 2: 'N'; 3: 'A'; otherwise 8-row blocks 'S', 4-row 'A'.
fn gemm_kind(m: usize, n: usize, i: usize, j: usize) -> Kind {
    let rb = (i / 8) * 8;
    let h = if m - rb > 4 { 8 } else { 4 };
    let nfull = n - (n % 6);
    let c = if j < nfull { 6 } else { n % 6 };
    if c == 1 || c == 2 {
        return Kind::N;
    }
    if c == 3 {
        return Kind::A;
    }
    if h == 8 {
        Kind::S
    } else {
        Kind::A
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Acc {
    Chain,
    Par,
    Lane4,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum ColK {
    G4,
    C2,
    C1,
}

/// Accumulator of the 12-row-tile kernels (`dgemm_kernel_0`, `dtrsm_kernel_ll_0`): terms by `term(j)` for `j < k`, fma accumulation from 0.
/// `Chain` one accumulator; `Par` accumulators of the even and the odd terms over the first `4 (k / 4)` terms, their sum, then the remaining terms
/// fma'd into the sum; `Lane4` accumulators of the terms `j % 4`, `(a0 + a1) + (a2 + a3)`, then the tail fma'd into the sum.
fn acc_terms(mode: Acc, k: usize, term: &dyn Fn(usize) -> (f64, f64)) -> f64 {
    let m4 = k - k % 4;
    let mut acc;
    match mode {
        Acc::Chain => {
            acc = 0.0;
            for j in 0..k {
                let (x, y) = term(j);
                acc = fma(x, y, acc);
            }
            return acc;
        }
        Acc::Par => {
            let (mut e, mut o) = (0.0, 0.0);
            for j in 0..m4 {
                let (x, y) = term(j);
                if j % 2 == 0 {
                    e = fma(x, y, e);
                } else {
                    o = fma(x, y, o);
                }
            }
            acc = e + o;
        }
        Acc::Lane4 => {
            let mut a = [0.0f64; 4];
            for j in 0..m4 {
                let (x, y) = term(j);
                a[j % 4] = fma(x, y, a[j % 4]);
            }
            acc = (a[0] + a[1]) + (a[2] + a[3]);
        }
    }
    for j in m4..k {
        let (x, y) = term(j);
        acc = fma(x, y, acc);
    }
    acc
}

/// Accumulator structure per (row tile, column kernel) of `dgemm_kernel_0` / `dtrsm_kernel_ll_0` (same code template).
fn tile_mode(tile: usize, ck: ColK) -> Acc {
    use Acc::*;
    use ColK::*;
    match (tile, ck) {
        (12, _) => Chain,
        (8, G4) | (8, C2) => Chain,
        (8, C1) => Lane4,
        (4, G4) => Par,
        (4, C2) => Chain,
        (4, C1) => Lane4,
        (2, G4) | (2, C2) => Chain,
        (2, C1) => Par,
        _ => Par, // tile 1
    }
}

/// Column kernel of column `j` of `n` in the 12-row-tile kernels: groups of 4 columns, then a 2-column and a 1-column remainder.
fn col_kernel(n: usize, j: usize) -> ColK {
    let (q, r) = (n / 4, n % 4);
    if j < 4 * q {
        return ColK::G4;
    }
    match r {
        1 => ColK::C1,
        2 => ColK::C2,
        _ => {
            if j < 4 * q + 2 {
                ColK::C2
            } else {
                ColK::C1
            }
        }
    }
}

/// `mkl_blas_avx2_xdgemm` (`0xa74800`) for DGEMM('T','N', alpha != 0, beta != 0): true = copy-based `dgemm_kernel_0`, false = the nocopy driver.
fn gemm_path2(m: usize, n: usize, k: usize) -> bool {
    if k > 128 || (m * n) / (m + n) > 10 {
        return m >= 56 && n >= 56;
    }
    if m <= 100 {
        return false;
    }
    n > 3 || k > 3
}

/// Height of the row tile that holds row `i` of `m`: tiles of 12 rows, then the remainder `m % 12` split by its bits 8, 4, 2, 1.
fn k0_row_tile(m: usize, i: usize) -> usize {
    let mut r0 = 0;
    let mut tiles = vec![12; m / 12];
    for t in [8usize, 4, 2, 1] {
        if (m % 12) & t != 0 {
            tiles.push(t);
        }
    }
    for t in tiles {
        if i < r0 + t {
            return t;
        }
        r0 += t;
    }
    unreachable!()
}

/// Column blocks of the nocopy driver (`0xa7be50`): remaining >= 192: 96; > 144: ceil(rem / 2); else rem.
fn nocopy_n_blocks(n: usize) -> Vec<usize> {
    let mut out = Vec::new();
    let mut rem = n;
    while rem > 0 {
        let b = if rem >= 192 {
            96
        } else if rem > 144 {
            (rem + 1) / 2
        } else {
            rem
        };
        out.push(b);
        rem -= b;
    }
    out
}

/// DGEMM('T','N', m, n, k, alpha = 1, A, lda, B, ldb, beta = 1, C, ldc): `C(i,j) += sum_k A(k,i) B(k,j)`, Fortran column-major (`A(k,i) = a[i lda + k]`).
/// `m == 1` is DGEMV('T') of `B` against the row, `n == 1` of `A` against the column; everything else goes to `gemm_path2` (copy-based
/// 12-row tile kernel or the nocopy TN kernel per column block).
fn dgemm_t_n(
    m: usize,
    n: usize,
    k: usize,
    a: &[f64],
    lda: usize,
    b: &[f64],
    ldb: usize,
    c: &mut [f64],
    ldc: usize,
) {
    if m == 1 {
        dgemv_t(b, ldb, k, n, &a[..k], 1.0, c, ldc);
        return;
    }
    if n == 1 {
        dgemv_t(a, lda, k, m, &b[..k], 1.0, c, 1);
        return;
    }
    if gemm_path2(m, n, k) {
        for j in 0..n {
            let ck = col_kernel(n, j);
            for i in 0..m {
                let mode = tile_mode(k0_row_tile(m, i), ck);
                let (ar, br) = (&a[i * lda..i * lda + k], &b[j * ldb..j * ldb + k]);
                let s = acc_terms(mode, k, &|t| (1.0 * ar[t], br[t]));
                c[j * ldc + i] = c[j * ldc + i] + s;
            }
        }
        return;
    }
    let mut j0 = 0;
    for nb in nocopy_n_blocks(n) {
        for j in 0..nb {
            for i in 0..m {
                let kind = gemm_kind(m, nb, i, j);
                c[(j0 + j) * ldc + i] = dgemm_tn_elem(
                    kind,
                    &a[i * lda..i * lda + k],
                    &b[(j0 + j) * ldb..(j0 + j) * ldb + k],
                    c[(j0 + j) * ldc + i],
                );
            }
        }
        j0 += nb;
    }
}

const TRSM_LL0_MIN: usize = 13; // m >= 13, except m % 4 == 0 && m <= 88 (dtrsm_left_ker): dtrsm_kernel_ll_0
const TRSM_LL0_BLOCK: usize = 12;

/// DTRSM with a single right-hand-side column (`n == 1`) is routed by `mkl_blas_dtrsm` to `mkl_blas_avx2_dtrsv_utu` (0x1ced5b0, incx = 1):
/// rows in pairs `(p, p+1)` (a last odd row alone). For a row with `p` solved terms: `p >= 16`: blocks of 16 terms; the single row has four
/// ymm accumulators (fused, block `s` of four terms into accumulator `s`, lane 0 starts at b), a pair has two for row `p` (fused, blocks 0, 2 and
/// 1, 3) and one for row `p+1` (rounded products, subtracted in order); then the remaining quads (rounded products, one accumulator, lane 0 starts
/// at the value), lanes summed `(l0 + l2) + (l1 + l3)`, then the remaining terms one by one with rounded products. Row `p+1` of a pair finally
/// takes the term of row `p` fused.
fn dtrsv_utu_block(l: &[f64], w: usize, b: &[f64]) -> Vec<f64> {
    let mut x = b.to_vec();
    let lanes_sum = |a: [f64; 4]| (a[0] + a[2]) + (a[1] + a[3]);
    // value of row `i` over the terms j < p (row `r` of the pair, `r = i` itself), `single`: the 4-accumulator form
    let reduce = |x: &[f64], i: usize, p: usize, b0: f64, mode: u8| -> f64 {
        // mode 0: single row, 1: first row of a pair, 2: second row of a pair
        let mut v = b0;
        let mut j0 = 0;
        if p >= 16 {
            let nb = p / 16;
            let mut acc = [[0.0f64; 4]; 4];
            acc[0][0] = b0;
            if mode == 2 {
                for c in 0..nb {
                    for s in 0..4 {
                        for k in 0..4 {
                            let j = 16 * c + 4 * s + k;
                            acc[0][k] = acc[0][k] - l[i * w + j] * x[j];
                        }
                    }
                }
                v = lanes_sum(acc[0]);
            } else {
                for c in 0..nb {
                    for s in 0..4 {
                        let a = if mode == 0 { s } else { s % 2 };
                        for k in 0..4 {
                            let j = 16 * c + 4 * s + k;
                            acc[a][k] = fma(-l[i * w + j], x[j], acc[a][k]);
                        }
                    }
                }
                v = if mode == 0 {
                    let t: Vec<f64> = (0..4)
                        .map(|k| (acc[0][k] + acc[1][k]) + (acc[2][k] + acc[3][k]))
                        .collect();
                    (t[0] + t[2]) + (t[1] + t[3])
                } else {
                    let t: Vec<f64> = (0..4).map(|k| acc[0][k] + acc[1][k]).collect();
                    (t[0] + t[2]) + (t[1] + t[3])
                };
            }
            j0 = 16 * nb;
        }
        let nq = (p - j0) / 4;
        if nq != 0 {
            let mut lanes = [v, 0.0, 0.0, 0.0];
            for c in 0..nq {
                for k in 0..4 {
                    let j = j0 + 4 * c + k;
                    lanes[k] = lanes[k] - l[i * w + j] * x[j];
                }
            }
            v = lanes_sum(lanes);
            j0 += 4 * nq;
        }
        for j in j0..p {
            v = v - l[i * w + j] * x[j];
        }
        v
    };
    let mut p = 0;
    while p + 1 < w {
        let (v0, v1) = if p > 0 {
            (reduce(&x, p, p, x[p], 1), reduce(&x, p + 1, p, x[p + 1], 2))
        } else {
            (x[p], x[p + 1])
        };
        x[p] = v0;
        x[p + 1] = fma(-v0, l[(p + 1) * w + p], v1);
        p += 2;
    }
    if p < w && p > 0 {
        x[p] = reduce(&x, p, p, x[p], 0);
    }
    x
}

/// `mkl_blas_avx2_xdtrsv` ('U', 'T', 'U') cuts the rows into blocks of 32: block `k >= 1` first takes `y = y - A[0..32k, block]' x[0..32k]`
/// (`mkl_blas_dgemv`, `dgemv_t`), then the block's own triangle is solved by `dtrsv_utu_block`.
fn dtrsv_utu_unit(l: &[f64], w: usize, b: &[f64]) -> Vec<f64> {
    let mut x = b.to_vec();
    let mut r0 = 0;
    while r0 < w {
        let s = (w - r0).min(32);
        if r0 > 0 {
            let (done, rest) = x.split_at_mut(r0);
            dgemv_t(&l[r0 * w..], w, r0, s, done, -1.0, &mut rest[..s], 1);
        }
        let sub: Vec<f64> = (0..s)
            .flat_map(|i| (0..s).map(move |j| (i, j)))
            .map(|(i, j)| l[(r0 + i) * w + r0 + j])
            .collect();
        let blk = dtrsv_utu_block(&sub, s, &x[r0..r0 + s]);
        x[r0..r0 + s].copy_from_slice(&blk);
        r0 += s;
    }
    x
}

/// DTRSM('L','U','T','U') on right-hand-side column `col` of `n` (B has ldb = w), `A = L` the unit lower triangle of the row-major `w x w` block:
/// the solution column `L^-1 b`. `w >= 13` unless `w % 4 == 0 && w <= 88`: `dtrsm_kernel_ll_0` (rows 0..11 an fnma chain from b; then tiles of 12, 4, 2, 1 rows
/// whose update from all solved rows is an fma accumulator from 0 (structure per tile and column kernel), `x_i = b_i + u_i`, then the in-tile chain).
/// Smaller `w`: `mkl_blas_avx2_dtrsm_left` (`dtrsm_left_ker`, see the three branches).
fn dtrsm_llu_col(l: &[f64], w: usize, b: &[f64], col: usize, n: usize) -> Vec<f64> {
    if n == 1 {
        return dtrsv_utu_unit(l, w, b);
    }
    let mut x = b.to_vec();
    if w >= TRSM_LL0_MIN && !(w % 4 == 0 && w <= 88) {
        let ck = col_kernel(n, col);
        for i in 0..TRSM_LL0_BLOCK {
            let mut v = x[i];
            for j in 0..i {
                v = fma(-l[i * w + j], x[j], v);
            }
            x[i] = v;
        }
        let mut r0 = TRSM_LL0_BLOCK;
        let mut rows_left = w - TRSM_LL0_BLOCK;
        let mut tiles = vec![12usize; rows_left / 12];
        rows_left %= 12;
        for t in [8usize, 4, 2, 1] {
            if rows_left & t != 0 {
                tiles.push(t);
            }
        }
        for t in tiles {
            let mode = tile_mode(t, ck);
            for i in r0..r0 + t {
                let acc = acc_terms(mode, r0, &|j| (-l[i * w + j], x[j]));
                x[i] = x[i] + acc;
            }
            for i in r0..r0 + t {
                let mut v = x[i];
                for j in r0..i {
                    v = fma(-l[i * w + j], x[j], v);
                }
                x[i] = v;
            }
            r0 += t;
        }
        return x;
    }
    if col >= 4 * (n / 4) {
        for i in 0..w {
            let mut v = b[i] * 1.0;
            let mut j0 = 0;
            if i >= 16 {
                // 16-wide loop: four ymm accumulators (j = 16 c + 4 a + l), fused, the first starts at [v, 0, 0, 0].
                let mut acc = [[v, 0.0, 0.0, 0.0], [0.0; 4], [0.0; 4], [0.0; 4]];
                for c in 0..i / 16 {
                    for a in 0..4 {
                        for lane in 0..4 {
                            let j = 16 * c + 4 * a + lane;
                            acc[a][lane] = fma(-l[i * w + j], x[j], acc[a][lane]);
                        }
                    }
                }
                let t: Vec<f64> = (0..4)
                    .map(|k| (acc[0][k] + acc[1][k]) + (acc[2][k] + acc[3][k]))
                    .collect();
                v = (t[0] + t[2]) + (t[1] + t[3]);
                j0 = 16 * (i / 16);
            }
            let nch = (i - j0) / 4;
            if nch != 0 {
                let mut lanes = [v, 0.0, 0.0, 0.0];
                for c in 0..nch {
                    for k in 0..4 {
                        lanes[k] = lanes[k] - l[i * w + j0 + 4 * c + k] * x[j0 + 4 * c + k];
                    }
                }
                v = (lanes[0] + lanes[2]) + (lanes[1] + lanes[3]);
            }
            for j in j0 + 4 * nch..i {
                v = v - l[i * w + j] * x[j];
            }
            x[i] = v;
        }
        return x;
    }
    let nb = w / 4;
    let r0 = 4 * nb;
    for i in 0..r0 {
        let mut v = x[i];
        for j in 0..i {
            v = fma(-l[i * w + j], x[j], v);
        }
        x[i] = v;
    }
    for i in r0..w {
        if r0 != 0 {
            let mut acc = [0.0f64; 4];
            for blk in 0..nb {
                for k in 0..4 {
                    acc[k] = l[i * w + 4 * blk + k] * x[4 * blk + k] + acc[k];
                }
            }
            let s = (acc[0] + acc[2]) + (acc[1] + acc[3]);
            x[i] = x[i] + (-s);
        }
    }
    for i in r0..w {
        let mut v = x[i];
        for j in r0..i {
            v = v - l[i * w + j] * x[j];
        }
        x[i] = v;
    }
    x
}

// ---------------------------------------------------------------------------------------------------------------------------------
// Dense LDL' of the diagonal block (0x2f1940, no FMA).

/// `M[i][j] + sum_{k<j} M[i][k] temp[k]` (no fma): `j < 8` sequentially from `M[i][j]`; otherwise eight lane accumulators over `k < j - (j & 7)`
/// (lane 0 starting at `M[i][j]`), `((a0 + a2) + (a4 + a6)) + ((a1 + a3) + (a5 + a7))`, then the tail sequentially.
fn dot_update(mm: &[f64], w: usize, i: usize, j: usize, temp: &[f64]) -> f64 {
    let acc0 = mm[i * w + j];
    if j < 8 {
        let mut a = acc0;
        for k in 0..j {
            a = a + mm[i * w + k] * temp[k];
        }
        return a;
    }
    let end = j - (j & 7);
    let mut a = [0.0f64; 8];
    a[0] = acc0;
    let mut k = 0;
    while k < end {
        for l in 0..8 {
            a[l] = a[l] + mm[i * w + k + l] * temp[k + l];
        }
        k += 8;
    }
    let s0 = (a[0] + a[2]) + (a[4] + a[6]);
    let s1 = (a[1] + a[3]) + (a[5] + a[7]);
    let mut t = s0 + s1;
    for k in end..j {
        t = t + mm[i * w + k] * temp[k];
    }
    t
}

/// The pivot test of `0x2f1940` (normal mode): `|orig| 1e-16 < |d|` keeps `d` with reciprocal `1 / d`, else `d = 0` and reciprocal 0.
fn pivot(dj: f64, orig: f64, reg: Option<(usize, usize)>) -> (f64, f64) {
    match reg {
        None => {
            if orig.abs() * 1e-16 < dj.abs() {
                (dj, 1.0 / dj)
            } else {
                (0.0, 0.0)
            }
        }
        // Regularising mode (`SYM+0x110 == 1`, dense columns): `tau = clamp(|orig| 1e-10, 1e-12, 1e-6)`; a pivot of magnitude `<= tau` becomes `tau`
        // (`-tau` when negative on a dense node, a column `>= n1`). `reg = (global column, n1)`.
        Some((gc, n1)) => {
            let tau = (orig.abs() * 1e-10).max(1e-12).min(1e-6);
            if tau < dj.abs() {
                (dj, 1.0 / dj)
            } else {
                let d = if gc < n1 || 0.0 <= dj { tau } else { -tau };
                (d, 1.0 / d)
            }
        }
    }
}

/// In-place `LDL'` of the `w x w` row-major block (left looking, column by column; `M[i][j] = (M[i][j] + sum_k L[i][k] temp[k]) / d_j` with
/// `temp[k] = (-L[j][k]) d[k]`); the diagonal becomes `d`. Returns `(d, reciprocals)`.
fn dense_ldl(
    mm: &mut [f64],
    w: usize,
    orig: &[f64],
    reg: Option<(usize, usize)>,
) -> (Vec<f64>, Vec<f64>) {
    let mut d = vec![0.0; w];
    let mut rcp = vec![0.0; w];
    let mut temp = vec![0.0; w];
    for j in 0..w {
        for i in j..w {
            if j > 0 {
                mm[i * w + j] = dot_update(mm, w, i, j, &temp);
            }
            if i == j {
                let (dj, rj) = pivot(mm[j * w + j], orig[j], reg.map(|(c0, n1)| (c0 + j, n1)));
                d[j] = dj;
                rcp[j] = rj;
                mm[j * w + j] = dj;
            } else {
                mm[i * w + j] = mm[i * w + j] * rcp[j];
            }
        }
        if j + 1 < w {
            for k in 0..=j {
                temp[k] = (-mm[(j + 1) * w + k]) * d[k];
            }
        }
    }
    (d, rcp)
}

/// The block has no nonzero below the diagonal: pivots only, `L = I`.
fn diag_only(
    mm: &mut [f64],
    w: usize,
    orig: &[f64],
    reg: Option<(usize, usize)>,
) -> (Vec<f64>, Vec<f64>) {
    let mut d = vec![0.0; w];
    let mut rcp = vec![0.0; w];
    for j in 0..w {
        let (dj, rj) = pivot(mm[j * w + j], orig[j], reg.map(|(c0, n1)| (c0 + j, n1)));
        d[j] = dj;
        rcp[j] = rj;
        mm[j * w + j] = dj;
    }
    (d, rcp)
}

/// After the updates: `LDL'` of `M`, then the rectangle `R` := `(L^-1 R')'` scaled by `1/d` (reciprocal 0 -> 1): DTRSM when `M` has a nonzero below the
/// diagonal, else only the scaling. Returns `(d, reciprocals)`; `M` (lower triangle, diagonal = d) and `R` are updated in place.
fn factor_block(
    mm: &mut [f64],
    r: &mut [f64],
    w: usize,
    c: usize,
    orig: &[f64],
    reg: Option<(usize, usize)>,
) -> (Vec<f64>, Vec<f64>) {
    let dense = (1..w).any(|i| (0..i).any(|j| mm[i * w + j] != 0.0));
    let (d, rcp) = if dense {
        dense_ldl(mm, w, orig, reg)
    } else {
        diag_only(mm, w, orig, reg)
    };
    let rows = c - w;
    if rows != 0 {
        if dense {
            for rr in 0..rows {
                let x = dtrsm_llu_col(mm, w, &r[rr * w..(rr + 1) * w], rr, rows);
                r[rr * w..(rr + 1) * w].copy_from_slice(&x);
            }
        }
        for rr in 0..rows {
            for k in 0..w {
                r[rr * w + k] = r[rr * w + k] * if rcp[k] != 0.0 { rcp[k] } else { 1.0 };
            }
        }
    }
    (d, rcp)
}

fn pack(mm: &[f64], r: &[f64], w: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity(w * (w + 1) / 2 + r.len());
    for i in 0..w {
        out.extend_from_slice(&mm[i * w..i * w + i + 1]);
    }
    out.extend_from_slice(r);
    out
}

// ---------------------------------------------------------------------------------------------------------------------------------
// Child -> parent update (0x2fcbd0).

struct Child<'a> {
    w: usize,
    c: usize,
    rows: &'a [usize],
    packed: &'a [f64],
    d: &'a [f64],
}

/// Child of width 1 (scalar branch of `0x2fcbd0`): `f_t = (-L_t) d`, every target gets exactly one product `L_u f_t`.
fn child_update_w1(
    ch: &Child,
    wp: usize,
    mp: &[usize],
    mpm: &mut [f64],
    rp: &mut [f64],
    cnt: usize,
    e: usize,
) {
    let (cc, l, dd) = (ch.c, ch.packed, ch.d[0]);
    let mpr: Vec<usize> = ch.rows.iter().map(|&r| mp[r]).collect();
    if wp < 5 {
        for t in cnt..e {
            let f = (-l[t]) * dd;
            for u in t..cc {
                if u < e {
                    let k = mpr[u] * wp + mpr[t];
                    mpm[k] = l[u] * f + mpm[k];
                } else {
                    let k = (mpr[u] - wp) * wp + mpr[t];
                    rp[k] = f * l[u] + rp[k];
                }
            }
        }
        return;
    }
    for t in cnt..cc {
        let f = (-l[t]) * dd;
        if t < e {
            for u in cnt..=t {
                let k = mpr[t] * wp + mpr[u];
                mpm[k] = l[u] * f + mpm[k];
            }
        } else {
            for u in cnt..e {
                let k = (mpr[t] - wp) * wp + mpr[u];
                rp[k] = f * l[u] + rp[k];
            }
        }
    }
}

/// `0x2fcbd0`, child `ch` into the front `(Mp, Rp)` of the parent (width `wp`, first column `col`). Rows of the child's list from `cnt` (first row inside
/// the parent's columns) to `e` (end of the parent's columns) are the `m` update rows, the rest the rect rows. `temp[t][k] = (-L[t][k]) d[k]`;
/// DGEMM('T','N', k = child width) blocks of 24 columns for the diagonal part (contiguous map: directly into `Mp`, else into a zero scratch and
/// scatter-add), then blocks of rect rows.
fn child_update(
    ch: &Child,
    wp: usize,
    col: usize,
    sym80: usize,
    mp: &[usize],
    mpm: &mut [f64],
    rp: &mut [f64],
) {
    let (wc, cc, rows_c, lp, d) = (ch.w, ch.c, ch.rows, ch.packed, ch.d);
    let cnt = (wc..cc)
        .find(|&t| rows_c[t] >= col)
        .expect("child without a row in the parent");
    let mut e = cnt;
    while e < cc && rows_c[e] < col + wp {
        e += 1;
    }
    if wc == 1 {
        child_update_w1(ch, wp, mp, mpm, rp, cnt, e);
        return;
    }
    let m = e - cnt;
    let tri = wc * (wc + 1) / 2;
    let lrow = |t: usize| &lp[tri + (t - wc) * wc..tri + (t - wc + 1) * wc];
    let mut temp = Vec::with_capacity(m * wc);
    for t in cnt..e {
        let lr = lrow(t);
        for k in 0..wc {
            temp.push((-lr[k]) * d[k]);
        }
    }
    let first = mp[rows_c[cnt]];
    let contiguous = mp[rows_c[e - 1]] == first + m - 1;
    let mut j0 = 0;
    while j0 < m {
        let nb = (m - j0).min(24);
        let mut b = Vec::with_capacity(nb * wc);
        for t in cnt + j0..cnt + j0 + nb {
            b.extend_from_slice(lrow(t));
        }
        if contiguous {
            let mmr = j0 + nb;
            let mut cm = Vec::with_capacity(mmr * nb);
            for jj in 0..nb {
                for ii in 0..mmr {
                    cm.push(mpm[(first + j0 + jj) * wp + first + ii]);
                }
            }
            dgemm_t_n(mmr, nb, wc, &temp, wc, &b, wc, &mut cm, mmr);
            for jj in 0..nb {
                for ii in 0..mmr {
                    mpm[(first + j0 + jj) * wp + first + ii] = cm[jj * mmr + ii];
                }
            }
        } else {
            let mut cn = vec![0.0; m * nb];
            dgemm_t_n(m, nb, wc, &temp, wc, &b, wc, &mut cn, m);
            for jj in 0..nb {
                let t = j0 + jj;
                for u in 0..=t {
                    let (r, cl) = (mp[rows_c[cnt + t]], mp[rows_c[cnt + u]]);
                    mpm[r * wp + cl] = mpm[r * wp + cl] + cn[jj * m + u];
                }
            }
        }
        j0 += nb;
    }
    let mut t0 = e;
    let blk = if m != 0 {
        (0xc0usize).min(sym80 / m).max(1)
    } else {
        1
    };
    while t0 < cc {
        let mut t1 = (t0 + blk).min(cc);
        let r_first = mp[rows_c[t0]];
        if contiguous && mp[rows_c[t1 - 1]] == r_first + (t1 - t0) - 1 {
            while t1 < cc && mp[rows_c[t1]] == r_first + (t1 - t0) {
                t1 += 1;
            }
            let nr = t1 - t0;
            let mut b = Vec::with_capacity(nr * wc);
            for t in t0..t1 {
                b.extend_from_slice(lrow(t));
            }
            let mut cm = Vec::with_capacity(m * nr);
            for jj in 0..nr {
                for ii in 0..m {
                    cm.push(rp[(r_first - wp + jj) * wp + first + ii]);
                }
            }
            dgemm_t_n(m, nr, wc, &temp, wc, &b, wc, &mut cm, m);
            for jj in 0..nr {
                for ii in 0..m {
                    rp[(r_first - wp + jj) * wp + first + ii] = cm[jj * m + ii];
                }
            }
        } else {
            let nr = t1 - t0;
            let mut b = Vec::with_capacity(nr * wc);
            for t in t0..t1 {
                b.extend_from_slice(lrow(t));
            }
            let mut cn = vec![0.0; m * nr];
            dgemm_t_n(m, nr, wc, &temp, wc, &b, wc, &mut cn, m);
            for jj in 0..nr {
                let r = mp[rows_c[t0 + jj]];
                for u in 0..m {
                    let cl = mp[rows_c[cnt + u]];
                    rp[(r - wp) * wp + cl] = rp[(r - wp) * wp + cl] + cn[jj * m + u];
                }
            }
        }
        t0 = t1;
    }
}

// ---------------------------------------------------------------------------------------------------------------------------------
// The numeric factor and the solves.

/// Supernodal factor as the solves read it: per supernode the packed `w (w + 1) / 2` diagonal block (row `i` holds `i + 1` entries, the strict lower part
/// `L[i][0..i)` is read) followed by the `(c - w) x w` rectangle row-major; `diag` (`SYM+0xe0`) and `rcp` (`SYM+0xe8`) per column.
#[derive(Debug, Clone)]
pub struct Factor {
    pub n: usize,
    pub vals: Vec<Vec<f64>>,
    pub diag: Vec<f64>,
    pub rcp: Vec<f64>,
    /// Parity of the element offset of each supernode's storage in the one contiguous, 16-byte aligned factor array (the solves' peels depend on it).
    pub par: Vec<usize>,
    /// `SYM+0x100`: the smallest `|pivot| / |assembled diagonal|` over the columns with a nonzero diagonal.
    pub minratio: f64,
}

/// `A D A' = L diag(d) L'` from the model and `D` (`n2 + m` entries: columns, then slacks): per supernode in column order the assembly, the children's
/// updates (ascending), the block factorization.
pub fn numeric_factor(
    f: &StdForm,
    rows: &[Vec<(usize, usize)>],
    plan: &Plan,
    d: &[f64],
    zx: Option<&[f64]>,
    mut on_sn: Option<&mut dyn FnMut(usize, &Supernode, &[f64], &[f64])>,
) -> Factor {
    let m = f.m;
    // Dense columns: the factor has `nt = m + nd` nodes (the dense nodes last), pivots regularise (`SYM+0x110`, `Q+0x218 > 0`), and every sparse supernode is
    // assembled by the generic routine; the dense supernode starts from the diagonal `-1/D[j]` of its columns, and a row's dense-node entries are the raw
    // coefficients `A[i][j]` (no `D`). The diagonal of a dense node is `-zx_j` (`z/x + v/u` of the iterate; the starting point has no `zx`: `-1/D_j`, as for a
    // clamped `D = 1e12`).
    let nt = plan.n;
    let reg = !f.dense.is_empty();
    let mut dnode_col = vec![NONE; nt];
    // Per dense column its coefficient in every row (0 where absent), indexed by the node `m + k`.
    let mut dval: Vec<Vec<f64>> = vec![Vec::new(); nt];
    for (r, &j) in f.dense.iter().enumerate() {
        dnode_col[m + f.dense_k[r]] = j;
        let mut v = vec![0.0; m];
        for q in f.col_start[j]..f.col_start[j + 1] {
            v[f.row_idx[q]] = f.val[q];
        }
        dval[m + f.dense_k[r]] = v;
    }
    let mut vals: Vec<Vec<f64>> = Vec::with_capacity(plan.sns.len());
    let mut dvec: Vec<Vec<f64>> = Vec::with_capacity(plan.sns.len());
    let mut diag = vec![0.0; nt];
    let mut rcp = vec![0.0; nt];
    let mut mp = vec![NONE; nt];
    let mut minratio = 1e300f64;
    for (k, s) in plan.sns.iter().enumerate() {
        let (col, w, c) = (s.col, s.w, s.c);
        for (t, &r) in s.rows.iter().enumerate() {
            mp[r] = t;
        }
        let (mut mm, mut r) = if reg && col >= m {
            let mut mm = vec![0.0; w * w];
            for t in 0..w {
                let jd = dnode_col[col + t];
                mm[t * w + t] = match zx {
                    Some(z) if z[jd] > 0.0 && d[jd] != 1e12 => -z[jd],
                    _ => -1.0 / d[jd],
                };
            }
            (mm, vec![0.0; (c - w) * w])
        } else if !reg && w >= 3 && c * w <= 30000 {
            assemble(f, rows, d, &mp, col, w, c)
        } else {
            assemble_generic(f, rows, d, &mp, col, w, c)
        };
        if reg && col < m {
            for (t, &rr) in s.rows.iter().enumerate() {
                if rr >= m && t >= w {
                    for a in 0..w {
                        r[(t - w) * w + a] = dval[rr][col + a];
                    }
                }
            }
        }
        let orig: Vec<f64> = (0..w).map(|i| mm[i * w + i]).collect();
        for &ci in &s.upd {
            let cs = &plan.sns[ci];
            let ch = Child {
                w: cs.w,
                c: cs.c,
                rows: &cs.rows,
                packed: &vals[ci],
                d: &dvec[ci],
            };
            child_update(&ch, w, col, plan.sym80, &mp, &mut mm, &mut r);
        }
        let (dd, rc) = factor_block(
            &mut mm,
            &mut r,
            w,
            c,
            &orig,
            if reg { Some((col, m)) } else { None },
        );
        let pk = pack(&mm, &r, w);
        for t in 0..w {
            diag[col + t] = dd[t];
            rcp[col + t] = rc[t];
            if orig[t] != 0.0 {
                minratio = minratio.min(dd[t].abs() / orig[t].abs());
            }
        }
        if let Some(cb) = on_sn.as_deref_mut() {
            cb(k, s, &pk, &dd);
        }
        vals.push(pk);
        dvec.push(dd);
        for &r_ in &s.rows {
            mp[r_] = NONE;
        }
    }
    let mut par = Vec::with_capacity(vals.len());
    let mut off = 0usize;
    for v in &vals {
        par.push(off & 1);
        off += v.len();
    }
    Factor {
        n: nt,
        vals,
        diag,
        rcp,
        par,
        minratio,
    }
}

/// Sum of the rounded products `term(i)`, `i < n`, as `0x2ec200` / `0x2ec480` add them, accumulator starting at 0: fewer than 8 terms (or fewer than
/// `peel + 8`): sequential; else `peel` (0, or 1 when the 16-byte alignment of the value pointer is 8) terms sequentially, blocks of 8 into eight
/// lanes (lane `l` of the block gets term `l`; lane 0 starts from the peeled sum), total `((l0 + l2) + (l4 + l6)) + ((l1 + l3) + (l5 + l7))`, then the scalar tail.
fn dot8(n: usize, peel: usize, term: &dyn Fn(usize) -> f64) -> f64 {
    if n < 8 || n < peel + 8 {
        let mut a = 0.0;
        for i in 0..n {
            a = a + term(i);
        }
        return a;
    }
    let mut a0 = 0.0;
    for i in 0..peel {
        a0 = a0 + term(i);
    }
    let end = n - ((n - peel) & 7);
    let mut l = [a0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let mut k = peel;
    while k < end {
        for q in 0..8 {
            l[q] = l[q] + term(k + q);
        }
        k += 8;
    }
    let mut tot = ((l[0] + l[2]) + (l[4] + l[6])) + ((l[1] + l[3]) + (l[5] + l[7]));
    for k in end..n {
        tot = tot + term(k);
    }
    tot
}

/// `x0 - sum_j a_j b_j` as the forward solve's single-row kernel (`0x307e50`): lane 0 starts from `x0`, the other lanes at 0 and subtract; peel and
/// blocks as [`dot8`], same pairing; the scalar tail subtracts.
fn dsub8(x0: f64, a: &[f64], b: &[f64], peel: usize) -> f64 {
    let n = a.len();
    if n < 8 || n < peel + 8 {
        let mut v = x0;
        for j in 0..n {
            v = v - a[j] * b[j];
        }
        return v;
    }
    let mut v = x0;
    for j in 0..peel {
        v = v - a[j] * b[j];
    }
    let end = n - ((n - peel) & 7);
    let mut l = [v, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let mut k = peel;
    while k < end {
        for q in 0..8 {
            l[q] = l[q] - a[k + q] * b[k + q];
        }
        k += 8;
    }
    let mut tot = ((l[0] + l[2]) + (l[4] + l[6])) + ((l[1] + l[3]) + (l[5] + l[7]));
    for j in end..n {
        tot = tot - a[j] * b[j];
    }
    tot
}

/// The 4-row group kernel of the forward solve: lane 0 (even terms) starts from `x0`, lane 1 (odd terms) from 0, both subtract; `peel` scalar terms
/// first; fewer than 2 (or `peel + 2`) terms: all scalar; total `lane0 + lane1`, then the scalar tail subtracts.
fn dsub2(x0: f64, a: &[f64], b: &[f64], peel: usize) -> f64 {
    let n = a.len();
    if n < 2 || n < peel + 2 {
        let mut v = x0;
        for j in 0..n {
            v = v - a[j] * b[j];
        }
        return v;
    }
    let mut v = x0;
    for j in 0..peel {
        v = v - a[j] * b[j];
    }
    let end = n - ((n - peel) & 1);
    let (mut l0, mut l1) = (v, 0.0);
    let mut k = peel;
    while k < end {
        l0 = l0 - a[k] * b[k];
        l1 = l1 - a[k + 1] * b[k + 1];
        k += 2;
    }
    let mut tot = l0 + l1;
    for j in end..n {
        tot = tot - a[j] * b[j];
    }
    tot
}

/// `0x3338e0` -> `0x307e50` over the whole range (one task): supernodes ascending; the diagonal block (`x_i = x_i - sum_{j<i} L_ij x_j`, [`dsub8`] with
/// the peel of the row's pointer parity `i (i + 1) / 2`; a column whose pivot and reciprocal are both 0 gets `x_i = 0`), then the rectangle rows (groups of
/// 4 consecutive rows with [`dsub2`], peel from the group's last row pointer; the remaining `(c - w) % 4` rows single [`dsub8`]).
fn forward(plan: &Plan, fac: &Factor, x: &mut [f64]) {
    for (k, s) in plan.sns.iter().enumerate() {
        let (col, w, c) = (s.col, s.w, s.c);
        let par = fac.par[k];
        let vals = &fac.vals[k];
        let mut xs: Vec<f64> = x[col..col + w].to_vec();
        for i in 0..w {
            let mut v = xs[i];
            if i > 0 {
                let off = i * (i + 1) / 2;
                v = dsub8(v, &vals[off..off + i], &xs[..i], (par + off) & 1);
            }
            if fac.diag[col + i] == 0.0 && fac.rcp[col + i] == 0.0 {
                v = 0.0;
            }
            xs[i] = v;
        }
        x[col..col + w].copy_from_slice(&xs);
        let base = w * (w + 1) / 2;
        let rr = c - w;
        let g = rr / 4;
        for gi in 0..g {
            let p0 = gi * 4;
            let peel = (par + base + (p0 + 3) * w) & 1;
            for r in p0..p0 + 4 {
                let row = s.rows[w + r];
                x[row] = dsub2(x[row], &vals[base + r * w..base + (r + 1) * w], &xs, peel);
            }
        }
        for r in g * 4..rr {
            let row = s.rows[w + r];
            x[row] = dsub8(
                x[row],
                &vals[base + r * w..base + (r + 1) * w],
                &xs,
                (par + base + r * w) & 1,
            );
        }
    }
}

/// `0x333bb0` -> `0x309270`: supernodes descending; rectangle rows from the last one in groups of 4
/// (`x_j -= ((L[p-3][j] x[p-3] + L[p-2][j] x[p-2]) + L[p-1][j] x[p-1]) + L[p][j] x[p]`), the remaining front rows one at a time (descending);
/// then the diagonal block, `i = w-1..1`: `x_j -= L_ij x_i` for `j < i`.
fn backward(plan: &Plan, fac: &Factor, x: &mut [f64]) {
    for (k, s) in plan.sns.iter().enumerate().rev() {
        let (col, w, c) = (s.col, s.w, s.c);
        let vals = &fac.vals[k];
        let base = w * (w + 1) / 2;
        let rr = c - w;
        let mut xs: Vec<f64> = x[col..col + w].to_vec();
        let mut p = rr as isize - 1;
        while p - 3 >= 0 {
            let pu = p as usize;
            let xr: [f64; 4] = std::array::from_fn(|q| x[s.rows[w + pu - 3 + q]]);
            let lr: [&[f64]; 4] =
                std::array::from_fn(|q| &vals[base + (pu - 3 + q) * w..base + (pu - 2 + q) * w]);
            for j in 0..w {
                let sm =
                    ((lr[0][j] * xr[0] + lr[1][j] * xr[1]) + lr[2][j] * xr[2]) + lr[3][j] * xr[3];
                xs[j] = xs[j] - sm;
            }
            p -= 4;
        }
        while p >= 0 {
            let pu = p as usize;
            let xr = x[s.rows[w + pu]];
            let lr = &vals[base + pu * w..base + (pu + 1) * w];
            for j in 0..w {
                xs[j] = xs[j] - lr[j] * xr;
            }
            p -= 1;
        }
        for i in (1..w).rev() {
            let xi = xs[i];
            let off = i * (i + 1) / 2;
            for j in 0..i {
                xs[j] = xs[j] - vals[off + j] * xi;
            }
        }
        x[col..col + w].copy_from_slice(&xs);
    }
}

/// The matrix as the solves see it (rows in the final supernodal order): the column storage of the [`StdForm`] and the row copy (`rstart`, `ridx`,
/// `rval`; a row lists its entries by ascending column). Value arrays are taken as 16-byte aligned at index 0 (the peel is the start parity).
pub struct Amat {
    m: usize,
    nc: usize,
    rstart: Vec<usize>,
    ridx: Vec<usize>,
    rval: Vec<f64>,
}

impl Amat {
    pub fn new(f: &StdForm) -> Amat {
        let rows = model_rows_all(f);
        let mut rstart = Vec::with_capacity(f.m + 1);
        let (mut ridx, mut rval) = (Vec::new(), Vec::new());
        for r in &rows {
            rstart.push(ridx.len());
            for &(j, q) in r {
                ridx.push(j);
                rval.push(f.val[q]);
            }
        }
        rstart.push(ridx.len());
        Amat {
            m: f.m,
            nc: f.n2,
            rstart,
            ridx,
            rval,
        }
    }

    /// `0x2ec200`: `out_i` = sum over the row ([`dot8`], peel = start parity) + `w[nc + i]` unless the row is '='.
    fn aw(&self, f: &StdForm, w: &[f64]) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.m);
        for i in 0..self.m {
            let (s, e) = (self.rstart[i], self.rstart[i + 1]);
            let mut a = if e > s {
                dot8(e - s, s & 1, &|q| self.rval[s + q] * w[self.ridx[s + q]])
            } else {
                0.0
            };
            if f.sense[i] != b'=' {
                a = a + w[self.nc + i];
            }
            out.push(a);
        }
        out
    }

    /// `0x2ec480` (unscaled): `out_j` = sum over the column ([`dot8`], peel = start parity), then `out[nc + i] = y_i` for rows other than '=', else 0.
    fn aty(&self, f: &StdForm, y: &[f64]) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.nc + self.m);
        for j in 0..self.nc {
            let (s, e) = (f.col_start[j], f.col_start[j + 1]);
            let sj = f.scale.as_ref().map_or(0.0, |v| v[j]);
            out.push(if sj != 0.0 {
                // scaled branch (`Q+0x648[j] != 0`): the common value replaces the entries, no alignment peel
                dot8(e - s, 0, &|q| y[f.row_idx[s + q]] * sj)
            } else if e > s {
                dot8(e - s, s & 1, &|q| f.val[s + q] * y[f.row_idx[s + q]])
            } else {
                0.0
            });
        }
        for i in 0..self.m {
            out.push(if f.sense[i] == b'=' { 0.0 } else { y[i] });
        }
        out
    }
}

/// `0x305f10` for the LP case: `w = rx D`, `t = A w` (`0x2ec200`), `rhs = ry + t`, forward solve, `x *= 1/d`, backward solve, `s = A' dy`
/// (`0x2ec480`), `dx = (s - rx) D`. Returns `(dx, dy)`.
fn solve305(
    a: &Amat,
    plan: &Plan,
    fac: &Factor,
    f: &StdForm,
    d: &[f64],
    rx: &[f64],
    ry: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let n = f.n();
    // Dense columns (`SYM+0x138 == 0`): `w_j = 0`; the right-hand side of node `m + k` is the raw `rx_j`; the nodes' solution is `dx_j`.
    let mut is_dense = vec![false; n];
    for &j in &f.dense {
        is_dense[j] = true;
    }
    let w: Vec<f64> = (0..n)
        .map(|j| if is_dense[j] { 0.0 } else { rx[j] * d[j] })
        .collect();
    let t = a.aw(f, &w);
    let mut x: Vec<f64> = (0..f.m).map(|i| ry[i] + t[i]).collect();
    x.resize(fac.n, 0.0);
    for (r, &j) in f.dense.iter().enumerate() {
        x[f.m + f.dense_k[r]] = rx[j];
    }
    forward(plan, fac, &mut x);
    for i in 0..fac.n {
        x[i] = x[i] * fac.rcp[i];
    }
    backward(plan, fac, &mut x);
    let s = a.aty(f, &x[..f.m]);
    let mut dx: Vec<f64> = (0..n)
        .map(|j| {
            if is_dense[j] {
                s[j]
            } else {
                (s[j] - rx[j]) * d[j]
            }
        })
        .collect();
    for (r, &j) in f.dense.iter().enumerate() {
        dx[j] = x[f.m + f.dense_k[r]];
    }
    x.truncate(f.m);
    (dx, x)
}

// ---------------------------------------------------------------------------------------------------------------------------------
// ntry of the first iteration.

/// `ntry` of the first iteration (`0x2e4f10` tail at `0x3040a0`): 5 when `rho = (4 nnz(A) + 2 NZ) / S <= 0.2`, else 3, with `NZ = SYM+0x38`,
/// `S = (e/5 + 8 NZ + 2 a (1 - b)) + ((OPS - e) + a b) / 30 + 75 n`, `a` = sum of the squared column counts of the non-dense columns (`b = 0`; the dense columns stay in `nnz(A)`),
/// `e` = sum over supernodes of width < 5 of `F(c) - F(c - w)`, `F(x) = x (x + 1) (2 x + 1) / 6` (30 is the constant of the CPU path of the reference machine).
pub fn ntry_first(f: &StdForm, sym: &Symbolic) -> usize {
    if sym_rho_s(f, sym).0 <= 0.2 {
        5
    } else {
        3
    }
}

/// `(rho, S)` of [`ntry_first`]; `S = SYM+0x168` also scales the work of a factorization (the dense columns are nodes after the rows: `n = m + dense`).
pub fn sym_rho_s(f: &StdForm, sym: &Symbolic) -> (f64, f64) {
    let ff = |t: i64| t * (t + 1) * (2 * t + 1) / 6;
    let a = (0..f.n2)
        .filter(|j| !f.dense.contains(j))
        .map(|j| ((f.col_start[j + 1] - f.col_start[j]) as i64).pow(2))
        .sum::<i64>() as f64;
    let e = sym
        .width
        .iter()
        .zip(&sym.count)
        .filter(|(&w, _)| w > 0 && w < 5)
        .map(|(&w, &c)| ff(c) - ff(c - w))
        .sum::<i64>() as f64;
    let (nz, ops, n) = (sym.nz as f64, sym.ops, (f.m + f.dense.len()) as f64);
    let mut t = (e / 5.0 + 8.0 * nz) + 2.0 * a * (1.0 - 0.0);
    t = ((ops - e) + a * 0.0) / 30.0 + t;
    let s = (1.0 * t + 0.0) + n * 75.0;
    let nnz = f.row_idx.len() as f64;
    ((nnz * 2.0 * 2.0 + 2.0 * nz) / s, s)
}

// ---------------------------------------------------------------------------------------------------------------------------------
// The NormalSolver.

/// The bit-exact solver: supernodal factor of `A D A'` and the solves of the reference.
pub struct CholSolver {
    plan: Plan,
    rows: Vec<Vec<(usize, usize)>>,
    a: Amat,
    d: Vec<f64>,
    fac: Option<Factor>,
    ntry0: usize,
    nfac: usize,
    /// `Q+0x800` (`z/x`, plus `v/u`) of the current iterate and the refinement threshold `Q+0x280`, set by [`NormalSolver::prepare`].
    zx: Vec<f64>,
    thr: f64,
    /// `zx` of the factor being computed (iteration factors only).
    fzx: Option<Vec<f64>>,
    /// `Q+0x770 / 0x780` of the last direction: the norms of its residual.
    norms: (f64, f64),
    /// `Q+0x790 / 0x798`: the norms of the residual of the homogeneous algorithm's second solve.
    raw_norms: (f64, f64),
    /// `NZ` and `S` of the symbolic factor (work model of the release profile) and the correction solves of the refinement since the last query.
    nz_s: (f64, f64),
    extra: std::cell::Cell<usize>,
}

impl CholSolver {
    /// `sym` is the symbolic factorization of the structure stage ([`crate::order::structure`]) the `StdForm` was permuted with.
    pub fn new(f: &StdForm, sym: &Symbolic) -> CholSolver {
        CholSolver {
            plan: Plan::new(sym),
            rows: model_rows(f),
            a: Amat::new(f),
            d: Vec::new(),
            fac: None,
            ntry0: ntry_first(f, sym),
            nfac: 0,
            zx: Vec::new(),
            thr: 0.0,
            fzx: None,
            norms: (0.0, 0.0),
            raw_norms: (0.0, 0.0),
            nz_s: (sym.nz as f64, sym_rho_s(f, sym).1),
            extra: std::cell::Cell::new(0),
        }
    }

    fn solve(&self, f: &StdForm, rx: &[f64], ry: &[f64]) -> (Vec<f64>, Vec<f64>) {
        solve305(
            &self.a,
            &self.plan,
            self.fac.as_ref().expect("factor before solve"),
            f,
            &self.d,
            rx,
            ry,
        )
    }

    /// `SYM+0x178 == 1` (`0x2f1720`): free columns remain, or the smallest pivot ratio is below `1e-10` (dense columns and scaling are not modelled).
    fn refining(&self, f: &StdForm) -> bool {
        (0..f.n2).any(|j| f.flags[j] == 0)
            || !f.dense.is_empty()
            || self.fac.as_ref().map_or(false, |x| x.minratio < 1e-10)
    }

    /// `0x311d70`, the residual of a solve `(dx, dy)` for `(rx, ry)`: `rows = ry - Af dx` and `cols = rx - (Af' dy - zx dx)` in blocks of 8 from index 0,
    /// `(rx - cols) + zx dx` in the scalar tail (`n mod 8` entries).
    fn residual(
        &self,
        f: &StdForm,
        rx: &[f64],
        ry: &[f64],
        dx: &[f64],
        dy: &[f64],
    ) -> (Vec<f64>, Vec<f64>) {
        let (m, n2, n) = (f.m, f.n2, f.n());
        let mut rows = vec![0.0; m];
        for j in 0..n2 {
            let sj = f.scale.as_ref().map_or(0.0, |v| v[j]);
            let t = dx[j] * sj;
            for k in f.col_start[j]..f.col_start[j + 1] {
                let i = f.row_idx[k];
                // scaled variant of `0x2ecab0`: `dx_j s_j` is added to every row of the column
                rows[i] = if sj != 0.0 {
                    rows[i] + t
                } else {
                    f.val[k] * dx[j] + rows[i]
                };
            }
        }
        let mut cols = vec![0.0; n];
        for j in 0..n2 {
            cols[j] = f.aty_col(j, dy);
        }
        for i in 0..m {
            if f.sense[i] != b'=' {
                rows[i] = rows[i] + dx[n2 + i];
                cols[n2 + i] = dy[i];
            }
        }
        let rows: Vec<f64> = (0..m).map(|i| ry[i] - rows[i]).collect();
        let end = if n >= 8 { n - (n & 7) } else { 0 };
        let cols = (0..n)
            .map(|j| {
                if j < end {
                    rx[j] - (cols[j] - self.zx[j] * dx[j])
                } else {
                    (rx[j] - cols[j]) + self.zx[j] * dx[j]
                }
            })
            .collect();
        (rows, cols)
    }

    /// The refinement loop of `0x311d70`: while `thr^2 <= |rows|^2 + |cols|^2` and fewer than 5 accepted corrections, solve for the residual and accept
    /// the corrected pair if `0.99 cur > cand`. Returns `(dx, dy, |rows|, |cols|)`.
    fn refine(
        &self,
        f: &StdForm,
        rx: &[f64],
        ry: &[f64],
        mut dx: Vec<f64>,
        mut dy: Vec<f64>,
    ) -> (Vec<f64>, Vec<f64>, f64, f64) {
        let norm = |v: &[f64]| crate::barrier::dot_lanes(v, v, 0).sqrt();
        let (mut rows, mut cols) = self.residual(f, rx, ry, &dx, &dy);
        let (mut nr, mut nc) = (norm(&rows), norm(&cols));
        let mut cur = nc * nc + nr * nr;
        let thr2 = self.thr * self.thr;
        let mut k = 0;
        while thr2 <= cur {
            let (cx, cy) = self.solve(f, &cols, &rows);
            self.extra.set(self.extra.get() + 1);
            let cdx: Vec<f64> = cx.iter().zip(&dx).map(|(a, b)| a + b).collect();
            let cdy: Vec<f64> = cy.iter().zip(&dy).map(|(a, b)| a + b).collect();
            let (r2, c2) = self.residual(f, rx, ry, &cdx, &cdy);
            let (nr2, nc2) = (norm(&r2), norm(&c2));
            let cand = nr2 * nr2 + nc2 * nc2;
            if !(0.99 * cur > cand) {
                break;
            }
            dx = cdx;
            dy = cdy;
            (rows, cols) = (r2, c2);
            nr = nr2;
            nc = nc2;
            cur = cand;
            k += 1;
            if k >= 5 {
                break;
            }
        }
        (dx, dy, nr, nc)
    }
}

impl NormalSolver for CholSolver {
    fn factor_iter(&mut self, f: &StdForm, d: &[f64], zx: &[f64]) -> Result<(), BarrierError> {
        self.fzx = Some(zx.to_vec());
        let r = self.factor(f, d);
        self.fzx = None;
        r
    }

    fn factor(&mut self, f: &StdForm, d: &[f64]) -> Result<(), BarrierError> {
        self.d = d.to_vec();
        let fno = self.nfac;
        self.nfac += 1;
        let on = dump::on("barrier_factor");
        let mut cb = |_k: usize, s: &Supernode, pk: &[f64], piv: &[f64]| {
            dump::emit_tag(
                "barrier_factor",
                "FSN",
                json!({"f": fno, "col": s.col, "w": s.w, "c": s.c, "rows": s.rows}),
            );
            dump::emit_tag(
                "barrier_factor",
                "FOUT",
                json!({"f": fno, "col": s.col, "packed": pk, "piv": piv}),
            );
        };
        let fac = if on {
            numeric_factor(
                f,
                &self.rows,
                &self.plan,
                d,
                self.fzx.as_deref(),
                Some(&mut cb),
            )
        } else {
            numeric_factor(f, &self.rows, &self.plan, d, self.fzx.as_deref(), None)
        };
        self.fac = Some(fac);
        Ok(())
    }

    /// Primal (`0x2e8ef0` stage 1): `rx = -ub` on the boxed columns, `ry = b`; returns `(dx, dy)` (the starting `x~`). Dual (stage 2): `rx = c`, `ry = 0`;
    /// returns `(dx, dy)` with `dx = (A'y - c) D` over all variables (the reference's `Q+0x7a8`).
    fn start(
        &mut self,
        f: &StdForm,
        kind: StartKind,
    ) -> Result<(Vec<f64>, Vec<f64>), BarrierError> {
        let n = f.n();
        match kind {
            StartKind::Primal { .. } => {
                let mut rx = vec![0.0; n];
                for j in 0..f.n2 {
                    if f.flags[j] & 2 != 0 {
                        rx[j] = -f.ub[j];
                    }
                }
                Ok(self.solve(f, &rx, &f.b))
            }
            StartKind::Dual => {
                let mut rx = f.c.clone();
                rx.resize(n, 0.0);
                Ok(self.solve(f, &rx, &vec![0.0; f.m]))
            }
        }
    }

    /// `0x30fc90`: `rx = rd - r3 / x` (flag bit 0), `+ (r7d0 - ru v) / u` (bit 1); the solve; `dz = (r3 - z dx) / x`; `du = (ru - dx) + 0 ub`,
    /// `dv = (r7d0 - v du) / u`. No refinement.
    fn direction(
        &mut self,
        f: &StdForm,
        s: &State,
        _kind: DirKind,
        rhs: &Rhs,
    ) -> Result<Direction, BarrierError> {
        let (n, n2) = (f.n(), f.n2);
        let mut rx = rhs.rd.clone();
        for j in 0..n {
            if f.flags[j] & 1 != 0 {
                rx[j] = rx[j] - rhs.rcx[j] / s.x[j];
            }
        }
        for j in f.ub_lo..f.ub_hi {
            if f.flags[j] & 2 != 0 {
                rx[j] = (rhs.rcu[j] - rhs.ru[j] * s.v[j]) / s.u[j] + rx[j];
            }
        }
        let (mut dx, mut dy) = self.solve(f, &rx, &rhs.rp);
        self.norms = (0.0, 0.0);
        if self.refining(f) {
            let (a, b, nr, nc) = self.refine(f, &rx, &rhs.rp, dx, dy);
            dx = a;
            dy = b;
            self.norms = (nr, nc);
        }
        let mut dz = vec![0.0; n];
        for j in 0..n {
            if f.flags[j] & 1 != 0 {
                dz[j] = (rhs.rcx[j] - s.z[j] * dx[j]) / s.x[j];
            }
        }
        let (mut du, mut dv) = (vec![0.0; n2], vec![0.0; n2]);
        for j in f.ub_lo..f.ub_hi {
            if f.flags[j] & 2 != 0 {
                let t = (rhs.ru[j] - dx[j]) + 0.0 * f.ub[j];
                du[j] = t;
                dv[j] = (rhs.rcu[j] - s.v[j] * t) / s.u[j];
            }
        }
        Ok(Direction { dx, dy, dz, du, dv })
    }

    fn prepare(&mut self, zx: &[f64], thr: f64) {
        self.zx = zx.to_vec();
        self.thr = thr;
    }

    fn norms(&self) -> (f64, f64) {
        self.norms
    }

    fn work_info(&self, _f: &StdForm) -> Option<(f64, f64)> {
        Some(self.nz_s)
    }

    fn refine_mode(&self, f: &StdForm) -> bool {
        self.refining(f)
    }

    fn take_extra_solves(&mut self) -> usize {
        self.extra.replace(0)
    }

    fn ntry(&mut self, it: i32, default: usize) -> usize {
        if it == 0 {
            self.ntry0
        } else {
            default
        }
    }

    fn solve_raw(
        &mut self,
        f: &StdForm,
        rx: &[f64],
        ry: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>), BarrierError> {
        // `0x3117e0`: the same refinement as the direction solves (`0x311d70`), norms to `Q+0x790 / 0x798`; its correction solves are not charged.
        let (mut dx, mut dy) = self.solve(f, rx, ry);
        self.raw_norms = (0.0, 0.0);
        if self.refining(f) {
            let extra = self.extra.get();
            let (a, b, nr, nc) = self.refine(f, rx, ry, dx, dy);
            self.extra.set(extra);
            dx = a;
            dy = b;
            self.raw_norms = (nr, nc);
        }
        Ok((dx, dy))
    }

    fn raw_norms(&self) -> (f64, f64) {
        self.raw_norms
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;

    fn floats(v: &Value) -> Vec<f64> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_f64().unwrap())
            .collect()
    }

    fn assert_bits(what: &str, want: &[f64], got: &[f64]) {
        assert_eq!(want.len(), got.len(), "{what}: length");
        for (i, (w, g)) in want.iter().zip(got).enumerate() {
            assert_eq!(
                w.to_bits(),
                g.to_bits(),
                "{what}[{i}]: model {w:e} rose {g:e}"
            );
        }
    }

    // Inputs and outputs of the Python model `tools/barrier_chol.py` (bit-exact against the reference's kernels), written by
    // `reverse-engineer/tools/make_chol_kernel_fixture.py`.
    fn fixture() -> Value {
        serde_json::from_str(include_str!("../tests/fixtures/chol-kernels.json")).unwrap()
    }

    #[test]
    fn dgemm_kernels_match_the_model() {
        for g in fixture()["gemm"].as_array().unwrap() {
            let u = |k: &str| g[k].as_u64().unwrap() as usize;
            let (m, n, k) = (u("m"), u("n"), u("k"));
            let mut c = floats(&g["C"]);
            dgemm_t_n(
                m,
                n,
                k,
                &floats(&g["A"]),
                u("lda"),
                &floats(&g["B"]),
                u("ldb"),
                &mut c,
                u("ldc"),
            );
            assert_bits(&format!("dgemm {m}x{n}x{k}"), &floats(&g["out"]), &c);
        }
    }

    #[test]
    fn dtrsm_single_column_matches_the_reference_dtrsv() {
        // Reference outputs of DTRSM with n = 1 (probe in the live library): blocks of 32 rows, pairs of rows, 16-term fused blocks.
        let v: Value =
            serde_json::from_str(include_str!("../tests/fixtures/dtrsv-probe.json")).unwrap();
        for c in v["cases"].as_array().unwrap() {
            let w = c["w"].as_u64().unwrap() as usize;
            assert_bits(
                &format!("dtrsv w {w}"),
                &floats(&c["x"]),
                &dtrsm_llu_col(&floats(&c["L"]), w, &floats(&c["b"]), 0, 1),
            );
        }
        // E226 factor 2, supernode 181 (m 2, n 1): the last row term is fused.
        let f = f64::from_bits;
        let l = [1.0, 0.0, f(0xbfe92d2a474b121c), 1.0];
        let b = [f(0xbfdb6170a82dfb5e), f(0x3fdb6170a82dfb5e)];
        assert_bits(
            "dtrsv e226",
            &[f(0xbfdb6170a82dfb5e), f(0x3fb75aadc5e50f0b)],
            &dtrsm_llu_col(&l, 2, &b, 0, 1),
        );
    }

    #[test]
    fn dtrsm_kernels_match_the_model() {
        for t in fixture()["trsm"].as_array().unwrap() {
            let (w, n) = (
                t["w"].as_u64().unwrap() as usize,
                t["n"].as_u64().unwrap() as usize,
            );
            if n == 1 && w > 1 {
                continue; // the fixture's model has no dtrsv path; see dtrsm_single_column_matches_the_reference_dtrsv
            }
            let l = floats(&t["L"]);
            for c in t["cols"].as_array().unwrap() {
                let col = c["col"].as_u64().unwrap() as usize;
                assert_bits(
                    &format!("dtrsm w {w} n {n} col {col}"),
                    &floats(&c["x"]),
                    &dtrsm_llu_col(&l, w, &floats(&c["b"]), col, n),
                );
            }
        }
    }

    #[test]
    fn dense_ldl_matches_the_model() {
        for t in fixture()["ldl"].as_array().unwrap() {
            let w = t["w"].as_u64().unwrap() as usize;
            let mut mm = floats(&t["M"]);
            let (d, rcp) = dense_ldl(&mut mm, w, &floats(&t["orig"]), None);
            let lower: Vec<f64> = (0..w)
                .flat_map(|i| (0..=i).map(move |j| (i, j)))
                .map(|(i, j)| mm[i * w + j])
                .collect();
            assert_bits(&format!("ldl w {w}"), &floats(&t["out"]), &lower);
            assert_bits("d", &floats(&t["d"]), &d);
            assert_bits("rcp", &floats(&t["rcp"]), &rcp);
        }
    }

    /// The factor reproduces `A D A'` (relative error of the order of the rounding) on a small model: a guard against a wrong assembly or update list that the
    /// bit-exact tests would also catch, but with a readable message.
    #[test]
    fn factor_reproduces_the_normal_matrix() {
        let model = crate::mps::parse(include_str!("../tests/fixtures/adlittle.mps")).unwrap();
        let s = crate::order::structure(&model, crate::order::DENSE_COL_MIN).unwrap();
        let f = crate::barrier::std_form(&model, &s.sym.perm).unwrap();
        let plan = Plan::new(&s.sym);
        let d: Vec<f64> = (0..f.n())
            .map(|j| {
                if f.flags[j] & 1 != 0 {
                    1.0 + 0.01 * j as f64
                } else {
                    0.0
                }
            })
            .collect();
        let fac = numeric_factor(&f, &model_rows(&f), &plan, &d, None, None);
        let m = f.m;
        // Dense A D A'.
        let mut a = vec![0.0; m * m];
        for j in 0..f.n2 {
            for p in f.col_start[j]..f.col_start[j + 1] {
                for q in f.col_start[j]..f.col_start[j + 1] {
                    a[f.row_idx[p] * m + f.row_idx[q]] += f.val[p] * d[j] * f.val[q];
                }
            }
        }
        for i in 0..m {
            if f.sense[i] != b'=' {
                a[i * m + i] += d[f.n2 + i];
            }
        }
        // L diag(d) L' from the packed supernodes.
        let mut l = vec![0.0; m * m];
        for (k, s) in plan.sns.iter().enumerate() {
            let (w, c, vals) = (s.w, s.c, &fac.vals[k]);
            for t in 0..w {
                l[(s.col + t) * m + s.col + t] = 1.0;
                for r in t + 1..c {
                    let v = if r < w {
                        vals[r * (r + 1) / 2 + t]
                    } else {
                        vals[w * (w + 1) / 2 + (r - w) * w + t]
                    };
                    l[s.rows[r] * m + s.col + t] = v;
                }
            }
        }
        let mut worst = 0.0f64;
        for i in 0..m {
            for j in 0..=i {
                let mut sum = 0.0;
                for k in 0..=j {
                    sum += l[i * m + k] * fac.diag[k] * l[j * m + k];
                }
                worst = worst.max((a[i * m + j] - sum).abs() / (1.0 + a[i * m + j].abs()));
            }
        }
        assert!(worst < 1e-10, "A D A' - L D L' = {worst:e}");
    }
}
