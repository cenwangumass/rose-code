//! Solver-owned geometry and scratch for the complete child-update schedule.
//! Numerical arrays and the single/double rounding mode remain live.
use super::*;
use std::ops::Range;

#[derive(Default)]
pub(super) struct Scratch {
    temp: Vec<f64>,
    product: Vec<f64>,
}

struct ScalarGroup {
    source: usize,
    coefficients: Range<usize>,
    addresses: Range<usize>,
    base: usize,
    rectangle: bool,
    factor_left: bool,
}

struct OutputRow {
    base: usize,
    count: usize,
    product: usize,
}

enum Output {
    Gather {
        rectangle: bool,
        base: usize,
        stride: usize,
    },
    Add {
        rectangle: bool,
        rows: Vec<OutputRow>,
    },
}

// These views preserve the original logical C matrix, including computed cells
// that the triangular scatter subsequently discards.
struct GatherDestination<'a> {
    values: &'a mut [f64],
    base: usize,
    stride: usize,
    leading: usize,
    columns: usize,
}
impl gemm_destination::Destination for GatherDestination<'_> {
    fn shape(&self) -> (usize, usize) {
        (self.leading, self.columns)
    }
    unsafe fn load(&self, i: gemm_destination::Index) -> f64 {
        *self
            .values
            .get_unchecked(self.base + i.col * self.stride + i.row)
    }
    unsafe fn store(&mut self, i: gemm_destination::Index, value: f64) {
        *self
            .values
            .get_unchecked_mut(self.base + i.col * self.stride + i.row) = value;
    }
    fn bounded(&self) -> bool {
        (0..self.columns).all(|j| {
            self.values[self.base + j * self.stride..self.base + j * self.stride + self.leading]
                .iter()
                .copied()
                .all(gemm_destination::parent_bound)
        })
    }
}
struct ScatterDestination<'a> {
    values: &'a mut [f64],
    rows: &'a [OutputRow],
    addresses: &'a [usize],
}
impl gemm_destination::Destination for ScatterDestination<'_> {
    fn shape(&self) -> (usize, usize) {
        (self.addresses.len(), self.rows.len())
    }
    unsafe fn load(&self, _: gemm_destination::Index) -> f64 {
        0.0
    }
    unsafe fn store(&mut self, i: gemm_destination::Index, value: f64) {
        let row = self.rows.get_unchecked(i.col);
        if i.row < row.count {
            let dst = row.base + *self.addresses.get_unchecked(i.row);
            let output = self.values.get_unchecked_mut(dst);
            *output = *output + value;
        } else {
            // The reference computes these products even though scatter drops
            // them. Retain their arithmetic and floating-point status effects.
            std::hint::black_box(value);
        }
    }
    fn bounded(&self) -> bool {
        // The original zero-C GEMM never inspects parent values. Its finite
        // products meet the parent only in one final addition, so even a NaN
        // parent cannot compete with another NaN operand. Validate addresses,
        // without scanning or restricting the live parent values.
        let last = self.addresses.iter().copied().max().unwrap_or(0);
        self.rows.iter().all(|r| {
            r.count <= self.addresses.len()
                && (r.count == 0
                    || r.base
                        .checked_add(last)
                        .is_some_and(|i| i < self.values.len()))
        })
    }
}

struct Batch {
    b: Range<usize>,
    rows_double: usize,
    rows_single: usize,
    columns: usize,
    leading: usize,
    output: Output,
}

pub(super) struct ChildProgram {
    wc: usize,
    cnt: usize,
    end: usize,
    addresses: Vec<usize>,
    scalar: Vec<ScalarGroup>,
    batches: Vec<Batch>,
}

impl ChildProgram {
    fn new(parent: &Supernode, child: &Supernode, sym80: usize, mp: &[usize]) -> Option<Self> {
        let (wp, wc, cc) = (parent.w, child.w, child.c);
        if wc == 0 || wp == 0 || child.rows.len() != cc {
            return None;
        }
        let cnt = (wc..cc).find(|&t| child.rows[t] >= parent.col)?;
        let mut end = cnt;
        while end < cc && child.rows[end] < parent.col + wp {
            end += 1;
        }
        if end == cnt {
            return None;
        }
        let mut mapped = Vec::with_capacity(cc);
        for &row in &child.rows {
            mapped.push(*mp.get(row)?);
        }
        if mapped[cnt..].iter().any(|&v| v >= parent.c)
            || mapped[cnt..].windows(2).any(|w| w[0] >= w[1])
        {
            return None;
        }
        if mapped[cnt..end].iter().any(|&v| v >= wp) || mapped[end..].iter().any(|&v| v < wp) {
            return None;
        }
        let mut result = Self {
            wc,
            cnt,
            end,
            addresses: Vec::new(),
            scalar: Vec::new(),
            batches: Vec::new(),
        };
        if wc == 1 {
            if wp < 5 {
                result.addresses = (cnt..cc)
                    .map(|u| {
                        if u < end {
                            mapped[u] * wp
                        } else {
                            (mapped[u] - wp) * wp
                        }
                    })
                    .collect();
                for t in cnt..end {
                    result.scalar.push(ScalarGroup {
                        source: t,
                        coefficients: t..end,
                        addresses: t - cnt..end - cnt,
                        base: mapped[t],
                        rectangle: false,
                        factor_left: false,
                    });
                    result.scalar.push(ScalarGroup {
                        source: t,
                        coefficients: end..cc,
                        addresses: end - cnt..cc - cnt,
                        base: mapped[t],
                        rectangle: true,
                        factor_left: true,
                    });
                }
            } else {
                result.addresses = mapped[cnt..end].to_vec();
                for t in cnt..cc {
                    let (coefficients, rectangle) = if t < end {
                        (cnt..t + 1, false)
                    } else {
                        (cnt..end, true)
                    };
                    let base = if rectangle {
                        (mapped[t] - wp) * wp
                    } else {
                        mapped[t] * wp
                    };
                    let addresses = coefficients.start - cnt..coefficients.end - cnt;
                    result.scalar.push(ScalarGroup {
                        source: t,
                        coefficients,
                        addresses,
                        base,
                        rectangle,
                        factor_left: false,
                    });
                }
            }
            return Some(result);
        }
        result.addresses = mapped[cnt..end].to_vec();
        let m = end - cnt;
        let tri = wc * (wc + 1) / 2;
        let first = mapped[cnt];
        let contiguous = mapped[end - 1] == first + m - 1;
        let mut j0 = 0;
        while j0 < m {
            let nb = (m - j0).min(24);
            let (rows_double, rows_single, leading, output) = if contiguous {
                let nr = j0 + nb;
                (
                    nr,
                    nr,
                    nr,
                    Output::Gather {
                        rectangle: false,
                        base: (first + j0) * wp + first,
                        stride: wp,
                    },
                )
            } else {
                let rows = (0..nb)
                    .map(|jj| OutputRow {
                        base: mapped[cnt + j0 + jj] * wp,
                        count: j0 + jj + 1,
                        product: jj * m,
                    })
                    .collect();
                (
                    m,
                    j0 + nb,
                    m,
                    Output::Add {
                        rectangle: false,
                        rows,
                    },
                )
            };
            result.batches.push(Batch {
                b: tri + (cnt + j0 - wc) * wc..tri + (cnt + j0 + nb - wc) * wc,
                rows_double,
                rows_single,
                columns: nb,
                leading,
                output,
            });
            j0 += nb;
        }
        let blk = 192.min(sym80 / m).max(1);
        let mut t0 = end;
        while t0 < cc {
            let mut t1 = (t0 + blk).min(cc);
            let r_first = mapped[t0];
            let output = if contiguous && mapped[t1 - 1] == r_first + t1 - t0 - 1 {
                while t1 < cc && mapped[t1] == r_first + t1 - t0 {
                    t1 += 1;
                }
                Output::Gather {
                    rectangle: true,
                    base: (r_first - wp) * wp + first,
                    stride: wp,
                }
            } else {
                let rows = (0..t1 - t0)
                    .map(|jj| OutputRow {
                        base: (mapped[t0 + jj] - wp) * wp,
                        count: m,
                        product: jj * m,
                    })
                    .collect();
                Output::Add {
                    rectangle: true,
                    rows,
                }
            };
            result.batches.push(Batch {
                b: tri + (t0 - wc) * wc..tri + (t1 - wc) * wc,
                rows_double: m,
                rows_single: m,
                columns: t1 - t0,
                leading: m,
                output,
            });
            t0 = t1;
        }
        Some(result)
    }

    pub(super) fn apply(
        &self,
        child: &Child<'_>,
        mm: &mut [f64],
        rect: &mut [f64],
        scratch: &mut Scratch,
    ) {
        if self.wc == 1 {
            let d = child.d[0];
            let mut source = NONE;
            let mut f = 0.0;
            for group in &self.scalar {
                if source != group.source {
                    f = (-child.packed[group.source]) * d;
                    source = group.source;
                }
                let coefficients = &child.packed[group.coefficients.clone()];
                let out = if group.rectangle {
                    &mut *rect
                } else {
                    &mut *mm
                };
                let addresses = &self.addresses[group.addresses.clone()];
                if group.factor_left {
                    for (&offset, &v) in addresses.iter().zip(coefficients) {
                        let dst = group.base + offset;
                        out[dst] = rf(f * v + out[dst]);
                    }
                } else {
                    for (&offset, &v) in addresses.iter().zip(coefficients) {
                        let dst = group.base + offset;
                        out[dst] = rf(v * f + out[dst]);
                    }
                }
            }
            return;
        }
        let wc = self.wc;
        let m = self.end - self.cnt;
        let tri = wc * (wc + 1) / 2;
        scratch.temp.clear();
        scratch.temp.reserve(m * wc);
        for t in self.cnt..self.end {
            for k in 0..wc {
                scratch
                    .temp
                    .push(rf((-child.packed[tri + (t - wc) * wc + k]) * child.d[k]));
            }
        }
        let single = SINGLE.with(|s| s.get());
        for batch in &self.batches {
            if !single {
                let direct = match &batch.output {
                    Output::Gather {
                        rectangle,
                        base,
                        stride,
                    } => {
                        let values = if *rectangle { &mut *rect } else { &mut *mm };
                        let mut dst = GatherDestination {
                            values,
                            base: *base,
                            stride: *stride,
                            leading: batch.leading,
                            columns: batch.columns,
                        };
                        gemm_destination::try_update(
                            batch.rows_double,
                            batch.columns,
                            wc,
                            &scratch.temp,
                            wc,
                            &child.packed[batch.b.clone()],
                            wc,
                            &mut dst,
                            batch.leading,
                        )
                    }
                    Output::Add { rectangle, rows } => {
                        let values = if *rectangle { &mut *rect } else { &mut *mm };
                        let mut dst = ScatterDestination {
                            values,
                            rows,
                            addresses: &self.addresses,
                        };
                        gemm_destination::try_update(
                            batch.rows_double,
                            batch.columns,
                            wc,
                            &scratch.temp,
                            wc,
                            &child.packed[batch.b.clone()],
                            wc,
                            &mut dst,
                            batch.leading,
                        )
                    }
                };
                if direct {
                    continue;
                }
            }
            scratch.product.clear();
            let len = batch.leading * batch.columns;
            match &batch.output {
                Output::Gather {
                    rectangle,
                    base,
                    stride,
                } => {
                    scratch.product.reserve(len);
                    let input = if *rectangle { &*rect } else { &*mm };
                    for jj in 0..batch.columns {
                        let start = base + jj * stride;
                        scratch
                            .product
                            .extend_from_slice(&input[start..start + batch.leading]);
                    }
                }
                Output::Add { .. } => scratch.product.resize(len, 0.0),
            }
            // Exact logical slice lengths also govern native finite dispatch and
            // whole-buffer single-precision conversion, including unused padding.
            dgemm_t_n(
                if single {
                    batch.rows_single
                } else {
                    batch.rows_double
                },
                batch.columns,
                wc,
                &scratch.temp,
                wc,
                &child.packed[batch.b.clone()],
                wc,
                &mut scratch.product,
                batch.leading,
            );
            match &batch.output {
                Output::Gather {
                    rectangle,
                    base,
                    stride,
                } => {
                    let output = if *rectangle { &mut *rect } else { &mut *mm };
                    for jj in 0..batch.columns {
                        for ii in 0..batch.leading {
                            output[base + jj * stride + ii] =
                                rf(scratch.product[jj * batch.leading + ii]);
                        }
                    }
                }
                Output::Add { rectangle, rows } => {
                    let output = if *rectangle { &mut *rect } else { &mut *mm };
                    for row in rows {
                        let product = &scratch.product[row.product..row.product + row.count];
                        for (u, &offset) in self.addresses[..row.count].iter().enumerate() {
                            let dst = row.base + offset;
                            output[dst] = rf(output[dst] + product[u]);
                        }
                    }
                }
            }
        }
    }
}

pub(super) struct NumericProgram {
    fronts: Vec<Vec<Option<ChildProgram>>>,
}
impl NumericProgram {
    pub(super) fn new(plan: &Plan) -> Self {
        let mut mp = vec![NONE; plan.n];
        let mut fronts = Vec::with_capacity(plan.sns.len());
        for parent in &plan.sns {
            for (i, &row) in parent.rows.iter().enumerate() {
                mp[row] = i;
            }
            fronts.push(
                parent
                    .upd
                    .iter()
                    .map(|&ci| ChildProgram::new(parent, &plan.sns[ci], plan.sym80, &mp))
                    .collect(),
            );
            for &row in &parent.rows {
                mp[row] = NONE;
            }
        }
        Self { fronts }
    }
    pub(super) fn child(&self, parent: usize, index: usize) -> Option<&ChildProgram> {
        self.fronts[parent][index].as_ref()
    }
}

#[cfg(all(test, target_arch = "x86_64"))]
mod tests {
    use super::*;
    use std::hint::black_box;

    fn bits(v: &[f64]) -> Vec<u64> {
        v.iter().map(|v| v.to_bits()).collect()
    }
    fn csr() -> u32 {
        let mut x = 0;
        unsafe {
            std::arch::asm!("stmxcsr [{p}]", p = in(reg) &mut x, options(nostack,preserves_flags));
        }
        x
    }
    fn set_csr(x: u32) {
        unsafe {
            std::arch::asm!("ldmxcsr [{p}]", p = in(reg) &x, options(nostack,preserves_flags));
        }
    }
    struct Restore(u32);
    impl Drop for Restore {
        fn drop(&mut self) {
            set_csr(self.0);
        }
    }

    #[test]
    fn discarded_product_keeps_inexact_flag() {
        let _restore = Restore(csr());
        let _single = SingleGuard::new(false);
        if !std::arch::is_x86_feature_detected!("fma") {
            return;
        }
        let q = f64::from_bits(1.0f64.to_bits() + 1);
        let a = [0.0, 0.0, q, 0.0];
        let b = [q, 0.0, 0.0, 0.0];
        let rows = [
            OutputRow {
                base: 0,
                count: 1,
                product: 0,
            },
            OutputRow {
                base: 2,
                count: 2,
                product: 2,
            },
        ];
        for mode in 0..16 {
            let controls =
                0x1f80 | ((mode & 3) << 13) | (((mode >> 2) & 1) << 15) | (((mode >> 3) & 1) << 6);
            let mut product = [0.0; 4];
            set_csr(controls);
            black_box(dgemm_t_n)(2, 2, 2, &a, 2, &b, 2, &mut product, 2);
            let mut old = [0.0; 4];
            for r in &rows {
                for i in 0..r.count {
                    old[r.base + i] = old[r.base + i] + product[r.product + i];
                }
            }
            let flags = csr();
            assert_ne!(flags & 32, 0, "unused output must be inexact");
            let mut actual = [0.0; 4];
            let mut dst = ScatterDestination {
                values: &mut actual,
                rows: &rows,
                addresses: &[0, 1],
            };
            set_csr(controls);
            assert!(gemm_destination::try_update(
                2, 2, 2, &a, 2, &b, 2, &mut dst, 2
            ));
            assert_eq!(flags, csr(), "mode={mode}");
            assert_eq!(bits(&old), bits(&actual));
        }
    }

    #[test]
    fn scatter_preserves_unrestricted_parent_values() {
        let _restore = Restore(csr());
        let _single = SingleGuard::new(false);
        if !std::arch::is_x86_feature_detected!("fma") {
            return;
        }
        let hi = f64::from_bits((1473u64 << 52) | 0x000f_ffff_ffff_ffff);
        let lo = f64::from_bits((573u64 << 52) | 0x000a_51fc_3461_8673);
        let a = [hi, lo, -hi, lo];
        let b = [hi, lo, lo, -hi];
        let rows = [
            OutputRow {
                base: 0,
                count: 1,
                product: 0,
            },
            OutputRow {
                base: 2,
                count: 2,
                product: 2,
            },
        ];
        for parent in [
            0.0,
            -0.0,
            f64::from_bits(1),
            -f64::from_bits(1),
            f64::MAX,
            -f64::MAX,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::from_bits(0x7ff8_0000_0000_0123),
            f64::from_bits(0xfff8_0000_0000_0456),
            f64::from_bits(0x7ff0_0000_0000_0789),
        ] {
            for mode in 0..16 {
                let controls = 0x1f80
                    | ((mode & 3) << 13)
                    | (((mode >> 2) & 1) << 15)
                    | (((mode >> 3) & 1) << 6);
                let mut product = [0.0; 4];
                let mut old = [parent; 4];
                set_csr(controls);
                black_box(dgemm_t_n)(2, 2, 2, &a, 2, &b, 2, &mut product, 2);
                for r in &rows {
                    for i in 0..r.count {
                        old[r.base + i] = old[r.base + i] + product[r.product + i];
                    }
                }
                let flags = csr();
                let mut actual = [parent; 4];
                let mut dst = ScatterDestination {
                    values: &mut actual,
                    rows: &rows,
                    addresses: &[0, 1],
                };
                set_csr(controls);
                assert!(gemm_destination::try_update(
                    2, 2, 2, &a, 2, &b, 2, &mut dst, 2
                ));
                assert_eq!(
                    flags,
                    csr(),
                    "flags mode={mode} parent={:x}",
                    parent.to_bits()
                );
                assert_eq!(
                    bits(&old),
                    bits(&actual),
                    "mode={mode} parent={:x}",
                    parent.to_bits()
                );
            }
        }
    }

    #[test]
    fn direct_gather_checks_logical_cells_and_full_operands() {
        let _restore = Restore(csr());
        let _single = SingleGuard::new(false);
        if !std::arch::is_x86_feature_detected!("fma") {
            return;
        }
        set_csr(0x1f80);
        let nan = f64::from_bits(0x7ff8_0000_0000_0123);
        let mut values = [0.0, 0.0, nan, nan, 0.0, 0.0, nan, nan];
        let mut dst = GatherDestination {
            values: &mut values,
            base: 0,
            stride: 4,
            leading: 2,
            columns: 2,
        };
        assert!(gemm_destination::try_update(
            2,
            2,
            1,
            &[1.0, 2.0],
            1,
            &[3.0, 4.0],
            1,
            &mut dst,
            2
        ));
        assert_eq!(
            bits(&values),
            bits(&[3.0, 6.0, nan, nan, 4.0, 8.0, nan, nan])
        );
        for (a, c) in [
            (vec![1.0, 2.0, nan], vec![0.0; 4]),
            (vec![1.0, 2.0], vec![0.0, 0.0, nan, 0.0]),
        ] {
            let mut actual = c.clone();
            let mut dst = GatherDestination {
                values: &mut actual,
                base: 0,
                stride: 2,
                leading: 2,
                columns: 2,
            };
            assert!(!gemm_destination::try_update(
                2,
                2,
                1,
                &a,
                1,
                &[3.0, 4.0],
                1,
                &mut dst,
                2
            ));
            assert_eq!(bits(&c), bits(&actual));
        }
    }

    #[test]
    #[cfg(target_arch = "x86_64")]
    fn complete_child_program_matches_frozen_kernel() {
        let _restore = Restore(csr());
        let mut seed = 0x3498a7b1_u64;
        let mut random = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let mut scratch = Scratch::default();
        let mut checks = 0;
        for wc in [1, 2, 4, 9] {
            for m in [1, 2, 3, 24, 25] {
                for diagonal_contiguous in [false, true] {
                    for rectangle_kind in 0..3 {
                        let wp = if diagonal_contiguous {
                            m + 2
                        } else {
                            2 * m + 1
                        };
                        let col = 64;
                        let rr = 17;
                        let rows: Vec<_> = (col..col + wp)
                            .chain((0..rr * 2 + 4).map(|i| col + wp + i * 2))
                            .collect();
                        let parent = Supernode {
                            col,
                            w: wp,
                            c: rows.len(),
                            rows: rows.clone(),
                            upd: vec![],
                        };
                        let mut child_rows: Vec<_> = (0..wc).collect();
                        for i in 0..m {
                            child_rows.push(rows[if diagonal_contiguous { i } else { i * 2 }]);
                        }
                        for i in 0..rr {
                            let j = match rectangle_kind {
                                0 => i,
                                1 => 2 * i,
                                _ => {
                                    if i < 5 {
                                        i
                                    } else {
                                        i + 2
                                    }
                                }
                            };
                            child_rows.push(rows[wp + j]);
                        }
                        let child_node = Supernode {
                            col: 0,
                            w: wc,
                            c: child_rows.len(),
                            rows: child_rows,
                            upd: vec![],
                        };
                        let mut mp = vec![NONE; *rows.last().unwrap() + 1];
                        for (i, &row) in rows.iter().enumerate() {
                            mp[row] = i;
                        }
                        let program = ChildProgram::new(&parent, &child_node, 17, &mp).unwrap();
                        let len = wc * (wc + 1) / 2 + (child_node.c - wc) * wc;
                        for pass in 0..3 {
                            let mut value = || match pass {
                                0 => (random() % 128) as f64 * 0.125 - 5.0,
                                1 => f64::from_bits(random()),
                                _ => match random() % 8 {
                                    0 => -0.0,
                                    1 => f64::from_bits(1),
                                    2 => 1e200,
                                    3 => -1e-200,
                                    4 => f64::INFINITY,
                                    5 => f64::from_bits(0x7ff8_0000_0000_0123),
                                    6 => f64::from_bits(0x7ff8_0000_0000_0456),
                                    _ => 1.0,
                                },
                            };
                            let packed: Vec<_> = (0..len).map(|_| value()).collect();
                            let d: Vec<_> = (0..wc).map(|_| value()).collect();
                            let mm: Vec<_> = (0..wp * wp).map(|_| value()).collect();
                            let rect: Vec<_> = (0..(parent.c - wp) * wp).map(|_| value()).collect();
                            let child = Child {
                                w: wc,
                                c: child_node.c,
                                rows: &child_node.rows,
                                packed: &packed,
                                d: &d,
                            };
                            for mode in 0..if m <= 3 { 16 } else { 1 } {
                                for single in [true, false] {
                                    let _single = SingleGuard::new(single);
                                    set_csr(
                                        0x1f80
                                            | ((mode & 3) << 13)
                                            | (((mode >> 2) & 1) << 15)
                                            | (((mode >> 3) & 1) << 6),
                                    );
                                    let (mut old_mm, mut old_r) = (mm.clone(), rect.clone());
                                    let (mut new_mm, mut new_r) = (mm.clone(), rect.clone());
                                    scratch.temp.resize(4096, f64::NAN);
                                    scratch.product.resize(8192, -0.0);
                                    let before = csr() & !63;
                                    set_csr(before);
                                    black_box(child_oracle::child_update)(
                                        &child,
                                        wp,
                                        col,
                                        17,
                                        &mp,
                                        &mut old_mm,
                                        &mut old_r,
                                    );
                                    let old_flags = csr();
                                    set_csr(before);
                                    black_box(ChildProgram::apply)(
                                        &program,
                                        &child,
                                        &mut new_mm,
                                        &mut new_r,
                                        &mut scratch,
                                    );
                                    let new_flags = csr();
                                    assert_eq!(bits(&old_mm), bits(&new_mm), "M wc={wc} m={m} contiguous={diagonal_contiguous} rectangle={rectangle_kind} pass={pass} single={single} mode={mode}");
                                    assert_eq!(
                                        bits(&old_r),
                                        bits(&new_r),
                                        "R wc={wc} m={m} pass={pass} single={single} mode={mode}"
                                    );
                                    assert_eq!(old_flags, new_flags, "flags wc={wc} m={m} pass={pass} single={single} mode={mode}");
                                    checks += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        println!("child_program_checks={checks}");
    }

    fn factor_bits(f: &Factor) -> Vec<u64> {
        let mut out = vec![f.n as u64, f.minratio.to_bits()];
        out.extend(f.par.iter().map(|&v| v as u64));
        out.extend(bits(&f.diag));
        out.extend(bits(&f.rcp));
        for v in &f.vals {
            out.push(v.len() as u64);
            out.extend(bits(v));
        }
        out
    }

    #[test]
    fn complete_factor_program_matches_original_small_models() {
        let _restore = Restore(csr());
        for name in ["afiro", "sc50a", "blend", "adlittle", "bore3d"] {
            let model = crate::mps::parse(
                &std::fs::read_to_string(format!("tests/fixtures/{name}.mps")).unwrap(),
            )
            .unwrap();
            let structure = crate::order::structure(&model, crate::order::DENSE_COL_MIN).unwrap();
            let mut form = crate::barrier::std_form(&model, &structure.sym.perm).unwrap();
            let mut plan = Plan::new(&structure.sym);
            let rows = model_rows(&form);
            let program = NumericProgram::new(&plan);
            let mut scratch = Scratch::default();
            let original_vals = form.val.clone();
            for pass in 0..6 {
                plan.single = pass % 2 == 0;
                for (i, v) in form.val.iter_mut().enumerate() {
                    *v = original_vals[i] * if pass >= 4 { 0.5 } else { 1.0 };
                }
                let d: Vec<_> = (0..form.n())
                    .map(|i| 0.125 * (1 + (i + pass) % 16) as f64)
                    .collect();
                for mode in 0..16 {
                    let controls = 0x1f80
                        | ((mode & 3) << 13)
                        | (((mode >> 2) & 1) << 15)
                        | (((mode >> 3) & 1) << 6);
                    let mut old_callbacks = Vec::new();
                    let mut old_cb = |k: usize, _: &Supernode, p: &[f64], d: &[f64]| {
                        old_callbacks.push((k, bits(p), bits(d)));
                    };
                    set_csr(controls);
                    let old = factor_oracle::numeric_factor(
                        &form,
                        &rows,
                        &plan,
                        &d,
                        None,
                        Some(&mut old_cb),
                    );
                    let old_flags = csr();
                    let mut new_callbacks = Vec::new();
                    let mut new_cb = |k: usize, _: &Supernode, p: &[f64], d: &[f64]| {
                        new_callbacks.push((k, bits(p), bits(d)));
                    };
                    set_csr(controls);
                    let new = numeric_factor_with(
                        &form,
                        &rows,
                        &plan,
                        &d,
                        None,
                        Some(&mut new_cb),
                        Some((&program, &mut scratch)),
                    );
                    assert_eq!(
                        old_flags,
                        csr(),
                        "factor flags {name} pass={pass} mode={mode}"
                    );
                    assert_eq!(
                        factor_bits(&old),
                        factor_bits(&new),
                        "factor {name} pass={pass}"
                    );
                    assert_eq!(old_callbacks, new_callbacks, "callbacks {name} pass={pass}");
                }
            }
        }
    }
}
