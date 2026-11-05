//! The ordinary iteration keeps short columns in arity order. Four independent
//! scalar-order dots feed four matching primal projections without staging Aᵀy.
use super::{gather_sum, Scaled};

enum Input {
    Broadcast(u32),
    Linear(u32),
    Indexed([u32; 4]),
}
struct PacketTerm {
    input: Input,
    value: [f64; 4],
}
struct Batch {
    start: usize,
    count: usize,
    len: usize,
    terms: Vec<PacketTerm>,
}

pub(super) struct Ordered {
    columns: Vec<usize>,
    mapped_col: Vec<u32>,
    forward: super::forward::Forward,
    batches: Vec<Batch>,
    loose: usize,
    c: Vec<f64>,
    lb: Vec<f64>,
    ub: Vec<f64>,
    x: Vec<f64>,
    anchor: Vec<f64>,
    xb: Vec<f64>,
    live: bool,
}
// AVX2 interprets every gather offset as a signed i32 on both matrix sides.
fn fits_gather_offsets(rows: usize, columns: usize) -> bool {
    rows.max(columns) <= i32::MAX as usize
}

#[cfg(test)]
#[test]
fn both_input_dimensions_fit_signed_gather_offsets() {
    let max = i32::MAX as usize;
    assert!(fits_gather_offsets(max, max));
    assert!(!fits_gather_offsets(max + 1, 1));
    assert!(!fits_gather_offsets(1, max + 1));
}

impl Ordered {
    pub(super) fn new(a: &Scaled) -> Self {
        assert!(a.crow.iter().all(|&i| (i as usize) < a.m));
        let mut small: [Vec<usize>; 9] = Default::default();
        let mut other = Vec::new();
        for j in 0..a.n {
            let b = a.cp[j];
            let n = a.cp[j + 1] - b;
            if n < 8 + (b & 1) {
                small[n].push(j);
            } else {
                other.push(j);
            }
        }
        let mut columns = Vec::with_capacity(a.n);
        let mut batches = Vec::new();
        for (n, rows) in small.into_iter().enumerate() {
            let (packets, tail) = rows.as_chunks::<4>();
            other.extend_from_slice(tail);
            if packets.is_empty() {
                continue;
            }
            let start = columns.len();
            let mut terms = Vec::with_capacity(packets.len() * n);
            for group in packets {
                columns.extend_from_slice(group);
                for k in 0..n {
                    let indices = group.map(|r| a.crow[a.cp[r] + k]);
                    let f = indices[0];
                    let input = if indices == [f; 4] {
                        Input::Broadcast(f)
                    } else if f.checked_add(3).is_some() && indices == [f, f + 1, f + 2, f + 3] {
                        Input::Linear(f)
                    } else {
                        Input::Indexed(indices)
                    };
                    terms.push(PacketTerm {
                        input,
                        value: group.map(|r| a.cval[a.cp[r] + k]),
                    });
                }
            }
            batches.push(Batch {
                start,
                count: packets.len(),
                len: n,
                terms,
            });
        }
        let loose = columns.len();
        columns.extend(other);
        let mut inverse = vec![0; a.n];
        for (i, &j) in columns.iter().enumerate() {
            inverse[j] = i as u32;
        }
        let mapped_col: Vec<u32> = a.rcol.iter().map(|&j| inverse[j as usize]).collect();
        let forward = super::forward::Forward::new(&a.rp, &mapped_col);
        let permute = |v: &[f64]| columns.iter().map(|&i| v[i]).collect();
        Self {
            c: permute(&a.c),
            lb: permute(&a.lb),
            ub: permute(&a.ub),
            x: vec![0.; a.n],
            anchor: vec![0.; a.n],
            xb: vec![0.; a.n],
            mapped_col,
            forward,
            columns,
            batches,
            loose,
            live: false,
        }
    }
    #[cfg(test)]
    pub(super) fn snapshot(&self, x: &[f64], xb: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let mut p = x.to_vec();
        let mut r = xb.to_vec();
        if self.live {
            for (i, &j) in self.columns.iter().enumerate() {
                p[j] = self.x[i];
                r[j] = self.xb[i];
            }
        }
        (p, r)
    }
    pub(super) fn observe(&mut self, x: &mut [f64], xb: &mut [f64]) {
        if self.live {
            for (i, &j) in self.columns.iter().enumerate() {
                x[j] = self.x[i];
                xb[j] = self.xb[i];
            }
            self.live = false;
        }
    }
    fn scalar(&mut self, i: usize, a: &Scaled, y: &[f64], tau: f64, w: f64) {
        let j = self.columns[i];
        let v = gather_sum(&a.cval, &a.crow, y, a.cp[j], a.cp[j + 1]);
        let old = self.x[i];
        let t = old - (self.c[i] - v) * tau;
        let xp = t.max(self.lb[i]).min(self.ub[i]);
        let xb = (xp + xp) - old;
        self.xb[i] = xb;
        self.x[i] = xb * w + self.anchor[i] * (1.0 - w);
    }
    pub(super) fn step(
        &mut self,
        a: &Scaled,
        y: &mut [f64],
        y0: &[f64],
        x: &[f64],
        anchor: &[f64],
        tau: f64,
        sigma: f64,
        w: f64,
    ) -> bool {
        #[cfg(target_arch = "x86_64")]
        if std::is_x86_feature_detected!("avx2") && fits_gather_offsets(a.m, a.n) {
            assert_eq!(y.len(), a.m);
            if !self.live {
                for (i, &j) in self.columns.iter().enumerate() {
                    self.x[i] = x[j];
                    self.anchor[i] = anchor[j];
                }
                self.live = true;
            }
            // Every packed input and output belongs to the checked CSC dimensions.
            unsafe {
                self.packed(a, y, tau, w);
            }
            for i in self.loose..self.columns.len() {
                self.scalar(i, a, y, tau, w);
            }
            unsafe {
                self.forward
                    .project(a, &self.mapped_col, &self.xb, y, y0, sigma, w);
            }
            return true;
        }
        false
    }
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    unsafe fn packed(&mut self, a: &Scaled, y: &[f64], tau: f64, w: f64) {
        for id in 0..self.batches.len() {
            unsafe {
                match self.batches[id].len {
                    0 => self.group::<0>(id, a, y, tau, w),
                    1 => self.group::<1>(id, a, y, tau, w),
                    2 => self.group::<2>(id, a, y, tau, w),
                    3 => self.group::<3>(id, a, y, tau, w),
                    4 => self.group::<4>(id, a, y, tau, w),
                    5 => self.group::<5>(id, a, y, tau, w),
                    6 => self.group::<6>(id, a, y, tau, w),
                    7 => self.group::<7>(id, a, y, tau, w),
                    8 => self.group::<8>(id, a, y, tau, w),
                    _ => unreachable!(),
                }
            }
        }
    }
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    unsafe fn group<const N: usize>(&mut self, id: usize, a: &Scaled, y: &[f64], tau: f64, w: f64) {
        use std::arch::x86_64::*;
        let start = self.batches[id].start;
        let count = self.batches[id].count;
        let wt = _mm256_set1_pd(w);
        let wt1 = _mm256_set1_pd(1.0 - w);
        let tv = _mm256_set1_pd(tau);
        for p in 0..count {
            let i = start + 4 * p;
            let mut sum = _mm256_setzero_pd();
            unsafe {
                for k in 0..N {
                    let term = self.batches[id].terms.get_unchecked(p * N + k);
                    let yy = match &term.input {
                        Input::Broadcast(j) => _mm256_set1_pd(*y.get_unchecked(*j as usize)),
                        Input::Linear(j) => _mm256_loadu_pd(y.as_ptr().add(*j as usize)),
                        Input::Indexed(ids) => _mm256_i32gather_pd::<8>(
                            y.as_ptr(),
                            _mm_loadu_si128(ids.as_ptr().cast()),
                        ),
                    };
                    sum =
                        _mm256_add_pd(sum, _mm256_mul_pd(_mm256_loadu_pd(term.value.as_ptr()), yy));
                }
                let pos = _mm256_loadu_pd(self.x.as_ptr().add(i));
                let rhs = _mm256_loadu_pd(self.c.as_ptr().add(i));
                let t = _mm256_sub_pd(pos, _mm256_mul_pd(_mm256_sub_pd(rhs, sum), tv));
                let lower = _mm256_max_pd(t, _mm256_loadu_pd(self.lb.as_ptr().add(i)));
                let lower = _mm256_blendv_pd(lower, t, _mm256_cmp_pd::<_CMP_UNORD_Q>(lower, lower));
                let xp = _mm256_min_pd(lower, _mm256_loadu_pd(self.ub.as_ptr().add(i)));
                let xp = _mm256_blendv_pd(xp, lower, _mm256_cmp_pd::<_CMP_UNORD_Q>(xp, xp));
                let xb = _mm256_sub_pd(_mm256_add_pd(xp, xp), pos);
                let next = _mm256_add_pd(
                    _mm256_mul_pd(xb, wt),
                    _mm256_mul_pd(_mm256_loadu_pd(self.anchor.as_ptr().add(i)), wt1),
                );
                let exception = _mm256_or_pd(
                    _mm256_cmp_pd::<_CMP_UNORD_Q>(t, t),
                    _mm256_cmp_pd::<_CMP_UNORD_Q>(next, next),
                );
                if _mm256_movemask_pd(exception) != 0 {
                    for j in i..i + 4 {
                        self.scalar(j, a, y, tau, w);
                    }
                } else {
                    _mm256_storeu_pd(self.x.as_mut_ptr().add(i), next);
                    _mm256_storeu_pd(self.xb.as_mut_ptr().add(i), xb);
                }
            }
        }
    }
}
