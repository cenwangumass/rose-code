//! Compute forward rows directly in the reflected primal layout. Each tile
//! feeds the dual projection before the next tile replaces its row products.
use super::{dual_project, gather_sum, Scaled};
pub(super) struct Forward {
    ends: Vec<usize>,
    spans: Vec<Span>,
}
struct Span {
    input: Option<u32>,
    blocks: usize,
}
impl Forward {
    pub(super) fn new(offsets: &[usize], index: &[u32]) -> Self {
        let mut result = Self {
            ends: vec![0],
            spans: Vec::new(),
        };
        for z in offsets.windows(2) {
            let (b, e) = (z[0], z[1]);
            let mark = result.spans.len();
            let mut linear = false;
            for block in index[(b + (b & 1)).min(e)..e].chunks_exact(8) {
                let start = block[0];
                let dense = block
                    .iter()
                    .enumerate()
                    .all(|(i, v)| start.checked_add(i as u32) == Some(*v));
                linear |= dense;
                if let Some(prev) = result.spans[mark..].last_mut() {
                    if (dense
                        && prev
                            .input
                            .is_some_and(|i| i as u64 + 8 * prev.blocks as u64 == start as u64))
                        || (!dense && prev.input.is_none())
                    {
                        prev.blocks += 1;
                        continue;
                    }
                }
                result.spans.push(Span {
                    input: dense.then_some(start),
                    blocks: 1,
                });
            }
            if !linear {
                result.spans.truncate(mark);
            }
            result.ends.push(result.spans.len());
        }
        result
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn project(
        &self,
        a: &Scaled,
        index: &[u32],
        x: &[f64],
        y: &mut [f64],
        anchor: &[f64],
        sigma: f64,
        w: f64,
    ) {
        assert_eq!(x.len(), a.n);
        assert_eq!(y.len(), a.m);
        assert_eq!(anchor.len(), a.m);
        let mut products = [0.; 32];
        for (tile, target) in y.chunks_mut(products.len()).enumerate() {
            let first = tile * products.len();
            let prod = &mut products[..target.len()];
            for (i, result) in prod.iter_mut().enumerate() {
                let r = first + i;
                let (b, e) = (a.rp[r], a.rp[r + 1]);
                let spans = &self.spans[self.ends[r]..self.ends[r + 1]];
                *result = if spans.is_empty() {
                    gather_sum(&a.rval, index, x, b, e)
                } else {
                    unsafe { Self::dot(spans, &a.rval, index, x, b, e) }
                };
            }
            dual_project(
                target,
                &anchor[first..],
                prod,
                &a.b[first..],
                &a.eq[first..],
                sigma,
                w,
            );
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    unsafe fn dot(
        spans: &[Span],
        value: &[f64],
        index: &[u32],
        x: &[f64],
        begin: usize,
        end: usize,
    ) -> f64 {
        use std::arch::x86_64::*;
        let mut k = begin;

        let mut seed = 0.0;
        if k & 1 != 0 {
            seed += value[k] * x[index[k] as usize];
            k += 1;
        }
        let mut lo = _mm256_set_pd(0.0, 0.0, 0.0, seed);
        let mut hi = _mm256_setzero_pd();
        for span in spans {
            let limit = k + 8 * span.blocks;
            // Construction consumed whole blocks from this row and verified
            // every input index. A dense span additionally proved start+i.
            unsafe {
                if let Some(start) = span.input {
                    let mut input = x.as_ptr().add(start as usize);
                    while k < limit {
                        let val = value.as_ptr().add(k);
                        lo = _mm256_add_pd(
                            lo,
                            _mm256_mul_pd(_mm256_loadu_pd(input), _mm256_loadu_pd(val)),
                        );
                        hi = _mm256_add_pd(
                            hi,
                            _mm256_mul_pd(
                                _mm256_loadu_pd(input.add(4)),
                                _mm256_loadu_pd(val.add(4)),
                            ),
                        );
                        k += 8;
                        input = input.add(8);
                    }
                } else {
                    while k < limit {
                        let ids = index.as_ptr().add(k);
                        let val = value.as_ptr().add(k);
                        let a = _mm256_i32gather_pd::<8>(x.as_ptr(), _mm_loadu_si128(ids.cast()));
                        let b = _mm256_i32gather_pd::<8>(
                            x.as_ptr(),
                            _mm_loadu_si128(ids.add(4).cast()),
                        );
                        lo = _mm256_add_pd(lo, _mm256_mul_pd(a, _mm256_loadu_pd(val)));
                        hi = _mm256_add_pd(hi, _mm256_mul_pd(b, _mm256_loadu_pd(val.add(4))));
                        k += 8;
                    }
                }
            }
        }
        let lower = _mm_add_pd(_mm256_castpd256_pd128(lo), _mm256_extractf128_pd::<1>(lo));
        let upper = _mm_add_pd(_mm256_castpd256_pd128(hi), _mm256_extractf128_pd::<1>(hi));
        let both = _mm_add_pd(lower, upper);
        let mut sum = _mm_cvtsd_f64(_mm_add_sd(both, _mm_unpackhi_pd(both, both)));
        for i in k..end {
            sum += value[i] * x[index[i] as usize];
        }

        if sum.is_nan() {
            gather_sum(value, index, x, begin, end)
        } else {
            sum
        }
    }
}
