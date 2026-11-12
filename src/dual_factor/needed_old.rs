//! The `ftran_needed` of the stack before F3a (tail loop per column with per-column bookkeeping, tail list built for the dense route,
//! full-length tail-entry count), frozen as the oracle of the differential tests in `needed_tests`.
use super::*;

impl<R: Real> DualFactor<R> {
    pub(super) fn ftran_needed_old(
        &self,
        rows: &[usize],
        vals: &[R],
        needed: &[usize],
        dense_in: bool,
    ) -> Vec<(usize, R)> {
        let m = self.m;
        let (lout, count_l, e_l) = self.lstage_num(rows, vals, dense_in);
        // A dense L-stage vector is the U stage's work vector itself (it implies `tau_dense`); only its counts are needed.
        let (spike, dense_x) = match lout {
            LOut::List(list) => (list, None),
            LOut::Dense(mut x) => {
                let nz = Self::normalize_count(&mut x);
                (Vec::new(), Some((x, nz)))
            }
        };
        self.tau_l.set(e_l);
        let tau_dense = rows.len() * self.dthr > self.m
            || dense_in
            || self.lev_dense_ft
            || (!self.etas.is_empty() && Self::count_dense(count_l, self.dthr, self.m));
        debug_assert!(dense_x.is_none() || tau_dense);
        let tau_n = if tau_dense { None } else { Some(spike.len()) };
        // The inline dense joint loop of two dense vectors (`0x128c56..0x129152`) tests the raw value before the quotient: a tail column whose
        // raw value is at most the tolerance is zeroed even when the quotient exceeds it.
        // The same holds for the joint routine `0x12c300` (column sparse, tau dense; PILOT.JA call 249): whenever tau does not take the single-vector U stage
        // (`ok(n2)` of `joint_u` fails), the tail columns test the raw value as well.
        let joint_dense = !tau_n.is_some_and(|n| {
            self.nuc.len() + n < self.i21 || self.dens * (n as f64) < 2.0 * m as f64
        });
        // The spike's entry count and its tail-entry count (read only where `m > h33^2 n`), before the tail loop changes `x`.
        let (spike_n, spike_t) = match &dense_x {
            Some((x, nz)) => (
                *nz,
                if self.tail_count_needed(*nz) { (0..m).filter(|&p| self.kind[p] == -2 && x[p] != R::ZERO).count() } else { 0 },
            ),
            None => (
                spike.len(),
                if self.tail_count_needed(spike.len()) { spike.iter().filter(|e| self.kind[e.0] == -2).count() } else { 0 },
            ),
        };
        let mut tail: Vec<usize> = Vec::new();
        let (mut x, mut touched) = match dense_x {
            Some((x, _)) => (x, Vec::new()),
            None => {
                let mut x = vec![R::ZERO; m];
                let mut touched = vec![false; m];
                for &(p, v) in &spike {
                    x[p] = v;
                    touched[p] = true;
                    if !joint_dense && self.kind[p] == -2 {
                        tail.push(p);
                    }
                }
                (x, touched)
            }
        };
        let keep_tiny = !tau_dense && !joint_dense && self.tau_u_dense_nt(spike_n, spike_t);
        if joint_dense || (!R::QUAD && keep_tiny) {
            // The ordinary dense tail loop (12403b..12405f) divides and stores
            // every slot before testing the quotient, including unreachable zero
            // slots. Their signed zeros are consumed by requested rows and F2.
            tail.clear();
            // Extra zero columns are zeroed by the joint dense kernel. Reached
            // columns keep the same descending slot order as the symbolic walk.
            tail.extend((self.i21..m).rev().map(|s| self.at[s]));
        } else {
            let mut k = 0;
            while k < tail.len() {
                let c = tail[k];
                k += 1;
                for &j in self.ab.a(c) {
                    let j = j as usize;
                    if !touched[j] {
                        touched[j] = true;
                        tail.push(j);
                    }
                }
            }
            self.sort_tail_desc(&mut tail);
        }
        #[cfg(any(test, feature = "factor-shadow"))]
        let mut tau_raw = vec![false; m];
        // The U stage `0x123690` solves the tail in one of three loops: the heap walk and the sparse-row loop zero a quotient of at most the tolerance, the dense
        // slot loop (`0x12403b`: `m <= h33^2 c` or `5 h33^2 t > m - i21`, see `tau_u_work`) stores the quotient and only skips its scatter (FIT2P call 14
        // keeps a tail value of 4e-16 that enters the row dots of the nucleus positions). A dense input (`tau_dense`) zeroes it (FIT2P call 49), and so does the
        // joint dense loop (`joint_dense`, fx_e2_r25 call 6: a tail quotient of -1.8e-14 is 0 in the reference).
        //
        // The work estimates of the U stage are booked here, where a tail column is finished (a scatter only hits lower slots, so its value
        // is final): every term below is an integer-valued count (`len`, `2 len + 4`), hence the sums are exact and order free.
        let tdense = self.tau_u_dense_nt(spike_n, spike_t);
        let (mut nz, mut abs_a, mut tp_a2, mut tp_only_a) = (0usize, 0usize, 0usize, 0usize);
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        let mut surv_nz = 0usize;
        {
            let f2 = self.f2.borrow();
            let x_col = &f2.x_col;
            for &c in &tail {
                let raw = x[c].abs() > R::c(TOL);
                #[cfg(any(test, feature = "factor-shadow"))]
                {
                    tau_raw[c] = raw;
                }
                let xc = x[c] / self.diag[c];
                let zeroed = joint_dense && !raw;
                let tiny = !zeroed && xc.abs() <= R::c(TOL);
                let fin = if zeroed {
                    R::ZERO
                } else if tiny {
                    if keep_tiny {
                        xc
                    } else {
                        R::ZERO
                    }
                } else {
                    xc
                };
                x[c] = fin;
                if fin != R::ZERO {
                    #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
                    {
                        surv_nz += 1;
                    }
                    if tdense {
                        if fin.abs() > R::c(TOL) {
                            abs_a += self.ab.na(c);
                        }
                    } else {
                        nz += 1;
                    }
                    if joint_dense && raw {
                        let len = self.ab.na(c);
                        tp_a2 += 2 * len + 4;
                        if x_col.get(c).map_or(true, |v| *v == R::ZERO) {
                            tp_only_a += len;
                        }
                    }
                }
                if zeroed || tiny {
                    continue;
                }
                // The tail targets of the column are its A list (distinct rows, so the order of the scatter is free).
                for (&j, &v) in self.ab.a(c).iter().zip(self.ab.a_val(c)) {
                    x[j as usize] -= v * xc;
                }
            }
        }
        // Invariant (P4-2): every nonzero tail slot was reached through the a-lists and processed, so the nonzero tail values left are exactly the
        // survivors counted in the loop.
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        assert_eq!(
            (self.i21..m).filter(|&s| x[self.at[s]] != R::ZERO).count(),
            surv_nz,
            "ftran_needed: nonzero tail slots differ from the processed survivors"
        );
        let (u_est, fin) = self.tau_u_finish(spike_n, spike_t, needed.len(), nz, abs_a);
        self.tau_u.set(u_est);
        self.tau_fin.set(fin);
        #[cfg(test)]
        {
            let (v, f) = self.tau_u_work_ref(spike_n, spike_t, &x, needed.len());
            assert_eq!(v.to_bits(), u_est.to_bits(), "tau_u");
            assert_eq!(f, fin, "tau_fin");
        }
        {
            let mut f = self.f2.borrow_mut();
            f.tau_n = tau_n;
            #[cfg(any(test, feature = "factor-shadow"))]
            {
                f.x_tau = x.clone();
                f.tau_raw = tau_raw;
            }
            f.tp_a2 = tp_a2;
            f.tp_only_a = tp_only_a;
        }
        needed
            .iter()
            .map(|&p| {
                let v = match self.kind[p] {
                    -2 => x[p],
                    -1 => x[p] / self.diag[p],
                    _ => self.nucleus_dot_old(p, &x) / self.diag[p],
                };
                (p, v)
            })
            .collect()
    }


    /// The row-wise dot product of nucleus position `p` (`0x1241ed`): `x[p]` minus the products of the row's values with `x` at its columns,
    /// in the row's physical order.
    
    fn nucleus_dot_old(&self, p: usize, x: &[R]) -> R {
        let (cols, vals) = self.urows.row(p);
        let prod = |i: usize| vals[i] * x[cols[i] as usize];
        if R::QUAD {
            row_dot_q_by(x[p], cols.len(), prod)
        } else {
            row_dot_by(x[p], cols.len(), self.row_par[p] as usize, prod)
        }
    }
}
