//! Repeated reduced costs at a retained numerical factor state.
use super::{Internal, NONE};
use crate::f80::F80;
use std::cell::OnceCell;

#[derive(Default)]
pub(super) struct CompactCosts(OnceCell<Costs>);

struct Costs {
    y: Vec<F80>,
    packed: Vec<[u8; 16]>,
    active: u64,
    costfree: bool,
    d: Vec<f64>,
    df: Vec<F80>,
    required_flags: u32,
    status: u16,
}

struct Environment {
    flags: u32,
    status: u16,
}

fn environment() -> Option<Environment> {
    // FXSAVE records both control words, the x87 stack occupancy and both
    // sticky exception masks without changing the caller's environment.
    #[repr(align(16))]
    struct Image([u8; 512]);
    let mut state = std::mem::MaybeUninit::<Image>::uninit();
    unsafe {
        std::arch::asm!("fxsave64 [{p}]",p=in(reg)state.as_mut_ptr(),options(nostack,preserves_flags));
    }
    let p = unsafe { std::ptr::addr_of!((*state.as_ptr()).0).cast::<u8>() };
    // These four fields are always written by FXSAVE. The reserved bytes of
    // the image may be uninitialized and must not be materialized as values.
    let (cw, sw, tag, csr) = unsafe {
        (
            p.cast::<u16>().read(),
            p.add(2).cast::<u16>().read(),
            p.add(4).read(),
            p.add(24).cast::<u32>().read(),
        )
    };
    if cw != 0x037f || csr & !63 != 0x1f80 || tag != 0 {
        return None;
    }
    Some(Environment {
        flags: (u32::from(sw & 63) << 6) | (csr & 63),
        status: sw & !63,
    })
}

impl CompactCosts {
    pub(super) fn apply(
        &self,
        p: &Internal,
        y: &[F80],
        pos: &[usize],
        excluded: Option<usize>,
        costfree: bool,
        d: &mut [f64],
        df: &mut [F80],
    ) -> bool {
        let n = p.n + p.m;
        if n > 64
            || y.len() != p.m
            || pos.len() != n
            || d.len() != n
            || df.len() != n
            || crate::dump::any()
        {
            return false;
        }
        let Some(before) = environment() else {
            return false;
        };
        let mut active = 0u64;
        for (j, &slot) in pos.iter().enumerate() {
            if slot == NONE && Some(j) != excluded {
                active |= 1 << j;
            }
        }
        if active.count_ones() <= 1 {
            return false;
        }
        if let Some(entry) = self.0.get() {
            if entry.active != active
                || entry.costfree != costfree
                || entry.required_flags & !before.flags != 0
                || entry.status != before.status
                || !entry.y.iter().zip(y).all(|(a, b)| a.same_bits(*b))
            {
                return false;
            }
            d.copy_from_slice(&entry.d);
            let last = (63 - active.leading_zeros()) as usize;
            let mut saved = active & !(1u64 << last);
            while saved != 0 {
                let j = saved.trailing_zeros() as usize;
                df[j] = entry.df[j];
                saved &= saved - 1;
            }
            p.write_reduced_costs(&entry.packed, pos, excluded, costfree, d, df, last);
            return true;
        }
        let packed = Internal::pack_dual(y);
        p.write_reduced_costs(&packed, pos, excluded, costfree, d, df, 0);
        if let Some(after) = environment() {
            if after.status == before.status {
                self.0
                    .set(Costs {
                        y: y.to_vec(),
                        packed,
                        active,
                        costfree,
                        d: d.to_vec(),
                        df: df.to_vec(),
                        required_flags: after.flags,
                        status: after.status,
                    })
                    .ok()
                    .unwrap();
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::arch::asm;

    fn set_cw(cw: u16) {
        unsafe {
            asm!("fldcw [{cw}]",cw=in(reg)&cw,options(nostack,preserves_flags));
        }
    }
    fn set_csr(csr: u32) {
        unsafe {
            asm!("ldmxcsr [{csr}]",csr=in(reg)&csr,options(nostack,preserves_flags));
        }
    }

    #[test]
    fn repeated_native_costs_preserve_partial_writes_and_environment() {
        let mut original_cw = 0u16;
        let mut original_csr = 0u32;
        unsafe {
            asm!("fnstcw [{cw}]", "stmxcsr [{csr}]",cw=in(reg)&mut original_cw,csr=in(reg)&mut original_csr,options(nostack,preserves_flags));
        }
        struct Restore(u16, u32);
        impl Drop for Restore {
            fn drop(&mut self) {
                set_cw(self.0);
                set_csr(self.1);
            }
        }
        let _restore = Restore(original_cw, original_csr);
        set_cw(0x037f);
        set_csr(0x1fbf);
        let model = crate::mps::parse(include_str!("../tests/fixtures/quad_2717.mps")).unwrap();
        let matrix = super::super::build(&model);
        let n = matrix.n + matrix.m;
        assert!(n <= 64);
        let mut seed = 0xa57ed47219ce3610u64;
        for round in 0..120 {
            let mut y = vec![F80::ZERO; matrix.m];
            for (i, value) in y.iter_mut().enumerate() {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                *value = match (i + round) % 12 {
                    0 => F80::ZERO,
                    1 => -F80::ZERO,
                    2 => F80::from_repr(1, 0),
                    3 => F80::from_repr(u64::MAX, 0x7ffe),
                    4 => F80::INFINITY,
                    5 => -F80::NAN,
                    _ => F80::from_repr(seed | (1 << 63), (0x3f80 + seed % 256) as u16),
                };
            }
            let pos: Vec<_> = (0..n)
                .map(|j| {
                    if (round + j) % 3 == 0 {
                        j % matrix.m
                    } else {
                        NONE
                    }
                })
                .collect();
            let ex = Some(round % n);
            let cf = round & 1 == 1;
            let owner = CompactCosts::default();
            let mut baseline_d = vec![-0.; n];
            let mut baseline_df = vec![F80::ONE; n];
            // Prime native flags and the condition codes before recording.
            assert!(matrix.reduced_costs_q_native(
                &y,
                &pos,
                ex,
                cf,
                &mut baseline_d,
                &mut baseline_df
            ));
            assert!(owner.apply(&matrix, &y, &pos, ex, cf, &mut baseline_d, &mut baseline_df));
            assert!(owner.0.get().is_some());
            for k in 0..3 {
                let mut new_d = vec![f64::NAN; n];
                let mut old_d = new_d.clone();
                let mut new_df: Vec<_> = (0..n)
                    .map(|j| F80::from_repr((1 << 63) + (j as u64) + k, 0x3fff))
                    .collect();
                let mut old_df = new_df.clone();
                assert!(matrix.reduced_costs_q_native(&y, &pos, ex, cf, &mut old_d, &mut old_df));
                let expected = environment().unwrap();
                assert!(owner.apply(&matrix, &y, &pos, ex, cf, &mut new_d, &mut new_df));
                let actual = environment().unwrap();
                assert_eq!(
                    (actual.flags, actual.status),
                    (expected.flags, expected.status)
                );
                assert!(old_d
                    .iter()
                    .zip(&new_d)
                    .all(|(a, b)| a.to_bits() == b.to_bits()));
                assert!(old_df.iter().zip(&new_df).all(|(a, b)| a.same_bits(*b)));
            }
            let mut d = vec![42.; n];
            let mut df = vec![F80::ONE; n];
            assert!(!owner.apply(&matrix, &y, &pos, ex, !cf, &mut d, &mut df));
            let mut other = y.clone();
            other[0] = -other[0];
            assert!(!owner.apply(&matrix, &other, &pos, ex, cf, &mut d, &mut df));
            let mut mask = pos.clone();
            let j = (0..n).find(|&j| Some(j) != ex).unwrap();
            mask[j] = if mask[j] == NONE { 0 } else { NONE };
            assert!(!owner.apply(&matrix, &y, &mask, ex, cf, &mut d, &mut df));
            unsafe {
                asm!("fld1", options(nostack, preserves_flags));
            }
            let occupied = owner.apply(&matrix, &y, &pos, ex, cf, &mut d, &mut df);
            unsafe {
                asm!("fstp st(0)",out("st(0)") _,options(nostack,preserves_flags));
            }
            assert!(!occupied);
            for csr in [0x3fbf, 0x9fbf, 0x1f7f] {
                set_csr(csr);
                assert!(!owner.apply(&matrix, &y, &pos, ex, cf, &mut d, &mut df));
            }
            set_csr(0x1fbf);
            set_cw(0x077f);
            assert!(!owner.apply(&matrix, &y, &pos, ex, cf, &mut d, &mut df));
            set_cw(0x037f);
            set_csr(0x1f80);
            unsafe {
                asm!("fnclex", options(nostack, preserves_flags));
            }
            assert!(!owner.apply(&matrix, &y, &pos, ex, cf, &mut d, &mut df));
            assert!(d.iter().all(|&v| v == 42.) && df.iter().all(|v| v.same_bits(F80::ONE)));
            set_csr(0x1fbf);
        }
    }
}
