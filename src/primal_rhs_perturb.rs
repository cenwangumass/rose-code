//! x87 primal infeasible-basic repair (`0x27d7b0`).
use crate::f80::F80;

/// Applies the reference's correction in basis order. `random` supplies the
/// generator outputs, one for each A/B row; columns use native row indices.
/// `rhs` and `x` are extended precision, while the objective offset is stored
/// as double after each basic variable.
#[allow(dead_code)] // Integrated when the primal controller carries the x87 RHS state.
pub(crate) fn repair(
    basis: &[usize],
    types: &mut [u8],
    x: &mut [F80],
    dist: &[F80],
    rhs: &mut [F80],
    columns: &[Vec<(usize, f64)>],
    costs: &[f64],
    offset: &mut f64,
    size: f64,
    tol: f64,
    mut random: impl FnMut() -> u32,
) {
    let m = basis.len();
    for i in 0..m {
        let typ = types[i];
        if typ != b'A' && typ != b'B' {
            continue;
        }
        let half_range = dist[i] * F80::from_f64(0.5);
        let r = random() as usize % m;
        let k = if F80::from_f64(size) < half_range {
            F80::from_f64((r as f64 * size) / m as f64)
        } else {
            (F80::from_i64(r as i64) * half_range) / F80::from_i64(m as i64)
        };
        let old = x[i];
        let delta = if typ == b'B' {
            k - old
        } else {
            (dist[i] - old) - k
        };
        x[i] = old + delta;
        types[i] = if dist[i] == F80::ZERO {
            b'X'
        } else if typ == b'B' {
            if x[i] <= F80::from_f64(tol) {
                b'L'
            } else {
                b'N'
            }
        } else if dist[i] - x[i] <= F80::from_f64(tol) {
            b'U'
        } else {
            b'N'
        };
        if delta == F80::ZERO {
            continue;
        }
        let j = basis[i];
        *offset = (F80::from_f64(*offset) - F80::from_f64(costs[j]) * delta).to_f64();
        if j >= columns.len() {
            rhs[j - columns.len()] += delta;
        } else {
            for &(row, value) in &columns[j] {
                rhs[row] += F80::from_f64(value) * delta;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    fn run(data: &str) {
        let v: Value = serde_json::from_str(data).unwrap();
        let input = &v["PCIN"];
        let output = &v["PCOUT"];
        let ints = |v: &Value| {
            v.as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_u64().unwrap() as usize)
                .collect::<Vec<_>>()
        };
        let doubles = |v: &Value| {
            v.as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_f64().unwrap())
                .collect::<Vec<_>>()
        };
        let f80s = |v: &Value| {
            let s = v.as_str().unwrap();
            (0..s.len() / 32)
                .map(|i| {
                    let mut b = [0; 10];
                    for j in 0..10 {
                        b[j] =
                            u8::from_str_radix(&s[i * 32 + j * 2..i * 32 + j * 2 + 2], 16).unwrap();
                    }
                    F80::from_repr(
                        u64::from_le_bytes(b[..8].try_into().unwrap()),
                        u16::from_le_bytes(b[8..].try_into().unwrap()),
                    )
                })
                .collect::<Vec<_>>()
        };
        let basis = ints(&input["basis"]);
        let starts = ints(&input["starts"]);
        let lens = ints(&input["lens"]);
        let inds = ints(&input["inds"]);
        let vals = doubles(&input["vals"]);
        let costs = doubles(&input["cost"]);
        let cols = starts
            .iter()
            .zip(lens)
            .map(|(&s, n)| (s..s + n).map(|i| (inds[i], vals[i])).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let mut types = input["types"].as_str().unwrap().as_bytes().to_vec();
        let mut x = f80s(&input["c238"]);
        let dist = f80s(&input["dist"]);
        let mut rhs = f80s(&input["rhs"]);
        let mut offset = input["c130"].as_f64().unwrap();
        let draws = ints(&v["random"]);
        let mut draws = draws.into_iter();
        repair(
            &basis,
            &mut types,
            &mut x,
            &dist,
            &mut rhs,
            &cols,
            &costs,
            &mut offset,
            input["c340"].as_f64().unwrap(),
            input["c380"].as_f64().unwrap(),
            || draws.next().unwrap() as u32,
        );
        assert!(draws.next().is_none());
        assert_eq!(types, output["types"].as_str().unwrap().as_bytes());
        assert_eq!(
            x.iter().map(|x| x.repr()).collect::<Vec<_>>(),
            f80s(&output["c238"])
                .iter()
                .map(|x| x.repr())
                .collect::<Vec<_>>()
        );
        assert_eq!(
            rhs.iter().map(|x| x.repr()).collect::<Vec<_>>(),
            f80s(&output["rhs"])
                .iter()
                .map(|x| x.repr())
                .collect::<Vec<_>>()
        );
        assert_eq!(offset.to_bits(), output["c130"].as_f64().unwrap().to_bits());
    }
    #[test]
    fn restored_5959() {
        run(include_str!("fixtures/primal_quad_perturb_capture.json"));
    }
    #[test]
    fn injected_branches() {
        run(include_str!("fixtures/primal_quad_perturb_injected.json"));
    }
}
