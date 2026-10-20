//! Continuous-LP MPS reader (fixed and free format, tokenized). Integer markers, SOS and quadratic sections are rejected.

use crate::model::{LpModel, Sense, INF};
use std::collections::HashMap;

#[derive(Debug)]
pub struct MpsError(pub String);

impl std::fmt::Display for MpsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "MPS error: {}", self.0)
    }
}

impl std::error::Error for MpsError {}

fn err<T>(msg: impl Into<String>) -> Result<T, MpsError> {
    Err(MpsError(msg.into()))
}

fn num(tok: &str) -> Result<f64, MpsError> {
    tok.parse::<f64>().or_else(|_| tok.replace(['d', 'D'], "e").parse::<f64>()).map_err(|_| MpsError(format!("bad number {tok:?}")))
}

#[derive(Clone, Copy, PartialEq)]
enum RowKind {
    E,
    L,
    G,
}

pub fn parse(text: &str) -> Result<LpModel, MpsError> {
    let mut name = String::new();
    let mut section = String::new();
    let mut sense = Sense::Minimize;
    let mut obj_row: Option<String> = None;
    let mut rows: Vec<(String, RowKind)> = Vec::new();
    let mut row_ix: HashMap<String, usize> = HashMap::new();
    let mut ignored_rows: HashMap<String, ()> = HashMap::new();
    let mut col_names: Vec<String> = Vec::new();
    let mut col_ix: HashMap<String, usize> = HashMap::new();
    let mut cols: Vec<Vec<(usize, f64)>> = Vec::new();
    let mut obj: Vec<f64> = Vec::new();
    let mut rhs: Vec<f64> = Vec::new();
    let mut obj_rhs = 0.0;
    let mut range: Vec<Option<f64>> = Vec::new();
    let mut lo: Vec<f64> = Vec::new();
    let mut up: Vec<f64> = Vec::new();
    let mut objsense_next = false;

    for raw in text.lines() {
        if raw.trim().is_empty() || raw.starts_with('*') {
            continue;
        }
        let toks: Vec<&str> = raw.split_whitespace().collect();
        if !raw.starts_with(' ') && !raw.starts_with('\t') {
            section = toks[0].to_ascii_uppercase();
            match section.as_str() {
                "NAME" => name = toks.get(1).map(|s| s.to_string()).unwrap_or_default(),
                "OBJSENSE" => {
                    if let Some(t) = toks.get(1) {
                        sense = if t.to_ascii_uppercase().starts_with("MAX") { Sense::Maximize } else { Sense::Minimize };
                    } else {
                        objsense_next = true;
                    }
                }
                "ROWS" | "COLUMNS" | "RHS" | "RANGES" | "BOUNDS" | "ENDATA" | "OBJSENSE_" => {}
                other => return err(format!("unsupported section {other}")),
            }
            continue;
        }
        match section.as_str() {
            "OBJSENSE" => {
                if objsense_next || !toks.is_empty() {
                    sense = if toks[0].to_ascii_uppercase().starts_with("MAX") { Sense::Maximize } else { Sense::Minimize };
                    objsense_next = false;
                }
            }
            "ROWS" => {
                if toks.len() != 2 {
                    return err(format!("bad ROWS line {raw:?}"));
                }
                let (kind, rname) = (toks[0].to_ascii_uppercase(), toks[1].to_string());
                match kind.as_str() {
                    "N" => {
                        if obj_row.is_none() {
                            obj_row = Some(rname);
                        } else {
                            ignored_rows.insert(rname, ());
                        }
                    }
                    "E" | "L" | "G" => {
                        row_ix.insert(rname.clone(), rows.len());
                        rows.push((rname, match kind.as_str() { "E" => RowKind::E, "L" => RowKind::L, _ => RowKind::G }));
                        rhs.push(0.0);
                        range.push(None);
                    }
                    _ => return err(format!("bad row type {kind}")),
                }
            }
            "COLUMNS" => {
                if toks.iter().any(|t| t.eq_ignore_ascii_case("'MARKER'")) {
                    return err("integer markers are not supported (continuous LP only)");
                }
                if toks.len() != 3 && toks.len() != 5 {
                    return err(format!("bad COLUMNS line {raw:?}"));
                }
                let cname = toks[0];
                let j = *col_ix.entry(cname.to_string()).or_insert_with(|| {
                    col_names.push(cname.to_string());
                    cols.push(Vec::new());
                    obj.push(0.0);
                    lo.push(0.0);
                    up.push(INF);
                    cols.len() - 1
                });
                for pair in toks[1..].chunks(2) {
                    let v = num(pair[1])?;
                    if Some(pair[0]) == obj_row.as_deref() {
                        obj[j] += v;
                    } else if let Some(&i) = row_ix.get(pair[0]) {
                        cols[j].push((i, v));
                    } else if !ignored_rows.contains_key(pair[0]) {
                        return err(format!("unknown row {}", pair[0]));
                    }
                }
            }
            "RHS" | "RANGES" => {
                let body: &[&str] = if toks.len() % 2 == 1 { &toks[1..] } else { &toks[..] };
                for pair in body.chunks(2) {
                    if pair.len() != 2 {
                        return err(format!("bad {section} line {raw:?}"));
                    }
                    let v = num(pair[1])?;
                    if Some(pair[0]) == obj_row.as_deref() {
                        if section == "RHS" {
                            obj_rhs = v;
                        }
                    } else if let Some(&i) = row_ix.get(pair[0]) {
                        if section == "RHS" {
                            rhs[i] = v;
                        } else {
                            range[i] = Some(v);
                        }
                    } else if !ignored_rows.contains_key(pair[0]) {
                        return err(format!("unknown row {}", pair[0]));
                    }
                }
            }
            "BOUNDS" => {
                let kind = toks[0].to_ascii_uppercase();
                let valued = matches!(kind.as_str(), "UP" | "LO" | "FX" | "LI" | "UI" | "SC");
                let want = if valued { 3 } else { 2 };
                let (cname, val) = if toks.len() == want + 1 {
                    (toks[2], if valued { Some(num(toks[3])?) } else { None })
                } else if toks.len() == want {
                    (toks[1], if valued { Some(num(toks[2])?) } else { None })
                } else {
                    return err(format!("bad BOUNDS line {raw:?}"));
                };
                let j = *col_ix.get(cname).ok_or_else(|| MpsError(format!("unknown column {cname}")))?;
                match kind.as_str() {
                    "UP" => {
                        let v = val.unwrap();
                        up[j] = v;
                        if v < 0.0 && lo[j] == 0.0 {
                            lo[j] = -INF;
                        }
                    }
                    "LO" => lo[j] = val.unwrap(),
                    "FX" => {
                        lo[j] = val.unwrap();
                        up[j] = val.unwrap();
                    }
                    "FR" => {
                        lo[j] = -INF;
                        up[j] = INF;
                    }
                    "MI" => lo[j] = -INF,
                    "PL" => up[j] = INF,
                    "BV" => {
                        return err("binary variables are not supported (continuous LP only)");
                    }
                    "LI" | "UI" | "SC" => return err(format!("bound type {kind} is not supported (continuous LP only)")),
                    _ => return err(format!("bad bound type {kind}")),
                }
            }
            _ => return err(format!("data line outside a section: {raw:?}")),
        }
    }

    let m = rows.len();
    let mut row_lo = vec![-INF; m];
    let mut row_up = vec![INF; m];
    for (i, (_, kind)) in rows.iter().enumerate() {
        let r = rhs[i];
        let rg = range[i];
        match kind {
            RowKind::E => match rg {
                None => (row_lo[i], row_up[i]) = (r, r),
                Some(x) if x >= 0.0 => (row_lo[i], row_up[i]) = (r, r + x),
                Some(x) => (row_lo[i], row_up[i]) = (r + x, r),
            },
            RowKind::L => (row_lo[i], row_up[i]) = (rg.map_or(-INF, |x| r - x.abs()), r),
            RowKind::G => (row_lo[i], row_up[i]) = (r, rg.map_or(INF, |x| r + x.abs())),
        }
    }
    let mut col_start = vec![0usize];
    let mut row_idx = Vec::new();
    let mut val = Vec::new();
    for c in &cols {
        let mut c = c.clone();
        c.sort_by_key(|e| e.0);
        for (i, v) in c {
            row_idx.push(i);
            val.push(v);
        }
        col_start.push(row_idx.len());
    }
    Ok(LpModel {
        name,
        sense,
        offset: -obj_rhs,
        obj,
        col_lo: lo,
        col_up: up,
        row_lo,
        row_up,
        col_start,
        row_idx,
        val,
        col_names,
        row_names: rows.into_iter().map(|r| r.0).collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_ranges_bounds_and_offset() {
        let text = "NAME T\nROWS\n N obj\n L r1\n G r2\n E r3\nCOLUMNS\n x obj 1 r1 1\n x r2 1 r3 1\n y obj 2 r1 1\n y r3 1\nRHS\n rhs obj -5 r1 10\n rhs r2 1 r3 4\nRANGES\n rng r1 3 r2 2 r3 -1\nBOUNDS\n UP bnd x 4\n FR bnd y\nENDATA\n";
        let m = parse(text).unwrap();
        assert_eq!((m.num_rows(), m.num_cols()), (3, 2));
        assert_eq!(m.offset, 5.0);
        assert_eq!((m.row_lo[0], m.row_up[0]), (7.0, 10.0));
        assert_eq!((m.row_lo[1], m.row_up[1]), (1.0, 3.0));
        assert_eq!((m.row_lo[2], m.row_up[2]), (3.0, 4.0));
        assert_eq!((m.col_lo[1], m.col_up[1]), (-INF, INF));
        assert_eq!(m.col_up[0], 4.0);
    }

    #[test]
    fn rejects_integer_markers() {
        let text = "NAME T\nROWS\n N obj\nCOLUMNS\n MARKER 'MARKER' 'INTORG'\nENDATA\n";
        assert!(parse(text).is_err());
    }
}
