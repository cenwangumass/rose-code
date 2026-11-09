//! Self-parity gate: checks that a rose build makes bit for bit the same solves as a frozen baseline build.
//!
//! Rose agrees with the reference on the runs of `parity/corpus.txt`, so a build that reproduces the baseline's every dump record and result
//! on them keeps that agreement. This is the check to run after a change meant to leave the numbers alone (performance work, refactors);
//! it needs neither GDB nor the reference library.
//!
//! ```text
//! rose-parity check  [--tier fast|full] [--only SUBSTR,..] [--bin PATH] [--jobs N] [--ignore TAG,..] [--baseline-bin PATH] [--no-locate]
//! rose-parity record [--tier fast|full] [--only SUBSTR,..] [--bin PATH] [--jobs N] [--source TEXT]
//! rose-parity list   [--tier fast|full]
//! ```
//!
//! Each run of the corpus is `rose solve MODEL ARGS` with `ROSE_DUMP=all`. Every dump record (`TAG {json}` on stderr) is hashed into the
//! stream of its tag, in chunks of `CHUNK` records, and the tag sequence into an order stream; the result file (without its wall time) and the
//! exit code are hashed too. `record` stores these digests in `parity/baseline/ID.json` and copies the binary to `target/parity-baseline/rose`;
//! `check` recomputes them and compares. The fast tier is a subset of the full tier (`full` runs every corpus line).
//!
//! On a difference `check` names the differing streams and their first differing chunk, then (unless `--no-locate`) reruns the run with the
//! frozen baseline binary and the checked binary side by side and prints the first differing record with the path of the first differing
//! field. `--ignore TAG` leaves a tag's records out of both, for a change that deliberately alters a diagnostic-only record.
use serde_json::{json, Map, Value};
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::{mpsc, Arc, Mutex};
use std::time::Instant;

const CHUNK: u64 = 1024;

const ROOT: &str = env!("CARGO_MANIFEST_DIR");

struct Run {
    tier: String,
    id: String,
    model: String,
    args: Vec<String>,
}

fn usage() -> ! {
    eprintln!("usage: rose-parity check|record|list [--tier fast|full] [--only SUBSTR,..] [--bin PATH] [--jobs N] [--ignore TAG,..] [--baseline-bin PATH] [--no-locate] [--source TEXT]");
    std::process::exit(2)
}

fn root() -> PathBuf {
    PathBuf::from(ROOT)
}

fn corpus() -> Vec<Run> {
    let p = root().join("parity/corpus.txt");
    let text = std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
    let mut runs = Vec::new();
    for line in text.lines() {
        let line = line.split('#').next().unwrap().trim();
        if line.is_empty() {
            continue;
        }
        let t: Vec<&str> = line.split_whitespace().collect();
        assert!(t.len() >= 3, "corpus line `{line}`: tier id model [args]");
        assert!(
            t[0] == "fast" || t[0] == "full",
            "corpus line `{line}`: tier is fast or full"
        );
        runs.push(Run {
            tier: t[0].into(),
            id: t[1].into(),
            model: t[2].into(),
            args: t[3..].iter().map(|s| s.to_string()).collect(),
        });
    }
    let mut ids: Vec<&str> = runs.iter().map(|r| r.id.as_str()).collect();
    ids.sort();
    for w in ids.windows(2) {
        assert!(w[0] != w[1], "corpus: duplicate id {}", w[0]);
    }
    runs
}

/// A fast 64-bit hash of a byte string (change detection, not adversarial input).
fn hash(bytes: &[u8]) -> u64 {
    const K: u64 = 0x9e37_79b9_7f4a_7c15;
    let mut h = 0x243f_6a88_85a3_08d3 ^ (bytes.len() as u64).wrapping_mul(K);
    let mut chunks = bytes.chunks_exact(8);
    for c in &mut chunks {
        let w = u64::from_le_bytes(c.try_into().unwrap());
        h = (h ^ w).wrapping_mul(K).rotate_left(29);
    }
    let mut tail = [0u8; 8];
    tail[..chunks.remainder().len()].copy_from_slice(chunks.remainder());
    h = (h ^ u64::from_le_bytes(tail)).wrapping_mul(K);
    h ^ (h >> 32)
}

#[derive(Default)]
struct Stream {
    n: u64,
    acc: u64,
    chunks: Vec<String>,
}

impl Stream {
    fn push(&mut self, h: u64) {
        self.acc = (self.acc ^ h).wrapping_mul(0x100_0000_01b3).rotate_left(17) ^ self.n;
        self.n += 1;
        if self.n % CHUNK == 0 {
            self.flush();
        }
    }

    fn flush(&mut self) {
        self.chunks.push(format!("{:016x}", self.acc));
        self.acc = 0;
    }

    fn finish(mut self) -> Value {
        if self.n % CHUNK != 0 {
            self.flush();
        }
        json!({ "n": self.n, "chunks": self.chunks })
    }
}

/// The dump tag of a stderr line, when it is a dump record (`TAG {json}`).
fn tag(line: &str) -> Option<&str> {
    let (t, rest) = line.split_once(' ')?;
    (!t.is_empty()
        && t.bytes()
            .all(|b| b.is_ascii_uppercase() || b.is_ascii_digit() || b == b'_')
        && rest.starts_with('{'))
    .then_some(t)
}

fn spawn(bin: &Path, run: &Run, out: &Path) -> Child {
    Command::new(bin)
        .current_dir(root())
        .arg("solve")
        .arg(&run.model)
        .args(&run.args)
        .arg("--output")
        .arg(out)
        .env("ROSE_DUMP", "all")
        .env_remove("ROSE_FROM")
        .env_remove("ROSE_WORK")
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap_or_else(|e| panic!("{}: {e}", bin.display()))
}

/// The result file without its wall-clock fields, canonically serialized.
fn result_digest(out: &Path) -> String {
    let Ok(text) = std::fs::read_to_string(out) else {
        return "none".into();
    };
    let mut v: Value = match serde_json::from_str(&text) {
        Ok(v) => v,
        Err(_) => return format!("{:016x}", hash(text.as_bytes())),
    };
    strip_time(&mut v);
    format!(
        "{:016x}",
        hash(serde_json::to_string(&v).unwrap().as_bytes())
    )
}

fn strip_time(v: &mut Value) {
    match v {
        Value::Object(m) => {
            m.remove("seconds");
            m.values_mut().for_each(strip_time);
        }
        Value::Array(a) => a.iter_mut().for_each(strip_time),
        _ => {}
    }
}

fn scratch(tag: &str) -> PathBuf {
    let d = root().join("target/parity-tmp");
    std::fs::create_dir_all(&d).unwrap();
    d.join(format!(
        "{tag}-{}-{:?}.json",
        std::process::id(),
        std::thread::current().id()
    ))
}

/// Runs `run` with `bin` and returns its digest.
fn digest(bin: &Path, run: &Run) -> Value {
    let out = scratch(&run.id.replace('/', "_"));
    let _ = std::fs::remove_file(&out);
    let t0 = Instant::now();
    let mut child = spawn(bin, run, &out);
    let mut tags: BTreeMap<String, Stream> = BTreeMap::new();
    let mut order = Stream::default();
    let mut rd = BufReader::with_capacity(1 << 20, child.stderr.take().unwrap());
    let mut line = String::new();
    let mut records = 0u64;
    loop {
        line.clear();
        if rd.read_line(&mut line).unwrap() == 0 {
            break;
        }
        let l = line.trim_end_matches('\n');
        let Some(t) = tag(l) else { continue };
        let h = hash(l.as_bytes());
        tags.entry(t.to_string()).or_default().push(h);
        order.push(hash(t.as_bytes()));
        records += 1;
    }
    let status = child.wait().unwrap();
    let ms = t0.elapsed().as_millis() as u64;
    let result = result_digest(&out);
    let _ = std::fs::remove_file(&out);
    json!({
        "id": run.id, "model": run.model, "args": run.args,
        "exit": status.code(), "records": records, "ms": ms, "result": result,
        "order": order.finish(),
        "tags": tags.into_iter().map(|(k, s)| (k, s.finish())).collect::<Map<_, _>>(),
    })
}

/// The differences between a baseline digest and a new one: (stream, first differing record index or None for whole-value fields).
fn differences(base: &Value, new: &Value, ignore: &[String]) -> Vec<String> {
    let mut d = Vec::new();
    for k in ["exit", "result"] {
        if base[k] != new[k] {
            d.push(format!("{k}: {} -> {}", base[k], new[k]));
        }
    }
    let stream = |name: &str, a: &Value, b: &Value, d: &mut Vec<String>| {
        let (na, nb) = (a["n"].as_u64().unwrap_or(0), b["n"].as_u64().unwrap_or(0));
        let empty = vec![];
        let ca = a["chunks"].as_array().unwrap_or(&empty);
        let cb = b["chunks"].as_array().unwrap_or(&empty);
        let first = (0..ca.len().max(cb.len())).find(|&i| ca.get(i) != cb.get(i));
        if let Some(i) = first {
            let lo = i as u64 * CHUNK;
            d.push(format!(
                "{name}: {na} -> {nb} records, first difference in records {lo}..{}",
                lo + CHUNK
            ));
        }
    };
    if ignore.is_empty() {
        stream("record order", &base["order"], &new["order"], &mut d);
    }
    let null = Value::Null;
    let (bt, nt) = (base["tags"].as_object(), new["tags"].as_object());
    let mut names: Vec<&String> = bt.iter().chain(nt.iter()).flat_map(|m| m.keys()).collect();
    names.sort();
    names.dedup();
    for t in names {
        if ignore.contains(t) {
            continue;
        }
        let a = bt.and_then(|m| m.get(t)).unwrap_or(&null);
        let b = nt.and_then(|m| m.get(t)).unwrap_or(&null);
        stream(t, a, b, &mut d);
    }
    d
}

/// The path and the two values of the first difference between two JSON values.
fn first_diff(a: &Value, b: &Value, path: String) -> Option<(String, String, String)> {
    match (a, b) {
        (Value::Object(x), Value::Object(y)) => {
            let mut keys: Vec<&String> = x.keys().chain(y.keys()).collect();
            keys.sort();
            keys.dedup();
            keys.into_iter().find_map(|k| match (x.get(k), y.get(k)) {
                (Some(p), Some(q)) => first_diff(p, q, format!("{path}.{k}")),
                (p, q) => Some((format!("{path}.{k}"), short(p), short(q))),
            })
        }
        (Value::Array(x), Value::Array(y)) => {
            let found = x
                .iter()
                .zip(y)
                .enumerate()
                .find_map(|(i, (p, q))| first_diff(p, q, format!("{path}[{i}]")));
            found.or_else(|| {
                (x.len() != y.len()).then(|| {
                    (
                        format!("{path}.len"),
                        x.len().to_string(),
                        y.len().to_string(),
                    )
                })
            })
        }
        (Value::Number(x), Value::Number(y)) if x.is_f64() || y.is_f64() => {
            (x.as_f64().map(f64::to_bits) != y.as_f64().map(f64::to_bits))
                .then(|| (path, short(Some(a)), short(Some(b))))
        }
        _ => (a != b).then(|| (path, short(Some(a)), short(Some(b)))),
    }
}

fn short(v: Option<&Value>) -> String {
    let s = v.map_or("(absent)".to_string(), |v| v.to_string());
    if s.len() > 200 {
        format!("{}...", &s[..s.floor_char_boundary(200)])
    } else {
        s
    }
}

/// Runs both binaries side by side and reports the first differing dump record (or the result when all records agree).
fn locate(base_bin: &Path, bin: &Path, run: &Run, ignore: &[String]) -> String {
    let (oa, ob) = (scratch("loc-a"), scratch("loc-b"));
    let (mut ca, mut cb) = (spawn(base_bin, run, &oa), spawn(bin, run, &ob));
    let records = |c: &mut Child| {
        let ignore = ignore.to_vec();
        BufReader::with_capacity(1 << 20, c.stderr.take().unwrap())
            .lines()
            .map_while(Result::ok)
            .filter(move |l| tag(l).is_some_and(|t| !ignore.iter().any(|i| i == t)))
    };
    let (mut ra, mut rb) = (records(&mut ca), records(&mut cb));
    let mut k = 0u64;
    let report = loop {
        match (ra.next(), rb.next()) {
            (None, None) => break None,
            (Some(a), None) => {
                break Some(format!(
                    "record {k}: baseline has {}, checked build has no more records",
                    head(&a)
                ))
            }
            (None, Some(b)) => {
                break Some(format!(
                    "record {k}: checked build has {}, baseline has no more records",
                    head(&b)
                ))
            }
            (Some(a), Some(b)) if a != b => {
                let (ta, tb) = (tag(&a).unwrap().to_string(), tag(&b).unwrap().to_string());
                if ta != tb {
                    break Some(format!(
                        "record {k}: baseline {}, checked build {}",
                        head(&a),
                        head(&b)
                    ));
                }
                let pa: Value = serde_json::from_str(&a[ta.len() + 1..]).unwrap_or(Value::Null);
                let pb: Value = serde_json::from_str(&b[tb.len() + 1..]).unwrap_or(Value::Null);
                let (p, x, y) = first_diff(&pa, &pb, String::new()).unwrap_or_default();
                break Some(format!(
                    "record {k} ({ta}): first difference at `{p}`: baseline {x}, checked build {y}\n      baseline record: {}",
                    head(&a)
                ));
            }
            _ => k += 1,
        }
    };
    drop((ra, rb));
    let _ = (ca.kill(), cb.kill(), ca.wait(), cb.wait());
    let r = match report {
        Some(r) => r,
        None => {
            let read = |p: &Path| {
                let mut v: Value = std::fs::read_to_string(p)
                    .ok()
                    .and_then(|t| serde_json::from_str(&t).ok())
                    .unwrap_or(Value::Null);
                strip_time(&mut v);
                v
            };
            match first_diff(&read(&oa), &read(&ob), String::new()) {
                Some((p, x, y)) => format!("all {k} records agree; result differs at `{p}`: baseline {x}, checked build {y}"),
                None => format!("all {k} records and the result agree on rerun (exit code differs?)"),
            }
        }
    };
    let _ = (std::fs::remove_file(&oa), std::fs::remove_file(&ob));
    r
}

fn head(line: &str) -> String {
    let t = tag(line).unwrap_or("?");
    format!(
        "{t} {}",
        short(Some(&Value::String(line[t.len()..].trim().to_string())))
    )
}

/// Runs `f` over `runs` on `jobs` threads, longest expected first, calling `done` on the main thread as results arrive.
fn pool<T: Send + 'static>(
    runs: Vec<Arc<Run>>,
    jobs: usize,
    f: impl Fn(&Run) -> T + Send + Sync + 'static,
    mut done: impl FnMut(&Run, T),
) {
    let queue = Arc::new(Mutex::new(runs));
    let f = Arc::new(f);
    let (tx, rx) = mpsc::channel();
    let threads: Vec<_> = (0..jobs)
        .map(|_| {
            let (queue, f, tx) = (queue.clone(), f.clone(), tx.clone());
            std::thread::spawn(move || loop {
                let Some(run) = queue.lock().unwrap().pop() else {
                    break;
                };
                let r = f(&run);
                tx.send((run, r)).unwrap();
            })
        })
        .collect();
    drop(tx);
    for (run, r) in rx {
        done(&run, r);
    }
    threads.into_iter().for_each(|t| t.join().unwrap());
}

fn baseline_path(id: &str) -> PathBuf {
    root().join("parity/baseline").join(format!("{id}.json"))
}

fn load(p: &Path) -> Option<Value> {
    serde_json::from_str(&std::fs::read_to_string(p).ok()?).ok()
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some(cmd) = args.first().cloned() else {
        usage()
    };
    let mut tier = "fast".to_string();
    let mut only: Vec<String> = Vec::new();
    let mut bin = root().join("target/release/rose");
    let mut base_bin = root().join("target/parity-baseline/rose");
    let mut jobs = std::thread::available_parallelism().map_or(4, |n| n.get());
    let mut ignore: Vec<String> = Vec::new();
    let mut locate_diffs = true;
    let mut source = String::new();
    let mut it = args[1..].iter();
    while let Some(a) = it.next() {
        let mut val = || it.next().cloned().unwrap_or_else(|| usage());
        match a.as_str() {
            "--tier" => tier = val(),
            "--only" => only = val().split(',').map(String::from).collect(),
            "--bin" => bin = PathBuf::from(val()),
            "--baseline-bin" => base_bin = PathBuf::from(val()),
            "--jobs" => jobs = val().parse().unwrap_or_else(|_| usage()),
            "--ignore" => ignore = val().split(',').map(String::from).collect(),
            "--no-locate" => locate_diffs = false,
            "--source" => source = val(),
            _ => usage(),
        }
    }
    if tier != "fast" && tier != "full" {
        usage();
    }
    // runs start in the crate root, so a relative path is made absolute against the current directory first
    let absolute = |p: PathBuf| std::path::absolute(&p).unwrap_or(p);
    let (bin, base_bin) = (absolute(bin), absolute(base_bin));
    let selected: Vec<Run> = corpus()
        .into_iter()
        .filter(|r| tier == "full" || r.tier == "fast")
        .filter(|r| only.is_empty() || only.iter().any(|o| r.id.contains(o.as_str())))
        .collect();
    // longest first: the pool pops from the end
    let expected = |r: &Run| {
        load(&baseline_path(&r.id))
            .and_then(|v| v["ms"].as_u64())
            .unwrap_or(u64::MAX)
    };
    let mut runs: Vec<Arc<Run>> = selected.into_iter().map(Arc::new).collect();
    runs.sort_by_key(|r| expected(r));
    let total = runs.len();
    match cmd.as_str() {
        "list" => {
            for r in runs.iter().rev() {
                println!(
                    "{} {} {} {} (baseline {} ms)",
                    r.tier,
                    r.id,
                    r.model,
                    r.args.join(" "),
                    expected(r)
                );
            }
        }
        "record" => {
            assert!(
                bin.exists(),
                "{}: build it first (cargo build --release)",
                bin.display()
            );
            std::fs::create_dir_all(root().join("parity/baseline")).unwrap();
            let t0 = Instant::now();
            let b = bin.clone();
            let mut n = 0;
            pool(
                runs,
                jobs,
                move |r| digest(&b, r),
                |run, d| {
                    n += 1;
                    std::fs::write(
                        baseline_path(&run.id),
                        serde_json::to_string_pretty(&d).unwrap() + "\n",
                    )
                    .unwrap();
                    println!(
                        "[{n}/{total}] {} {} records, {} ms",
                        run.id, d["records"], d["ms"]
                    );
                },
            );
            let same = |a: &Path, b: &Path| {
                a.canonicalize()
                    .ok()
                    .is_some_and(|a| b.canonicalize().ok() == Some(a))
            };
            if !same(&bin, &base_bin) {
                std::fs::create_dir_all(base_bin.parent().unwrap()).unwrap();
                std::fs::copy(&bin, &base_bin).unwrap();
            }
            let meta = json!({ "source": source, "binary_hash": format!("{:016x}", hash(&std::fs::read(&base_bin).unwrap())) });
            std::fs::write(
                root().join("parity/baseline/META.json"),
                serde_json::to_string_pretty(&meta).unwrap() + "\n",
            )
            .unwrap();
            println!(
                "recorded {total} runs in {:.1} s; baseline binary {}",
                t0.elapsed().as_secs_f64(),
                base_bin.display()
            );
        }
        "check" => {
            assert!(
                bin.exists(),
                "{}: build it first (cargo build --release)",
                bin.display()
            );
            let t0 = Instant::now();
            let (b, ig) = (bin.clone(), ignore.clone());
            let mut failed: Vec<(Arc<Run>, Vec<String>)> = Vec::new();
            let mut n = 0;
            pool(
                runs.clone(),
                jobs,
                move |r| {
                    let base = load(&baseline_path(&r.id));
                    let d = digest(&b, r);
                    let diffs = match &base {
                        Some(base) => differences(base, &d, &ig),
                        None => vec!["no baseline (run `rose-parity record --only ID`)".into()],
                    };
                    (diffs, d["ms"].as_u64().unwrap_or(0))
                },
                |run, (diffs, ms)| {
                    n += 1;
                    if diffs.is_empty() {
                        println!("[{n}/{total}] ok   {} ({ms} ms)", run.id);
                    } else {
                        println!("[{n}/{total}] DIFF {}", run.id);
                        for d in &diffs {
                            println!("       {d}");
                        }
                        let r = runs.iter().find(|r| r.id == run.id).unwrap().clone();
                        failed.push((r, diffs));
                    }
                },
            );
            let secs = t0.elapsed().as_secs_f64();
            if failed.is_empty() {
                println!("parity: all {total} runs agree with the baseline ({secs:.1} s)");
                return;
            }
            failed.sort_by(|a, b| a.0.id.cmp(&b.0.id));
            println!(
                "parity: {} of {total} runs differ ({secs:.1} s)",
                failed.len()
            );
            if locate_diffs {
                if base_bin.exists() {
                    for (r, _) in &failed {
                        println!("  {}: {}", r.id, locate(&base_bin, &bin, r, &ignore));
                    }
                } else {
                    println!(
                        "  (no {}: cannot locate the first differing record)",
                        base_bin.display()
                    );
                }
            }
            std::process::exit(1);
        }
        _ => usage(),
    }
}

#[cfg(test)]
mod strict_diff_tests {
    use super::*;

    #[test]
    fn reports_signed_zero_field_and_preserves_integer_identity() {
        let a: Value = serde_json::from_str(r#"{"tau":[[32,0.0]]}"#).unwrap();
        let b: Value = serde_json::from_str(r#"{"tau":[[32,-0.0]]}"#).unwrap();
        let diff = first_diff(&a, &b, String::new()).unwrap();
        assert_eq!(diff.0, ".tau[0][1]");
        assert!(first_diff(&a, &a, String::new()).is_none());
        let a = json!(9007199254740992u64);
        let b = json!(9007199254740993u64);
        assert!(first_diff(&a, &b, String::new()).is_some());
    }
}
