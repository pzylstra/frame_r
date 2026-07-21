# Migrating frame's fire engine off the JVM: the `framecpp` bridge

_A proposal for discussion, not a merge request to act on unilaterally. This
branch is prepared so the work is ready to evaluate. — July 2026_

## The idea in one line

Replace frame's Scala-JVM + SQLite compute core with **`framecpp`**, a
self-contained C++/Rcpp reimplementation of the FRaME mechanics that runs
**in process** — same numbers, no `java`, no `.db` round-trip, ~500× faster per
run. frame keeps everything else.

## Why

frame's entire dependence on the Scala engine funnels through **one function**:
`ffm_run()`, which writes a parameter CSV, shells out to
`system("java -cp … ffm.runner.CSVRunner …")`, and reads a seven-table SQLite
database back. Everything else in frame (Building, Flora dynamics, Ensembles,
Risk, Fauna, Earth, Summaries) reads results from that database. So the seam to
replace is small and well defined — and replacing it removes the JVM/JDK/SQLite
dependency and makes the whole ecosystem far faster, which is what unlocks the
large-sample analyses (global sensitivity, ensembles, continuous-vs-discrete
studies) that are slow today.

`framecpp` reproduces the FRaME model **bit-exact** against frame v0.5.3 (worst
relative error ~3e-15 across the golden corpus — floating-point noise, not model
drift). It is a faithful port, not a reinterpretation: divergences from frame are
treated as bugs unless deliberately flagged.

## The package boundary we propose

**`framecpp` owns the physics and the data types; frame owns the domains.**

- **`framecpp`** = a slim, deterministic, single-site fire *engine* plus two typed
  objects: `fire_scenario` (validated input) and `fire_result` (the full
  seven-table run output). Nothing stochastic, nothing temporal, nothing that
  loops the engine, nothing domain-secondary.
- **frame** stays the home for the empirical ecosystem — flora dynamics,
  ensembles, risk, fauna, earth, FMC, weather, survey/trait sourcing — rebuilt on
  framecpp's fast engine instead of the JVM.

Litmus test, applied per function: *does it need to know fire physics (→ framecpp)
or a domain — weather, traits, fauna, soil, time, probability (→ frame)?* The
moment a function varies inputs, iterates over time, loops the engine, or maps
fire output onto a secondary effect, it is an empirical assessment and stays in
frame. Two consequences worth calling out:

- **Survey stratification (`frameSurvey`/`frameStratify`) stays in frame.** It is
  a domain model (k-means over transects, stochastic), so it belongs on the frame
  side and emits a scenario for the engine — not in the deterministic core.
- **The `plant` demographic model meets the engine directly.** framecpp provides a
  `scenario_from_plant()` front door, so `plant` needs no frame bridge and no
  k-means.

Dependency direction is one-way: both `plant` and frame depend on framecpp; framecpp
depends on neither.

## What this branch adds (the validation bridge)

A **flag-gated, additive** bridge — the JVM path is untouched and remains the
default:

```r
ffm_run(params, db.path = "out.db", db.recreate = TRUE, engine = "framecpp")
```

`engine = "framecpp"` routes to a new `ffm_run_framecpp()` (`R/framecpp_bridge.R`) that
builds a framecpp `fire_scenario` from the same `paramBuilder` table, runs the C++
engine in process, and writes the **identical** seven-table SQLite database at
`db.path`. Species in `IgnitionPaths` are resolved back to names from the
parameter table; `repId` append semantics match the JVM runner. So every
downstream analysis that calls `ffm_run()` then `ffm_db_load()` runs unchanged.

Requires the `framecpp` package (traitecoevo/framecpp).

## What is validated

Driving the bridge over frame's golden corpus (four communities: mild, extreme,
a hand-authored round-leaf community, and the survey-derived tuart community) and
diffing every one of the seven result tables, cell by cell, against the captured
**JVM** database:

- All seven tables bit-exact, **worst relative error ~2e-14**, including the
  `IgnitionPaths` species *names* and segment geometry.
- Downstream `frameSummary()` reads the bridge database back cleanly.

Harness: `framecpp/tools/parity/frame_bridge_parity.R` (`make bridge-parity` in the
framecpp repo).

## What is NOT yet tested — the honest gap

The corpus validates *single* runs. It does **not** yet validate frame's
higher-level workflows end to end:

- **Ensembles** (Monte Carlo looping the engine over many trait/weather draws),
- **Flora dynamics** (growth over time),
- **Risk / spotFire / probFire**, **Fauna**, **Earth** (secondary effects).

The next step is to run each of these whole workflows twice — `engine = "java"`
vs `engine = "framecpp"`, with a fixed seed so both see identical scenario draws —
and diff the *final* outputs. That is a far broader test than four communities:
it exercises the corners an ensemble reaches that the corpus never does. This is
scoped but not done; happy to do it.

## Suggested path (staged, your call)

1. **Broaden the diff** to live ensemble/risk/fauna/earth runs vs the JVM;
   classify any divergence (real vs sub-ULP boundary rounding). Parity thresholds
   are yours to set.
2. **Only once that is clean:** flip `engine` to default to `"framecpp"`, then
   retire the Java path (drop the bundled jars, `ffm_run_command`,
   `ffm_check_java`, and Java/SQLite from `SystemRequirements`). frame becomes
   pure R + framecpp.

Whether Scala is retired or kept as an independent reference is a project
decision — the bridge is designed so the JVM can stay as an oracle for as long as
useful.

## Companion housekeeping

An earlier version of this branch (targeting `master` / v0.5.3) carried two
behaviour-preserving source fixes so frame would run on a current
macOS/R stack: an `all()` wrap of the mandatory-columns check (R ≥ 4.3 errors on
`&&` with a length>1 operand) and a repair of the survey ingestion path broken by
commit `452bc91`. **`developframe` already resolves both** — the columns check is
rewritten as `prod(as.numeric(...)) == 1`, and `frameStratify`/`frameSurvey` use
the canonical column names throughout — so this PR drops them and targets
`developframe` directly. This PR is now purely additive.

Remaining portability items are environment/harness workarounds, not frame source
changes: the Windows `;` classpath separator, the Java-version gate, unqualified
dplyr verbs (need dplyr attached), and the bundled `sqlite-jdbc` lacking an
Apple-Silicon native lib.
