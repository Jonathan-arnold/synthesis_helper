# synthesis_helper

A bottom-up retrobiosynthesis tool that builds a synthesis hypergraph from a
reaction corpus (MetaCyc) and enumerates **cascades, pathways, and enzyme
compositions** for any reachable target chemical. The whole thing is wrapped
as a [Model Context Protocol](https://modelcontextprotocol.io/) server so
Claude Code can drive it as an interactive lab notebook.

This document is organized around the two pieces of work that took the most
effort:

1. **What the MCP layer does and how it is implemented.** — *Benson Li*
2. **How we validated that `synthesize` and the cascade traceback are correct.** — *Benson Li*

The underlying algorithm (BFS waveform expansion → cascade traceback →
pathway enumeration → enzyme composition) is documented in `CLAUDE.md`.

> **Status:** sections 1 and 2 below were written and implemented by
> Benson Li. This file is being handed to the rest of the team — please
> append your sections (algorithm internals, data preparation, downstream
> integration, etc.) below section 2 and add yourselves to the relevant
> headers.

---

## 1. The MCP layer
*Author: Benson Li*

### 1.1 What the server exposes

The `synthesis_helper.mcp` package wraps the core API as a FastMCP server
over stdio. Once registered with Claude Code, the model gets the following
**tools** and **resources**.

#### Tools

| Tool | What it does |
|---|---|
| `find_chemical(query, limit)` | Look up a chemical by id, name substring, or InChI. |
| `find_by_structure(smiles, limit, similarity_threshold)` | Tanimoto search over Morgan fingerprints (radius 2, 2048 bits) for "find me chemicals shaped like this SMILES". |
| `get_shell(chemical_ref)` | Minimum number of reaction steps from native metabolites. `0` = native/universal, `null` = unreachable. |
| `list_reachables(max_shell, name_contains, limit)` | Filter the reachables table. |
| `get_cascade(chemical_ref, max_producers_per_chemical)` | The full producer tree (all reactions transitively reaching the target). |
| `enumerate_pathways_for(chemical_ref, max_pathways)` | Flatten the cascade into individual linear pathways. |
| `describe_pathway(chemical_ref, pathway_index)` | Markdown narrative of a single pathway. |
| `pathway_to_composition(chemical_ref, pathway_index)` | Per-step enzyme list with engineering flags (`is_orphan`, `is_p450`, `is_heme`, `ec_class`). |
| `compare_pathways(chemical_ref, n)` | Multi-dimensional scorecard ranking up to N pathways by step count, heme dependence, toxic intermediates, orphan steps, and net cofactor balance (NADH / NADPH / ATP). |
| `open_pathway_interactive(chemical_ref, pathway_index)` | Render one pathway as a self-contained interactive HTML page (Cytoscape.js, offline) and open it in the browser. |
| `open_cascade_interactive(chemical_ref, max_reactions)` | Same, for the full cascade. |
| `resynthesize_with_fed(fed_chemical_refs)` | Re-run BFS with extra fed chemicals added to shell 0; returns the delta vs. baseline. |

#### Resources

| URI | Contents |
|---|---|
| `synthesis://stats` | Counts, load time, data directory. |
| `synthesis://reachables` | TSV dump of every reachable chemical. |
| `synthesis://chemical/{id}` | Chemical detail + `produced_by` / `consumed_by` reaction ids. |
| `synthesis://reaction/{id}` | Reaction with expanded substrate / product names. |

### 1.2 How it is implemented

The MCP layer is split into a thin server and a few helpers so the heavy
work happens once per process and the tool functions stay shallow.

```
synthesis_helper/mcp/
├── server.py        # FastMCP @tool / @resource definitions
├── state.py         # Lazy singletons: chemicals, reactions, hypergraph,
│                    #   fingerprint index, EC-name table, cofactor sets,
│                    #   toxic-intermediate map, fed-hypergraph cache
├── lookup.py        # resolve_one() — id / exact-name / InChI / fuzzy match
├── similarity.py    # RDKit Morgan fingerprints + Tanimoto top-N search
├── serializers.py   # Domain object → JSON-friendly DTOs
├── html_render.py   # Cytoscape.js graph generation for the viewer tools
└── assets/
    └── cytoscape.min.js   # Inlined offline copy (fetched on first install)
```

Key design choices:

- **Lazy singleton state.** The hypergraph takes a few seconds to build from
  ~9k chemicals and ~9k reactions. `state.py` builds it on the first tool
  call and caches it in process memory; every subsequent tool call is
  effectively a dictionary lookup. The fingerprint index for
  `find_by_structure` follows the same pattern (built only when the first
  structure query arrives).
- **Fed-hypergraph cache.** `resynthesize_with_fed` keys its cache on the
  sorted tuple of fed chemical ids, so re-asking the same "what if I feed
  PABA?" question costs nothing after the first call.
- **Reference resolution at the boundary.** Every tool that takes a
  `chemical_ref` runs it through `lookup.resolve_one`, which accepts
  integer ids, exact names, InChI strings, or substring matches with a
  single confident hit. This means callers (and Claude) never have to
  pre-canonicalize identifiers.
- **DTO serializers, not raw dataclasses.** `serializers.py` converts
  `Chemical` / `Reaction` / `Cascade` / `Pathway` into plain dicts with
  shell numbers attached. The MCP wire format is JSON; keeping that
  conversion in one place means the domain objects can evolve without
  breaking the MCP contract.
- **Offline interactive viewer.** `open_pathway_interactive` and
  `open_cascade_interactive` write a self-contained HTML file (Cytoscape.js
  is inlined) to `/tmp` and shell out to `open` / `xdg-open` / `start` to
  launch the user's default browser. The page supports five layout
  algorithms (Layered, Breadth-first, Concentric, Force, Grid), three
  group-by themes (shell, EC class, role), a cofactor mode toggle (Inline /
  Shared / Hide), click-to-focus on a 1-hop neighborhood, and double-click
  to copy a node's structured details to the clipboard. All of this runs
  client-side — no server round-trip on layout / theme switch.
- **Engineering-aware pathway scoring.** `compare_pathways` is more than a
  step counter. It flags P450s (`EC 1.14.*`), broader heme-dependent
  enzymes (`1.11.1.*`, `1.13.11.*`), orphan steps with no EC mapping, and
  reactions that touch a curated blacklist of toxic / reactive
  intermediates (catechol, hydroquinone, reactive aldehydes, …). It also
  computes a crude net balance for ATP / NADH / NADPH so the user can spot
  pathways that drain biosynthetic reducing power.
- **Side-effect honesty in tool docstrings.** The `open_*_interactive`
  tools spell out, in the docstring the LLM sees, that they write a file
  and open a browser tab on the user's machine — and that retrying to
  "verify" the result will just open another tab. This is a small fix but
  a critical one for usability when the model can't see the user's
  desktop.

### 1.3 Install and register

```bash
python3 -m venv ~/.venvs/synthesis_helper
~/.venvs/synthesis_helper/bin/pip install -r requirements.txt
~/.venvs/synthesis_helper/bin/python scripts/fetch_cytoscape.py   # offline asset
```

> **macOS + iCloud note.** If the project lives under `~/Documents` or
> `~/Desktop` with iCloud "Optimize Mac Storage" on, files inside `.venv/`
> can get evicted to cloud stubs and Python imports stall for minutes
> during the MCP handshake. Keep the venv outside any synced folder
> (`~/.venvs/...` is fine).

Register with Claude Code (preferred — the launcher script materializes any
evicted source files and finds the venv automatically):

```bash
claude mcp add synthesis-helper --scope user -- /absolute/path/to/synthesis_helper/scripts/start_mcp.sh
```

Verify with `/mcp` inside Claude Code — you should see
`synthesis-helper: connected`.

#### Environment variables

- `SYNTHESIS_DATA_DIR` — point at a different `data/` directory (default:
  the repo's own `data/`).
- `SYNTHESIS_MAX_PATHWAYS` — hard cap on `enumerate_pathways_for` (default
  `1000`).
- `SYNTHESIS_HELPER_VENV` — override the venv path used by
  `scripts/start_mcp.sh`.

#### Optional: refresh the EC-name table

`data/ec_names.tsv` ships with a full EC → enzyme-name table derived from
[ExPASy ENZYME](https://enzyme.expasy.org/) (8k+ entries, transferred
entries resolved). To rebuild:

```bash
curl -sS https://ftp.expasy.org/databases/enzyme/enzyme.dat -o /tmp/enzyme.dat
.venv/bin/python scripts/build_ec_names.py /tmp/enzyme.dat
```

The file is redistributed under the same
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) license as the
upstream ENZYME database (© SIB Swiss Institute of Bioinformatics).

### 1.4 Example prompts

- "How many steps from native E. coli metabolites to vanillin?" →
  `get_shell`.
- "List the top 3 pathways to 2′-dehydrokanamycin a, ranked by feasibility,
  and avoid P450s." → `compare_pathways` + `pathway_to_composition`.
- "If I feed PABA into the cell, what becomes newly reachable?" →
  `resynthesize_with_fed`.
- "Open an interactive cascade view for chemical 317157." →
  `open_cascade_interactive`.

---

## 2. Validating `synthesize` and the cascade traceback
*Author: Benson Li*

The hard part of any retrobiosynthesis tool is convincing yourself the
algorithm is right. Empirical correctness on a 9k-reaction MetaCyc dump is
hard to inspect by eye, so we built a layered validation strategy: a
hand-drawn synthetic cell where every shell number is computable on paper,
implementation-agnostic invariant tests, and a real-world cross-reference
against a published MetaCyc pathway.

The full test suite is **114 tests** across 10 files, run with
`.venv/bin/pytest tests/ -q`.

### 2.1 The hand-drawn small cell (`validation/small_cell.md`)

A purely synthetic 12-chemical, 11-reaction "cell" defined inline in
`tests/fixtures/small_cell.py`. It does **not** read anything under
`data/`, so changes to the real MetaCyc corpus cannot affect these tests
(and conversely, a regression in the algorithm cannot hide behind corpus
churn).

The fixture is engineered to exercise every structural scenario the BFS
must handle:

| Scenario | What it pins down |
|---|---|
| (a) Linear chain N1 → M1 → M2 → M_deep | Shells propagate correctly past depth 2. |
| (b) Multi-substrate `R3: N1 + N2 → M4` | Reaction shell = `max(substrate shell) + 1`. |
| (c) Universal cofactor in `R4: M4 + U1 → M5` | A shell-0 universal cannot inflate a downstream shell. |
| (d) Multi-producer target T (via R7 *or* R6) | Cascade contains both routes; pathway enumeration must surface both. |
| (e) Unreachable target `M7 ← X_missing` | BFS stops cleanly; traceback raises `ValueError`. |
| (f) Tied-shell producers of M3b (R8 at shell 1, R9 at shell 2) | The shorter route wins. |
| (g) Fed chemical `F1 ← R10 → M6` | `resynthesize_with_fed` enables the reaction without re-running the whole pipeline. |

Every expected shell, reaction, cascade member, and pathway is tabulated
in `validation/small_cell.md` and asserted by tests in
`tests/test_small_cell_synthesize.py`,
`tests/test_small_cell_traceback.py`, and
`tests/test_small_cell_pathways.py`.

This fixture also surfaces a real bug: in `pathways.py:78–80`,
`enumerate_pathways` only recurses on the first non-native substrate of a
multi-substrate producer reaction. The test
`test_pathways_enumerates_both_complete_routes` is decorated
`@pytest.mark.xfail(strict=True)` so the bug stays visible — when fixed,
the test flips to XPASS and the suite goes red, forcing a conscious
removal of the marker.

### 2.2 Implementation-agnostic invariants (`tests/test_invariants.py`)

The small cell tests pin specific values; the invariant tests pin the
**meaning** of `synthesize` and `traceback`. They apply to any correct
implementation, so a teammate could swap the BFS for a totally different
algorithm and these tests would still catch a spec violation.

`synthesize()` invariants asserted:

- Every native and universal sits at shell 0.
- Shells are non-negative integers; reactions start at 1, chemicals at 0.
- **BFS defining equation:** `shell(R) == max(shell(s) for s in R.substrates) + 1`.
- Every product of an enabled reaction is reachable (and likewise for
  substrates).
- Every non-native chemical has a producer at the matching shell.
- **Fixed-point completeness:** no unchosen reaction has all substrates
  already reachable — i.e., the BFS ran to closure.
- Every reachable chemical can be successfully tracebacked.

Cascade (traceback) invariants asserted:

- Cascade reactions ⊆ enabled reactions (no phantoms).
- No orphan reactions in the cascade — every reaction produces something
  that contributes to the target.
- **Closure under substrate backtracking:** every non-native substrate of
  a cascade reaction is itself produced by another cascade reaction.
- The cascade of a non-native target contains at least one reaction
  producing it.
- **Monotonicity:** `Cascade(M) ⊆ Cascade(T)` when M is an intermediate of
  T.
- Cascade of a shell-0 target is empty; cascade of an unreachable target
  raises `ValueError`.
- **Cross-algorithm:** every reaction in every enumerated pathway must
  appear in the cascade.

### 2.3 Edge cases (`tests/test_edge_cases.py`)

A second tier of inline-fixture tests that target structural failure
modes the small cell does not cover:

- Self-referential reaction `R: A → A` (matches a real MetaCyc rxn 88
  pattern) does not loop the BFS.
- A reaction producing a shell-0 chemical does not overwrite that
  chemical's native shell.
- Multi-product reactions — both products reachable.
- Deep chain (20 shells) terminates cleanly.
- Zero-substrate reactions enable at shell 1 (pins the current
  `all([])==True` behavior so any future change is conscious).
- Zero-product reactions are no-ops.
- Determinism — same inputs, same shell map every time.
- Scale smoke test — a 100-chemical synthetic corpus completes.
- Chemical-level cycle in a cascade (B ↔ C) — `enumerate_pathways` must
  not loop.
- `max_pathways` cap is enforced even on deep cascades.
- Model invariant: `Chemical` equality is by id only.
- Parser silently drops unknown substrate ids (documented risk).

### 2.4 MetaCyc smoke test (`tests/test_metacyc_smoke.py`)

A minimal regression check that runs the full pipeline on the real
MetaCyc data shipped in `data/` and asserts a small set of stable facts —
total reachable counts, a few well-known shell numbers, the lecture
example "L-tyrosine → 4-fumarylacetoacetate", etc. If anyone breaks the
parser or the BFS in a way that escapes the synthetic tests, this
catches it.

### 2.5 Real-world cross-reference: tyrosine degradation
(`validation/tyrosine_degradation_validation.md`)

A manual validation against MetaCyc pathway TYRFUMCAT-PWY ("L-tyrosine
degradation I"). We compared our shell numbers to the canonical
five-step biological pathway and got:

| Chemical | Our shell | MetaCyc step |
|---|---|---|
| L-tyrosine | 0 | (start) |
| 4-hydroxyphenylpyruvate | 1 | 1 |
| Homogentisic acid | 2 | 2 |
| 4-fumarylacetoacetate | 3 | 4 |

Two intermediate steps come back differently from canonical MetaCyc, and
both are explained by **corpus quality, not algorithm bugs**:

1. Our BFS finds an alternative (non-canonical) reaction for the
   tyrosine → 4-HPP step that does not require the 2-ketoglutarate
   cofactor (which sits at shell 1 in our model, not shell 0).
2. The intermediate 4-maleylacetoacetate is missing from
   `good_chems.txt`, which collapses MetaCyc's two-step
   dioxygenase + isomerase into a single reaction in our corpus.

The conclusion in the validation memo is explicit: the BFS is computing
the correct *minimum* given the available reactions; the discrepancies
are entirely upstream of the algorithm. This is a real example of using
domain ground truth to **localize a discrepancy** to corpus rather than
code.

### 2.6 MCP-layer tests

- `tests/test_mcp_server.py` — exercises every MCP tool / resource
  end-to-end against the real loaded hypergraph (not mocked), including
  reference resolution edge cases, structure search, fed-chemical
  resynthesis, and the interactive viewer return shape.
- `tests/test_html_render.py` — checks the generated Cytoscape JSON for
  pathways and cascades (node / edge counts, theme metadata, cofactor
  duplication mode).

### 2.7 Running everything

```bash
.venv/bin/pytest tests/ -q
```

Expected: all green, with exactly one **XFAIL** (the
`pathways.py:78–80` multi-substrate enumeration bug, kept visible on
purpose).

---

## 3. _(placeholder for the next teammate)_

> Add your section here — e.g., the BFS / synthesize implementation, the
> data-prep scripts, the composition → DNA-construct work, or anything
> else you owned. Add yourself with `*Author: <Your Name>*` under the
> heading.

---

## Run the standalone pipeline

If you just want the reachables file without going through MCP:

```bash
.venv/bin/python main.py
```

Writes `data/metacyc_L2_reachables.txt` (columns: id, name, inchi, shell).
