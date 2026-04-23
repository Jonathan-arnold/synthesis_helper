# synthesis_helper

A tool that builds a synthesis hypergraph from reaction data and enumerates all
cascades, pathways, and compositions for a reachable target.

## Run as a Claude Code MCP server

This project ships an MCP (Model Context Protocol) server so Claude Code can
call `synthesize`, `traceback`, `get_cascade`, `enumerate_pathways`, and
`resynthesize_with_fed` directly. The hypergraph is built lazily on first
call and cached in-process, so subsequent queries return in milliseconds.

### Install

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.txt
```

> **macOS + iCloud note.** If this repo lives under `~/Documents` or `~/Desktop`
> with iCloud "Optimize Mac Storage" on, `.venv/` files get evicted to cloud
> stubs (`dataless` in `ls -laO`), and Python imports stall for minutes while
> the MCP handshake is blocked. Put the venv **outside** any synced folder:
>
> ```bash
> python3 -m venv ~/.venvs/synthesis_helper
> ~/.venvs/synthesis_helper/bin/pip install -r requirements.txt
> ```
>
> Then register using `scripts/start_mcp.sh` (see below), which also warms any
> evicted project source via `brctl download` before launching.

**Enzyme names on reaction nodes.** `data/ec_names.tsv` is shipped with a
full EC → name table derived from [ExPASy ENZYME](https://enzyme.expasy.org/)
(8k+ entries, transferred entries resolved to their new ECs). The
interactive viewer uses the enzyme name as the node title and the EC
number as a small subtitle. To refresh the table:

```bash
curl -sS https://ftp.expasy.org/databases/enzyme/enzyme.dat -o /tmp/enzyme.dat
.venv/bin/python scripts/build_ec_names.py /tmp/enzyme.dat
```

The `ec_names.tsv` file is distributed under the same
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) license as the
upstream ENZYME database (© SIB Swiss Institute of Bioinformatics).
Delete the file to fall back to EC-only labels.

### Register with Claude Code

Preferred — launcher script (handles iCloud materialization and venv
resolution automatically):

```bash
claude mcp add synthesis-helper --scope user -- /absolute/path/to/synthesis_helper/scripts/start_mcp.sh
```

The launcher honors `SYNTHESIS_HELPER_VENV` if you keep the venv somewhere
other than `~/.venvs/synthesis_helper`.

Raw invocation (only safe when the venv is guaranteed to be outside iCloud):

```bash
claude mcp add synthesis-helper --scope user -- ~/.venvs/synthesis_helper/bin/python -m synthesis_helper
```

Verify with `/mcp` inside Claude Code — you should see `synthesis-helper: connected`.

### Tools

| Name | Purpose |
|---|---|
| `find_chemical(query, limit)` | Search by id / name / InChI. |
| `get_shell(chemical_ref)` | Min reaction steps from native metabolites. |
| `list_reachables(max_shell, name_contains, limit)` | Filter the reachables table. |
| `get_cascade(chemical_ref, max_producers_per_chemical)` | Full cascade (tree of reactions + leaf natives). |
| `enumerate_pathways_for(chemical_ref, max_pathways, …)` | Individual linear pathways. |
| `describe_pathway(chemical_ref, pathway_index)` | Markdown-rendered pathway. |
| `open_pathway_interactive(chemical_ref, pathway_index)` | Render an **interactive HTML page** (Cytoscape.js, offline) and open it in the default browser. Draggable nodes, 5 layout algorithms, group-by themes (shell / ec_class / role), cofactor toggle (inline / shared / hide), click-to-focus. Double-click any node to copy its info to the clipboard. |
| `open_cascade_interactive(chemical_ref, max_reactions)` | Same as above for the full cascade, with themes {shell, ec_class, role, producer_depth}. |
| `resynthesize_with_fed(fed_chemical_refs)` | Re-run BFS with additional fed chemicals; returns delta. |

#### Interactive HTML viewer

`open_pathway_interactive` and `open_cascade_interactive` emit a
**self-contained HTML file** (Cytoscape.js inlined — works offline) and
launch it in the OS default browser. The page supports:

- **5 layout algorithms** in the toolbar: Layered (shell-based preset),
  Breadth-first, Concentric, Force (cose), Grid.
- **Group-by themes** — independent from the layout. Pathway views offer
  {shell, ec_class, role}; cascade views offer {shell, ec_class, role,
  producer_depth}. Switching theme re-groups nodes into different
  dashed-box containers, recolors them from the theme's palette, and
  re-runs the active layout — all client-side, no server round-trip.
- **Cofactor mode** — Inline (default) duplicates currency cofactors
  (H2O, ATP, NAD(P)H, …) next to each reaction so no single node acts as
  a hub; Shared reverts to the classical shared-node graph; Hide removes
  cofactors entirely.
- **Click-to-focus** — single-click any node to dim everything except
  that node and its 1-hop neighborhood. Click empty space or press ESC
  to clear.
- **Drag** nodes around; hover for a full tooltip (name, InChI, substrates,
  products).
- **Double-click** any chemical or reaction node to copy its structured
  details to the clipboard, formatted so it pastes cleanly back into a
  Claude conversation ("Chemical #317157 …", "Reaction #42 — EC 1.1.1.1 …").

**First-time setup** (only needed once per clone): download the inlined
Cytoscape.js asset.

```bash
.venv/bin/python scripts/fetch_cytoscape.py
```

This writes `synthesis_helper/mcp/assets/cytoscape.min.js` (pinned
version, SHA-256 verified). Cytoscape.js is MIT-licensed.

### Resources

- `synthesis://stats` — counts, load time, data dir.
- `synthesis://reachables` — TSV dump of every reachable chemical.
- `synthesis://chemical/{id}` — chemical + produced_by / consumed_by reaction ids.
- `synthesis://reaction/{id}` — reaction with expanded substrate/product names.

### Env vars

- `SYNTHESIS_DATA_DIR` — point at a different `data/` directory (default: repo `data/`).
- `SYNTHESIS_MAX_PATHWAYS` — hard cap on `enumerate_pathways_for` (default 1000).

### Example prompts (once registered in Claude Code)

- "How many steps does it take to go from native E. coli metabolites to vanillin?" → `get_shell`.
- "List the top 3 pathways to 2′-dehydrokanamycin a." → `enumerate_pathways_for` + `describe_pathway`.
- "If I feed PABA into the cell, how many more chemicals become reachable?" → `resynthesize_with_fed`.
- "Give me a cascade summary for target 317157." → `get_cascade`.
- "Open an interactive pathway view for 317157." → `open_pathway_interactive` (Cytoscape.js in the browser; draggable nodes, group-by theme switcher, cofactor toggle, click-to-focus, double-click-to-copy).

## Run the standalone pipeline

```bash
.venv/bin/python main.py
```

Writes `data/metacyc_L2_reachables.txt`.

## Tests

```bash
.venv/bin/pytest tests/ -q
```
