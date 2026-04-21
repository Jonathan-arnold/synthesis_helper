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

**Optional — visualization tools.** `visualize_pathway` and `visualize_cascade`
shell out to Graphviz's `dot` binary. Install it once:

```bash
brew install graphviz       # macOS
# or: sudo apt install graphviz
```

The text tools work without Graphviz; only the two PNG tools require it.

### Register with Claude Code

Project scope (this repo only):

```bash
claude mcp add synthesis-helper --scope project -- .venv/bin/python -m synthesis_helper
```

User scope (available everywhere, but requires an absolute path to the venv):

```bash
claude mcp add synthesis-helper \
  --scope user \
  -- /absolute/path/to/synthesis_helper/.venv/bin/python -m synthesis_helper
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
| `visualize_pathway(chemical_ref, pathway_index)` | PNG diagram of one pathway (requires Graphviz). |
| `visualize_cascade(chemical_ref, max_reactions)` | PNG diagram of full cascade tree (requires Graphviz). |
| `resynthesize_with_fed(fed_chemical_refs)` | Re-run BFS with additional fed chemicals; returns delta. |

### Resources

- `synthesis://stats` — counts, load time, data dir.
- `synthesis://reachables` — TSV dump of every reachable chemical.
- `synthesis://chemical/{id}` — chemical + produced_by / consumed_by reaction ids.
- `synthesis://reaction/{id}` — reaction with expanded substrate/product names.

### Env vars

- `SYNTHESIS_DATA_DIR` — point at a different `data/` directory (default: repo `data/`).
- `SYNTHESIS_MAX_PATHWAYS` — hard cap on `enumerate_pathways_for` (default 1000).

### Example prompts (once registered in Claude Code)

- “從 native E. coli 到 vanillin 走幾步？” → `get_shell`.
- “列出通往 2′-dehydrokanamycin a 的前 3 條 pathway。” → `enumerate_pathways_for` + `describe_pathway`.
- “若我餵 PABA 進細胞，多出多少 reachable？” → `resynthesize_with_fed`.
- “給我 target 317157 的 cascade 摘要。” → `get_cascade`.
- “畫出 target 317157 的第一條 pathway。” → `visualize_pathway` (PNG inline).

## Run the standalone pipeline

```bash
.venv/bin/python main.py
```

Writes `data/metacyc_L2_reachables.txt`.

## Tests

```bash
.venv/bin/pytest tests/ -q
```
