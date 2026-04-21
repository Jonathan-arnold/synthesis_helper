# BioE234 Synthesis Project

## Overview

This project implements a **bottom-up retrobiosynthesis / pathway synthesis tool** inspired by the Act Synthesizer (UCB & 20n). Given a source strain's native metabolites and a corpus of known enzymatic reactions, the tool simulates a hypothetical "metagenome" cell, enumerates all reachable chemicals, traces cascades back to native metabolites, extracts discrete pathways, and converts those pathways into enzyme compositions that could be realized as DNA constructs.

The approach is **synthesis** (forward expansion from starting materials) rather than classical retrobiosynthesis (backward from target). We specify starting metabolites, run a breadth-first waveform expansion across the reaction corpus, and discover what can be produced.

## Core Concepts

### Native vs. Universal Metabolites
- **Native metabolites** (`minimal_metabolites.txt`): the starting pool — what exists in the hypothetical cell before any enzymes are added (glucose, ions, H2O, CO2, phosphate, etc.), plus any **fed chemicals** the user supplies (e.g., PABA).
- **Universal metabolites** (`universal_metabolites.txt`): cofactors, amino acids, and nucleotides that the target organism is assumed to produce (NADPH, ATP, Ala/Leu/Trp/Phe/Ser, A/T/C/G, etc.). These get merged into shell 0 so the expansion doesn't have to rediscover them from glucose.

### Shells
Each chemical and reaction is labeled with the **shell** in which it first becomes enabled. Shell 0 is the starting pool. Shell N+1 contains everything first reachable after one more round of enzymatic expansion. Shell number doubles as the minimum step count from native metabolites.

### HyperGraph
The primary output of the expansion:
```python
@dataclass
class HyperGraph:
    chemical_to_shell: Dict[Chemical, int]
    reaction_to_shell: Dict[Reaction, int]
```

### Cascade
The tree of all reactions that can produce a given chemical, recursively — a Cascade per Reachable is what lets us enumerate pathways later:
```
Cascade: Set<Reaction> rxnsThatFormPdt
HyperGraph also holds: Map<Chemical, Cascade> chemicalToCascade
```

### Pathway
A single directed graph (nodes = metabolites, edges = reactions) representing **one specific route** by which one equivalent of the target is made from stoichiometric starting material. A Cascade contains many Pathways mixed together.

### Composition
An ordered list of enzymes (e.g., `[AfGGPPS, OsCPS4, OsKSL4]`) derived from a Pathway — one enzyme per reaction. This is the handoff point to a DNA construct (promoter → CDS → CDS → CDS → terminator).

## Algorithm (Waveform Expansion / BFS)

1. **Initialize**: `chemical_to_shell` = {m: 0 for m in native ∪ universal ∪ fed}; `reaction_to_shell` = {}.
2. **Expand one shell**: for every reaction in the corpus not yet enabled, check whether all substrates are already in `chemical_to_shell`. If yes, add the reaction at the current shell number and add any new products at the current shell number.
3. **Repeat** until a full pass adds no new reactions or chemicals (fixed point).
4. **Output**: populated `HyperGraph` → dump to `metacyc_L2_reachables.txt` (columns: id, name, inchi, shell).

### Traceback
After expansion, for any specified target chemical:
- Look up the chemical in the HyperGraph.
- Walk backward: for each reaction producing it, recurse on that reaction's substrates until native metabolites (shell 0) are hit.
- Build the `Cascade` — the full tree of routes.

### Pathway Enumeration (stubbed in source slides, to be implemented)
```python
def enumerate_pathways(cascade: Any) -> list[Any]:
    """Enumerate all pathways implied by a cascade."""
```
Walk the Cascade and flatten it into every valid singular route from native metabolites to the target. Memory practicality constraints apply — the cascade can branch combinatorially.

### Pathway → Composition
For each Reaction in a Pathway, assign an enzyme (one that catalyzes that reaction — EC-number lookup or curated mapping). The ordered list of enzymes is the Composition.

## Data Model

```python
class Chemical:
    id: long
    inchi: str
    smiles: str
    name: str

class Reaction:
    id: long
    substrates: set[Chemical]
    products: set[Chemical]
    # optional: ecnum, observation flag
```

### `isObserved` filter (currently skipped)
If non-MetaCyc sources are added (e.g., BRENDA), an `Observation` table with an `isObserved` boolean must gate which reactions enter the corpus. All MetaCyc reactions are treated as observed=true in the slides' implementation.

## Input Files

| File | Contents |
|---|---|
| `metacyc_chemicals.txt` | id, name, inchi (one chemical per row) |
| `metacyc_reactions.txt` | rxnid, ecnum, space-separated substrate ids, space-separated product ids |
| `minimal_metabolites.txt` | native metabolites (starting pool) |
| `universal_metabolites.txt` | cofactors/AAs/NTs to merge into shell 0 |

Optional fed chemicals can be appended to the native set to simulate biodegradation/biotransformation of exogenous substrates.

## Output

| File | Contents |
|---|---|
| `metacyc_L2_reachables.txt` | id, name, inchi, shell — the full Reachables list |

From a target chemical in this file, the user can query the Cascade, enumerate Pathways, and derive Compositions.

## Public API

```python
def synthesize(
    all_reactions: list[Reaction],
    native_metabolites: set[Chemical],
    universal_metabolites: set[Chemical],
    verbose: bool = True,
) -> HyperGraph: ...

def traceback(hypergraph: HyperGraph, target: Chemical) -> Cascade: ...

def enumerate_pathways(cascade: Cascade) -> list[Pathway]: ...

def pathway_to_composition(pathway: Pathway, enzyme_map: dict) -> list[Enzyme]: ...
```

## MCP surface (`synthesis_helper.mcp`)

The `synthesis_helper.mcp` subpackage wraps the API above as a FastMCP server
(stdio transport) so Claude Code can call it directly. The baseline
hypergraph is lazy-loaded and cached in-process; see `README.md` for registration.

| Entry point | Exposed |
|---|---|
| Tools | `find_chemical`, `get_shell`, `list_reachables`, `get_cascade`, `enumerate_pathways_for`, `describe_pathway`, `resynthesize_with_fed` |
| Resources | `synthesis://stats`, `synthesis://reachables`, `synthesis://chemical/{id}`, `synthesis://reaction/{id}` |

Run locally: `python -m synthesis_helper`. Env: `SYNTHESIS_DATA_DIR`, `SYNTHESIS_MAX_PATHWAYS`.

## Example Query

> "What is the minimum number of reaction steps to convert native E. coli metabolites to 2′-dehydrokanamycin a?"

Answer: `hypergraph.chemical_to_shell[2'-dehydrokanamycin a]` → `8` (per the lecture example).

## Project Scope Checklist

- [ ] Parse `metacyc_chemicals.txt` and `metacyc_reactions.txt` into `Chemical` and `Reaction` objects.
- [ ] Load `minimal_metabolites.txt` + `universal_metabolites.txt` + optional fed chemicals.
- [ ] Implement `synthesize()` — BFS waveform expansion producing a HyperGraph.
- [ ] Dump reachables to `metacyc_L2_reachables.txt`.
- [ ] Implement `traceback()` — build a Cascade for any target chemical.
- [ ] Implement `enumerate_pathways()` — flatten a Cascade into individual Pathways.
- [ ] Implement `pathway_to_composition()` — map reactions to enzymes.
- [ ] (Stretch) Visualize Cascades as directed graphs (see cannabidiol example on slide 37).
- [ ] (Stretch) Export Compositions to DNA construct specs (promoter + CDSs + terminator).

## References in the Lecture

- Slides 7, 8–10: BFS waveform expansion visual.
- Slides 11–17: shell-by-shell state evolution of the HyperGraph.
- Slide 21: `HyperGraph` dataclass definition.
- Slide 22: `synthesize()` signature.
- Slides 23–25: file formats.
- Slide 27: example reachables output.
- Slides 30–38: Cascade / Pathway / Traceback concepts.
- Slide 40: `enumerate_pathways()` signature (not implemented in starter).
- Slides 41–42: Pathway → Composition → DNA construct.

## Notes

- This tool produces a **design objective**, not a production-ready strain. Concretizing each reaction with a specific enzyme (and dealing with promiscuity, kinetics, toxicity, balance) is a downstream concern.
- The Java reference implementation lives at the course GitHub classroom repo and depends on ChemAxon jars for InChI/SMILES handling. A Python port is acceptable and aligns with the `synthesize()` signature shown.
