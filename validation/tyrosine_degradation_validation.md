# Validation Record: Tyrosine Degradation Pathway (TYRFUMCAT-PWY)

**Date:** 2026-04-19  
**Validator:** Manual cross-reference with MetaCyc BioCyc API  
**Purpose:** Verify that the BFS shell assignment algorithm produces correct minimum step counts for a known biological pathway.

---

## Ground Truth Source

**MetaCyc pathway:** TYRFUMCAT-PWY — "L-tyrosine degradation I"  
**URL:** `https://websvc.biocyc.org/getxml?id=META:TYRFUMCAT-PWY&detail=full`  
**Starting compound:** L-tyrosine (shell 0 — universal metabolite)

### MetaCyc canonical reaction sequence

| Step | MetaCyc reaction ID | EC | Substrate(s) | Product(s) |
|---|---|---|---|---|
| 1 | TYROSINE-AMINOTRANSFERASE-RXN | 2.6.1.5 | L-tyrosine + 2-ketoglutarate | 4-hydroxyphenylpyruvate + L-glutamate |
| 2 | 4-HYDROXYPHENYLPYRUVATE-DIOXYGENASE-RXN | 1.13.11.27 | 4-HPP + O₂ | homogentisate + CO₂ |
| 3 | HOMOGENTISATE-12-DIOXYGENASE-RXN | 1.13.11.5 | homogentisate + O₂ | **4-maleylacetoacetate** + H⁺ |
| 4 | MALEYLACETOACETATE-ISOMERASE-RXN | 5.2.1.2 | 4-maleylacetoacetate | **4-fumarylacetoacetate** |
| 5 | FUMARYLACETOACETASE-RXN | — | 4-fumarylacetoacetate | fumarate + 3-ketobutyrate |

Expected minimum steps from L-tyrosine to 4-fumarylacetoacetate: **4**.

---

## Our Model's Output

| Chemical | Our shell | MetaCyc step |
|---|---|---|
| L-tyrosine | 0 | — (starting compound) |
| 4-hydroxyphenylpyruvate | **1** | 1 |
| Homogentisic acid | **2** | 2 |
| 4-fumarylacetoacetate | **3** | 4 |

Our model gives shell **3** for 4-fumarylacetoacetate; MetaCyc expects **4** steps.

---

## Step-by-step Comparison

### Step 2 — MATCH ✓

| | Our model (rxn 358773) | MetaCyc |
|---|---|---|
| EC | 1.13.11.27 | 1.13.11.27 |
| Substrates | 4-HPP + O₂ | 4-HPP + O₂ |
| Products | homogentisate + CO₂ | homogentisate + CO₂ |

Exact match on EC number, substrates, and products.

---

### Step 1 — DISCREPANCY (different reaction used)

| | Our model (rxn 308800) | MetaCyc canonical |
|---|---|---|
| EC | *(none)* | 2.6.1.5 |
| Substrates | L-tyrosine only | L-tyrosine + **2-ketoglutarate** |
| Products | 4-HPP + NH₃ | 4-HPP + **L-glutamate** |

**Root cause:** 2-ketoglutarate (alpha-ketoglutarate) is at shell **1** in our model, not shell 0. The canonical aminotransferase reaction requires it, which would place 4-HPP at shell 2. Our BFS found an alternative non-canonical reaction (rxn 308800) that deaminates L-tyrosine without a cofactor, placing 4-HPP at shell 1 instead.

**Implication for algorithm:** The BFS is working correctly — it found the *minimum* path given the available reactions. The discrepancy is a corpus issue (non-canonical reaction available + 2-ketoglutarate not at shell 0), not an algorithm bug.

---

### Step 3–4 — DISCREPANCY (missing intermediate)

| | Our model (rxn 1853) | MetaCyc |
|---|---|---|
| EC | 1.13.11.5 | 1.13.11.5 (step 3) + 5.2.1.2 (step 4) |
| Substrates | homogentisate + O₂ | homogentisate + O₂ |
| Products | **4-fumarylacetoacetate** | **4-maleylacetoacetate** → then isomerized to 4-fumarylacetoacetate |

**Root cause:** 4-maleylacetoacetate does not exist in `good_chems.txt`. Without the intermediate chemical, MetaCyc's two-step sequence (dioxygenase → isomerase) collapses into a single reaction in our corpus. Our rxn 1853 is labeled EC 1.13.11.5 but produces the isomerized product directly.

**Implication for algorithm:** Again, the BFS is correct — it computes the minimum given available reactions. The discrepancy is a corpus issue: the maleylacetoacetate intermediate is absent, collapsing two steps into one.

---

## Conclusion

| Claim | Result |
|---|---|
| BFS logic correctly assigns minimum shell from available reactions | ✓ CONFIRMED |
| Shell numbers match MetaCyc when corpus is complete | ✗ DIFFERS by 1 step for 4-fumarylacetoacetate |
| Discrepancy source is algorithm bug | ✗ NOT the case |
| Discrepancy source is corpus quality (missing intermediates, alternative reactions) | ✓ CONFIRMED |

**Verdict:** The BFS algorithm is correct. The shell-number differences from MetaCyc ground truth are explained entirely by two corpus-level issues:
1. Some intermediates are absent from `good_chems.txt`, collapsing multi-step sequences.
2. Alternative (non-canonical) reactions exist in the corpus that provide shorter routes than the canonical biological pathway.

---

## How to Run This Validation

```python
from synthesis_helper.parser import parse_chemicals, parse_metabolite_list, parse_reactions
from synthesis_helper.synthesize import synthesize
from pathlib import Path

DATA = Path("data")
chemicals = parse_chemicals(DATA / "good_chems.txt")
reactions = parse_reactions(DATA / "good_reactions.txt", chemicals)
native = parse_metabolite_list(DATA / "minimal_metabolites.txt", chemicals)
universal = parse_metabolite_list(DATA / "ubiquitous_metabolites.txt", chemicals)
hg = synthesize(reactions, native, universal, verbose=False)

ids = {
    "l-tyrosine":             122,
    "4-hydroxyphenylpyruvate": 3255,
    "homogentisic acid":       3202,
    "4-fumarylacetoacetate":   3912,
}
for name, cid in ids.items():
    c = chemicals[cid]
    print(f"{name}: shell {hg.chemical_to_shell.get(c, 'UNREACHABLE')}")
# Expected: 0, 1, 2, 3
```
