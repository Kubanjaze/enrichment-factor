# enrichment-factor — Phase 35

Enrichment Factor (EF@K) calculator for virtual screening evaluation.
Establishes the EF harness used throughout the ML track (Phases 35–54).

## Usage

```bash
PYTHONUTF8=1 python main.py --input data/compounds.csv --threshold 7.0
```

## Outputs

| File | Description |
|---|---|
| `output/ef_results.csv` | EF at 5%, 10%, 20%, 25% cutoffs |
| `output/enrichment_curve.png` | Cumulative actives vs rank vs random/ideal |
| `output/ef_bar.png` | EF@K bar chart |

## Formula

```
EF@K = (hits_in_topK / K) / (total_hits / N)
```

Random baseline = 1.0 at all cutoffs. Oracle scoring (pIC50) = theoretical maximum EF.
