# Phase 35 — Enrichment Factor Calculator (EF@K)

**Version:** 1.1 | **Tier:** Micro | **Date:** 2026-03-26

## Goal
Implement the Enrichment Factor (EF@K) metric for virtual screening evaluation.
Simulate a scoring scenario, compute EF at multiple cutoffs, and visualize
the cumulative enrichment curve.

CLI: `python main.py --input data/compounds.csv --threshold 7.0`

Outputs: ef_results.csv, enrichment_curve.png, ef_bar.png

## EF@K Formula
EF@K = (hits_in_topK / K) / (total_hits / N)

Where:
- K = number of top-ranked compounds inspected
- hits_in_topK = actives in the top K
- total_hits = all actives in the dataset
- N = total compounds
- Random baseline = 1.0 at all cutoffs

## Logic
- Define "active" as pIC50 >= threshold (default 7.0)
- Score = pIC50 directly (oracle scoring) — establishes maximum possible EF
- Rank compounds by score descending
- Compute EF at K = 5, 10, 20, 25% of N (i.e., 2, 4, 9, 11 compounds from 45)
- Compute BEDROC (alpha=20) as supplementary metric — rewards early enrichment more
- Plot 1: Cumulative hits vs rank (enrichment curve) vs random baseline
- Plot 2: EF@K bar chart at 5%, 10%, 20%, 25% cutoffs

## Notes
- Phase 35 is the ML track entry point — establishes the evaluation harness used in 36–54
- EF@K function will be reused across all classifier phases
- pIC50 as oracle score = theoretical maximum EF — actual models will score lower

## Key Concepts
- Enrichment Factor (EF@K) — virtual screening evaluation metric
- Oracle scoring (pIC50 as score) to establish theoretical maximum EF
- BEDROC (alpha=20) — Boltzmann-enhanced discrimination metric for early enrichment
- Cumulative enrichment curve vs random baseline

## Verification Checklist
- [x] EF@K constant at 1.50x for all cutoffs <= 30 compounds (theoretical maximum = 1/hit_rate)
- [x] ef_results.csv and enrichment_curve.png saved to output/
- [x] 30/45 actives (66.7% hit rate) at pIC50 >= 7.0 threshold
- [x] EF bar chart shows all cutoffs at ceiling value

## Risks
- Oracle scoring is unrealistically optimistic; real models will produce lower EF values
- High hit rate (66.7%) limits the dynamic range of EF (max = 1.50x); lower hit rates would show more differentiation

## Actual Results (v1.1)

| Cutoff | K | Hits@K | EF |
|---|---|---|---|
| Top 5% | 2 | 2 | 1.50× |
| Top 10% | 4 | 4 | 1.50× |
| Top 20% | 9 | 9 | 1.50× |
| Top 25% | 11 | 11 | 1.50× |

**Actives:** 30/45 compounds (66.7%) with pIC50 ≥ 7.0
**Key insight:** EF@K is constant at 1.50× for all cutoffs ≤ 30 compounds with oracle scoring. This is the theoretical maximum: EF_max = 1/hit_rate = 45/30 = 1.50. Oracle sorting places all actives first. Real models will score lower.
**EF harness ready for Phases 36–54**
