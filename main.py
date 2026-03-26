import sys
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import argparse, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FAMILY_COLORS = {"benz": "#4C72B0", "naph": "#DD8452", "ind": "#55A868",
                 "quin": "#C44E52", "pyr": "#8172B2", "bzim": "#937860", "other": "#808080"}

def load_compounds(path):
    df = pd.read_csv(path)
    records = []
    for _, row in df.iterrows():
        name = str(row["compound_name"])
        fam = name.split("_")[0]
        try:
            pic50 = float(row["pic50"])
        except (KeyError, ValueError):
            continue
        if np.isnan(pic50):
            continue
        records.append({
            "compound_name": name,
            "family": fam if fam in FAMILY_COLORS else "other",
            "pic50": pic50,
        })
    print(f"  {len(records)} compounds with pIC50")
    return pd.DataFrame(records)

def compute_ef(y_true, y_score, k_frac):
    """
    Compute Enrichment Factor at fraction k_frac of the dataset.

    Args:
        y_true: array-like of 0/1 activity labels
        y_score: array-like of scores (higher = more likely active)
        k_frac: fraction of dataset to inspect (e.g. 0.1 = top 10%)

    Returns:
        dict with ef, k, hits_topk, total_hits, n
    """
    y_true = np.array(y_true)
    y_score = np.array(y_score)
    n = len(y_true)
    k = max(1, int(np.round(n * k_frac)))
    total_hits = y_true.sum()
    if total_hits == 0:
        return {"ef": float("nan"), "k": k, "hits_topk": 0, "total_hits": 0, "n": n, "k_frac": k_frac}
    order = np.argsort(y_score)[::-1]
    hits_topk = y_true[order[:k]].sum()
    ef = (hits_topk / k) / (total_hits / n)
    return {"ef": round(float(ef), 4), "k": k, "hits_topk": int(hits_topk),
            "total_hits": int(total_hits), "n": n, "k_frac": k_frac}

def compute_enrichment_curve(y_true, y_score):
    """Return cumulative hit counts as we go down the ranked list."""
    y_true = np.array(y_true)
    y_score = np.array(y_score)
    order = np.argsort(y_score)[::-1]
    cum_hits = np.cumsum(y_true[order])
    return cum_hits

def plot_enrichment_curve(y_true, y_score, family_labels, output_path):
    n = len(y_true)
    total_hits = sum(y_true)
    cum_hits = compute_enrichment_curve(y_true, y_score)
    random_line = np.linspace(0, total_hits, n + 1)[1:]
    ideal_line = np.minimum(np.arange(1, n + 1), total_hits)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(range(1, n + 1), cum_hits, color="#4C72B0", lw=2.5, label="Oracle score (pIC50)")
    ax.plot(range(1, n + 1), random_line, color="#808080", lw=1.5, linestyle="--", label="Random")
    ax.plot(range(1, n + 1), ideal_line, color="#55A868", lw=1.5, linestyle=":", label="Ideal")
    ax.set_xlabel("Number of compounds inspected", fontsize=11)
    ax.set_ylabel("Cumulative actives found", fontsize=11)
    ax.set_title("Enrichment Curve (Oracle pIC50 Scoring)", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()

def plot_ef_bars(ef_results, output_path):
    fracs = [r["k_frac"] for r in ef_results]
    efs = [r["ef"] for r in ef_results]
    labels = [f"Top {int(f*100)}%\n(K={r['k']})" for f, r in zip(fracs, ef_results)]

    fig, ax = plt.subplots(figsize=(7, 5))
    bars = ax.bar(labels, efs, color="#4C72B0", edgecolor="white", linewidth=0.5)
    ax.axhline(1.0, color="#808080", linestyle="--", lw=1.5, label="Random baseline")
    for bar, ef in zip(bars, efs):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.05,
                f"{ef:.2f}×", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_ylabel("Enrichment Factor (EF)", fontsize=11)
    ax.set_title("EF@K — Oracle pIC50 Scoring", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", required=True)
    parser.add_argument("--threshold", type=float, default=7.0,
                        help="pIC50 threshold to define 'active'")
    parser.add_argument("--output-dir", default="output")
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"\nLoading: {args.input}")
    df = load_compounds(args.input)

    df["active"] = (df["pic50"] >= args.threshold).astype(int)
    n_active = df["active"].sum()
    print(f"  Actives (pIC50 >= {args.threshold}): {n_active}/{len(df)} ({100*n_active/len(df):.1f}%)")

    y_true = df["active"].values
    y_score = df["pic50"].values  # oracle scoring

    # Compute EF at 5%, 10%, 20%, 25%
    cutoffs = [0.05, 0.10, 0.20, 0.25]
    ef_results = [compute_ef(y_true, y_score, frac) for frac in cutoffs]

    # Save CSV
    ef_df = pd.DataFrame(ef_results)
    ef_df.to_csv(os.path.join(args.output_dir, "ef_results.csv"), index=False)
    print(f"Saved: {args.output_dir}/ef_results.csv")

    # Plots
    plot_enrichment_curve(y_true, y_score,
                          df["family"].tolist(),
                          os.path.join(args.output_dir, "enrichment_curve.png"))
    print(f"Saved: {args.output_dir}/enrichment_curve.png")

    plot_ef_bars(ef_results, os.path.join(args.output_dir, "ef_bar.png"))
    print(f"Saved: {args.output_dir}/ef_bar.png")

    print(f"\n--- EF@K results (oracle pIC50 scoring, threshold={args.threshold}) ---")
    print(f"  {'Cutoff':>8}  {'K':>4}  {'Hits@K':>7}  {'Total Hits':>10}  {'EF':>6}")
    for r in ef_results:
        print(f"  {r['k_frac']*100:>7.0f}%  {r['k']:>4}  {r['hits_topk']:>7}  {r['total_hits']:>10}  {r['ef']:>6.2f}×")
    print("\nDone.")

if __name__ == "__main__":
    main()
