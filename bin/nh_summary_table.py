#!/usr/bin/env python3
"""
nh_summary_table.py
------------------
Generate two MultiQC-compatible summary tables of NewHybrids assignments:
 1. Default summary of P0, P1, F1, F2, Bx0, Bx1, Unassigned, Total_Hybrids, and N
 2. Masked summary (optional) where samples listed in --mask are forced to Unassigned
    before computing proportions.

If a mask file is provided but none of the listed samples are present in the run,
the masked output is still written as a valid placeholder table with a single row:
"No samples masked".
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


def read_nh_results(path):
    cats = ["P0", "P1", "F1", "F2", "Bx0", "Bx1"]
    names = ["Index", "Individual"] + cats
    return pd.read_csv(path, sep=r"\s+", skiprows=1, header=None, names=names)


def load_maps(nh_map, popmap, speciesmap):
    df_map = (
        pd.read_csv(nh_map, sep="\t", header=0)
        .rename(columns={"Sample": "Individual"})
    )
    df_pop = pd.read_csv(popmap, sep="\t", header=None, names=["Individual", "Population"])
    df_spc = pd.read_csv(speciesmap, sep="\t", header=None, names=["Individual", "Species"])
    return df_map, df_pop, df_spc


def parse_html_header(path):
    meta = {}
    with open(path) as fh:
        for line in fh:
            txt = line.strip().lstrip("<!--").rstrip("-->").strip()
            m = re.match(r'([A-Za-z0-9_]+):\s*"(.*)"', txt)
            if m:
                meta[m.group(1)] = m.group(2)
    return meta


def write_mqc_json(df, metadata, output):
    data = {
        row["Group"]: {k: v for k, v in row.items() if k not in ("Group", "N")}
        for _, row in df.iterrows()
    }
    pconfig = {
        "id": metadata.get("id"),
        "ylab": "Proportion",
        "xlab": "Group",
        "xDecimals": False,
        "min": 0.0,
        "max": 1.0,
        "scale": "YlGnBu",
    }
    out = {"data": data, "pconfig": pconfig}
    out.update({k: v for k, v in metadata.items() if k != "id"})
    with open(output, "w") as f:
        json.dump(out, f, indent=2)


def compute_summary(df, cats, hybrids):
    all_cats = cats + ["Unassigned"]

    sp_counts = df.groupby("Species")["AssignedCategory"].value_counts().unstack(fill_value=0)
    sp_n = sp_counts.sum(axis=1)
    sp_prop = sp_counts.div(sp_n, axis=0)
    for c in all_cats:
        if c not in sp_prop.columns:
            sp_prop[c] = 0.0
    sp_prop["Total_Hybrids"] = sp_prop[hybrids].sum(axis=1)
    sp_prop = sp_prop.reset_index().rename(columns={"Species": "Group"})
    sp_prop["N"] = sp_n.values.astype(int)
    sp_prop = sp_prop[["Group", "N"] + cats + ["Unassigned", "Total_Hybrids"]]

    pop_counts = (
        df.groupby(["Species", "Population"])["AssignedCategory"]
        .value_counts()
        .unstack(fill_value=0)
    )
    pop_n = pop_counts.sum(axis=1)
    pop_prop = pop_counts.div(pop_n, axis=0)
    for c in all_cats:
        if c not in pop_prop.columns:
            pop_prop[c] = 0.0
    pop_prop["Total_Hybrids"] = pop_prop[hybrids].sum(axis=1)
    pop_prop = pop_prop.reset_index()
    pop_prop["Group"] = pop_prop["Species"] + "|" + pop_prop["Population"]
    pop_n_df = pop_n.reset_index(name="N")
    pop_prop = pop_prop.merge(pop_n_df, on=["Species", "Population"])
    pop_prop["N"] = pop_prop["N"].astype(int)
    pop_prop = pop_prop[["Group", "N"] + cats + ["Unassigned", "Total_Hybrids"]]

    summary = pd.concat([sp_prop, pop_prop], ignore_index=True)
    summary[cats + ["Unassigned", "Total_Hybrids"]] = summary[
        cats + ["Unassigned", "Total_Hybrids"]
    ].round(2)
    return summary


def make_empty_masked_summary(cats):
    row = {"Group": "No samples masked", "N": 0}
    for c in cats:
        row[c] = 0.0
    row["Unassigned"] = 1.0
    row["Total_Hybrids"] = 0.0
    return pd.DataFrame([row], columns=["Group", "N"] + cats + ["Unassigned", "Total_Hybrids"])


def write_summary(df, template, outpath):
    if template:
        meta = parse_html_header(template)
        write_mqc_json(df, meta, outpath)
    else:
        df.to_csv(outpath, sep="\t", index=False, float_format="%.2f")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--result", required=True, help="NH posterior file")
    p.add_argument("--result_map", required=True, help="Index→sample TSV")
    p.add_argument("--popmap", required=True, help="Sample→population TSV")
    p.add_argument("--speciesmap", required=True, help="Sample→species TSV")
    p.add_argument("--template", help="MultiQC HTML header for default table")
    p.add_argument("--out", required=True, help="Output path (tsv or JSON) for default table")
    p.add_argument(
        "--threshold",
        type=float,
        default=0.0,
        help="Minimum posterior probability; else Unassigned",
    )
    p.add_argument("--list", help="Optional output path for list of hybrids")
    p.add_argument("--mask", help="File listing individuals to mask (one per line)")
    p.add_argument("--masked_template", help="MultiQC HTML header for masked table")
    p.add_argument("--out_mask", help="Output path (tsv or JSON) for masked table")
    p.add_argument("--list_masked", help="Optional output path for list of hybrids after masking")
    args = p.parse_args()

    df_nh = read_nh_results(args.result)
    df_map, df_pop, df_spc = load_maps(args.result_map, args.popmap, args.speciesmap)

    df = (
        df_nh.drop(columns="Individual")
        .merge(df_map, on="Index")
        .merge(df_pop, on="Individual", how="left")
        .merge(df_spc, on="Individual", how="left")
    )
    df = df.loc[:, ~df.columns.duplicated()]

    cats = ["P0", "P1", "F1", "F2", "Bx0", "Bx1"]
    hybrids = ["F1", "F2", "Bx0", "Bx1"]
    cat2idx = {c: i for i, c in enumerate(cats)}

    df["MaxProb"] = df[cats].max(axis=1)
    df["AssignedCategory"] = df[cats].idxmax(axis=1)
    df.loc[df["MaxProb"] <= args.threshold, "AssignedCategory"] = "Unassigned"

    if args.list:
        is_hybrid = df["AssignedCategory"].isin(hybrids)
        hybrid_df = df.loc[is_hybrid].copy()
        idx = hybrid_df["AssignedCategory"].map(cat2idx).to_numpy()
        vals = hybrid_df[cats].to_numpy()
        hybrid_df["Prob"] = vals[np.arange(len(hybrid_df)), idx]
        hybrid_df["Individual"].to_csv(args.list, index=False, header=False)
        print(f"✅ Hybrid list → {args.list}")

    summary = compute_summary(df, cats, hybrids)
    write_summary(summary, args.template, args.out)
    print(f"✅ Default summary table → {args.out}")

    if args.mask:
        out_m = args.out_mask or str(
            Path(args.out).with_name(Path(args.out).stem + ".masked" + Path(args.out).suffix)
        )

        mask_ids = set(Path(args.mask).read_text().split())
        present_mask_ids = mask_ids & set(df["Individual"])

        if present_mask_ids:
            df_masked = df.copy()
            df_masked.loc[df_masked["Individual"].isin(present_mask_ids), "AssignedCategory"] = "Unassigned"
            summary_m = compute_summary(df_masked, cats, hybrids)
            write_summary(summary_m, args.masked_template, out_m)
            print(f"✅ Masked summary table → {out_m}")

            if args.list_masked:
                is_hybrid_m = df_masked["AssignedCategory"].isin(hybrids)
                hybrid_m = df_masked.loc[is_hybrid_m].copy()
                idx_m = hybrid_m["AssignedCategory"].map(cat2idx).to_numpy()
                vals_m = hybrid_m[cats].to_numpy()
                hybrid_m["Prob"] = vals_m[np.arange(len(hybrid_m)), idx_m]
                hybrid_m["Individual"].to_csv(args.list_masked, index=False, header=False)
                print(f"✅ Masked hybrid list → {args.list_masked}")
        else:
            summary_m = make_empty_masked_summary(cats)
            write_summary(summary_m, args.masked_template, out_m)
            print(f"✅ No samples from mask present in this run; wrote placeholder masked summary → {out_m}")

            if args.list_masked:
                Path(args.list_masked).write_text("")
                print(f"✅ Wrote empty masked hybrid list → {args.list_masked}")


if __name__ == "__main__":
    main()
