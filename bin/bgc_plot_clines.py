#!/usr/bin/env python3

import re
import math
import argparse
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.io as pio


def parse_template(template_file):
    with open(template_file, "r") as f:
        content = f.read()
    match = re.search(r"<!--(.*?)-->", content, re.DOTALL)
    if not match:
        raise ValueError(f"No metadata block found in template: {template_file}")
    raw_block = match.group(1)
    meta = {}
    for line in raw_block.strip().splitlines():
        if ":" in line:
            key, value = line.strip().split(":", 1)
            meta[key.strip()] = value.strip().strip("\"'")
    return meta


def build_comment(meta_dict):
    lines = ["<!--"]
    for key, value in meta_dict.items():
        lines.append(f'{key}: "{value}"')
    lines.append("-->")
    return "\n".join(lines)


def write_html(fig, output_file, metadata):
    header_comment = build_comment(metadata)
    html = pio.to_html(fig, full_html=False, include_plotlyjs="cdn")
    with open(output_file, "w") as f:
        f.write(header_comment + "\n" + html)


def read_locus_order(locus_file):
    with open(locus_file, "r") as f:
        loci = [line.strip() for line in f if line.strip()]
    if not loci:
        raise ValueError(f"No locus names found in locus order file: {locus_file}")
    return pd.DataFrame({"row_id": np.arange(1, len(loci) + 1, dtype=int), "locus": loci})


def read_param_table(path, param_name):
    """
    Read posterior summary TSV and automatically detect:
      - median column: '50%' (required)
      - lower CI column: smallest percentile column other than 50%
      - upper CI column: largest percentile column other than 50%

    Supports e.g.:
      50%, 5%, 95%
      50%, 2.5%, 97.5%
    """
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df.reset_index()
    df = df.rename(columns={df.columns[0]: "row_id"})
    df["row_id"] = pd.to_numeric(df["row_id"], errors="raise")

    def pct_value(col):
        m = re.match(r"^\s*([0-9]*\.?[0-9]+)\s*%\s*$", str(col))
        return float(m.group(1)) if m else None

    pct_cols = {col: pct_value(col) for col in df.columns if pct_value(col) is not None}

    median_col = None
    for col, val in pct_cols.items():
        if abs(val - 50.0) < 1e-9:
            median_col = col
            break

    if median_col is None:
        raise ValueError(f"Required median column '50%' not found in {param_name} file: {path}")

    non_median = {col: val for col, val in pct_cols.items() if abs(val - 50.0) >= 1e-9}
    if len(non_median) < 2:
        raise ValueError(f"Could not detect lower/upper CI columns in {param_name} file: {path}")

    lower_col = min(non_median, key=non_median.get)
    upper_col = max(non_median, key=non_median.get)

    df = df.rename(
        columns={
            median_col: f"{param_name}_50",
            lower_col: f"{param_name}_low",
            upper_col: f"{param_name}_high",
        }
    )

    return df


# Basically does the CI NOT overlap with neutral expectation transition point
def classify_outlier(center_low, center_high, grad_low, grad_high):
    labels = []

    # center < 0.5 => left-shift => excess P1 ancestry
    if center_high < 0.5:
        labels.append("P1-bias")
    # center > 0.5 => right-shift => excess P0 ancestry
    elif center_low > 0.5:
        labels.append("P0-bias")

    if grad_low > 1.0:
        labels.append("Steeper")
    elif grad_high < 1.0:
        labels.append("Shallower")

    if not labels:
        return "Non-outlier"
    return ";".join(labels)


def write_summary_table(df, outfile):
    cols = [
        "row_id",
        "locus",
        "center_50",
        "center_low",
        "center_high",
        "gradient_50",
        "gradient_low",
        "gradient_high",
        "outlier_type",
        "plot_class",
    ]
    df[cols].to_csv(outfile, sep="\t", index=False, float_format="%.6f")


def canonical_plot_class(label):
    """
    Keep full multi-type combinations rather than collapsing to a single priority.
    """
    if label == "Non-outlier":
        return "Non-outlier"

    parts = [x.strip() for x in label.split(";") if x.strip()]
    order = ["P0-bias", "P1-bias", "Shallower", "Steeper"]
    parts = [x for x in order if x in parts]
    return ";".join(parts)


def color_for_class(cls):
    colmap = {
        "Non-outlier": "rgba(160,160,160,0.60)",
        "Steeper": "rgba(178,34,34,0.85)",
        "Shallower": "rgba(30,144,255,0.85)",
        "P1-bias": "rgba(255,140,0,0.90)",
        "P0-bias": "rgba(106,27,154,0.90)",
        "P0-bias;Steeper": "rgba(160,32,240,0.92)",
        "P1-bias;Steeper": "rgba(255,99,71,0.92)",
        "P0-bias;Shallower": "rgba(72,61,139,0.92)",
        "P1-bias;Shallower": "rgba(218,165,32,0.92)",
    }
    return colmap.get(cls, "rgba(0,0,0,0.65)")


def prepare_dataframe(center_file, gradient_file, locus_order_file):
    loci = read_locus_order(locus_order_file)
    center = read_param_table(center_file, "center")
    gradient = read_param_table(gradient_file, "gradient")

    df = loci.merge(center, on="row_id", how="inner").merge(gradient, on="row_id", how="inner")

    if df.shape[0] != loci.shape[0]:
        raise ValueError(
            f"Locus order file has {loci.shape[0]} loci, but merged parameter tables yielded {df.shape[0]} rows."
        )

    for col in [
        "center_50", "center_low", "center_high",
        "gradient_50", "gradient_low", "gradient_high",
    ]:
        df[col] = pd.to_numeric(df[col], errors="raise")

    eps = 1e-12
    df["center_clipped"] = df["center_50"].clip(eps, 1 - eps)
    df["u_50"] = np.log(df["center_clipped"] / (1 - df["center_clipped"])) * df["gradient_50"]
    df["log10_gradient_50"] = np.log10(np.maximum(df["gradient_50"], eps))

    df["outlier_type"] = df.apply(
        lambda r: classify_outlier(
            r["center_low"], r["center_high"],
            r["gradient_low"], r["gradient_high"],
        ),
        axis=1,
    )
    df["plot_class"] = df["outlier_type"].apply(canonical_plot_class)
    df["color"] = df["plot_class"].apply(color_for_class)

    return df


def compute_phi(h, center, gradient):
    u = math.log(center / (1.0 - center)) * gradient
    numerator = np.power(h, gradient)
    denominator = numerator + np.power(1.0 - h, gradient) * np.exp(u)
    return numerator / denominator


def rgba_with_alpha(rgba_string, alpha):
    m = re.match(r"rgba\((\d+),(\d+),(\d+),([0-9.]+)\)", rgba_string.replace(" ", ""))
    if not m:
        return rgba_string
    r, g, b, _ = m.groups()
    return f"rgba({r},{g},{b},{alpha})"


def make_hover_customdata_line(row, n):
    return np.column_stack([
        np.repeat(row["locus"], n),
        np.repeat(row["row_id"], n),
        np.repeat(row["outlier_type"], n),
        np.repeat(row["plot_class"], n),
        np.repeat(row["center_50"], n),
        np.repeat(row["center_low"], n),
        np.repeat(row["center_high"], n),
        np.repeat(row["gradient_50"], n),
        np.repeat(row["gradient_low"], n),
        np.repeat(row["gradient_high"], n),
        np.repeat(row["u_50"], n),
    ])


def cline_hovertemplate():
    return (
        "<b>%{customdata[0]}</b><br>"
        "Row: %{customdata[1]}<br>"
        "Outlier type: %{customdata[2]}<br>"
        "Plot class: %{customdata[3]}<br><br>"
        "Center (50%): %{customdata[4]:.4f}<br>"
        "Center CI: [%{customdata[5]:.4f}, %{customdata[6]:.4f}]<br>"
        "Gradient (50%): %{customdata[7]:.4f}<br>"
        "Gradient CI: [%{customdata[8]:.4f}, %{customdata[9]:.4f}]<br>"
        "u (50%): %{customdata[10]:.4f}<br><br>"
        "Hybrid index: %{x:.3f}<br>"
        "Ancestry probability: %{y:.3f}<extra></extra>"
    )


def param_hovertemplate():
    return (
        "<b>%{customdata[0]}</b><br>"
        "Row: %{customdata[1]}<br>"
        "Outlier type: %{customdata[2]}<br>"
        "Plot class: %{customdata[3]}<br><br>"
        "Center (50%): %{customdata[4]:.4f}<br>"
        "Center CI: [%{customdata[5]:.4f}, %{customdata[6]:.4f}]<br>"
        "Gradient (50%): %{customdata[7]:.4f}<br>"
        "Gradient CI: [%{customdata[8]:.4f}, %{customdata[9]:.4f}]<br>"
        "u (50%): %{customdata[10]:.4f}<extra></extra>"
    )


def make_cline_overlay_plot(df):
    h = np.linspace(0.001, 0.999, 400)
    fig = go.Figure()
    n = df.shape[0]

    # Colored traces
    seen_classes = set()
    for _, row in df.iterrows():
        phi = compute_phi(h, row["center_clipped"], row["gradient_50"])
        customdata = make_hover_customdata_line(row, len(h))

        show_legend = row["plot_class"] not in seen_classes
        if show_legend:
            seen_classes.add(row["plot_class"])

        fig.add_trace(
            go.Scatter(
                x=h,
                y=phi,
                mode="lines",
                line=dict(color=row["color"], width=1.6),
                name=row["plot_class"],
                legendgroup=row["plot_class"],
                showlegend=show_legend,
                visible=True,
                customdata=customdata,
                hovertemplate=cline_hovertemplate(),
            )
        )

    # Black traces
    for _, row in df.iterrows():
        phi = compute_phi(h, row["center_clipped"], row["gradient_50"])
        customdata = make_hover_customdata_line(row, len(h))
        fig.add_trace(
            go.Scatter(
                x=h,
                y=phi,
                mode="lines",
                line=dict(color="rgba(0,0,0,0.45)", width=1.2),
                name=row["locus"],
                showlegend=False,
                visible=False,
                customdata=customdata,
                hovertemplate=cline_hovertemplate(),
            )
        )

    fig.add_trace(
        go.Scatter(
            x=[0, 1],
            y=[0, 1],
            mode="lines",
            line=dict(color="black", width=1, dash="dash"),
            name="y=x",
            hovertemplate="Hybrid index: %{x:.3f}<br>Ancestry probability: %{y:.3f}<extra></extra>",
            visible=True,
        )
    )

    fig.update_layout(
        title="Genomic Clines",
        template="plotly_white",
        xaxis_title="Hybrid index",
        yaxis_title="Ancestry probability",
        xaxis=dict(range=[0, 1]),
        yaxis=dict(range=[0, 1]),
        hovermode="closest",
        showlegend=True,
        legend=dict(
            title="Outlier class",
            x=1.02,
            y=1.0,
            xanchor="left",
            yanchor="top",
            bgcolor="rgba(255,255,255,0.85)",
        ),
        margin=dict(t=120, r=220),
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                x=0.5,
                y=1.20,
                xanchor="center",
                yanchor="top",
                buttons=[
                    dict(
                        label="Color by outlier type",
                        method="update",
                        args=[
                            {"visible": ([True] * n) + ([False] * n) + [True]},
                            {"title": "Genomic Clines (colored by outlier type)"},
                        ],
                    ),
                    dict(
                        label="All black",
                        method="update",
                        args=[
                            {"visible": ([False] * n) + ([True] * n) + [True]},
                            {"title": "Genomic Clines (all loci)"},
                        ],
                    ),
                ],
            )
        ],
    )
    return fig


def monotone_chain(points):
    """
    Returns hull points in order as an array of coordinates.
    """
    pts = np.array(sorted(set(map(tuple, points))))
    if len(pts) < 3:
        return None

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(tuple(p))

    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(tuple(p))

    hull = lower[:-1] + upper[:-1]
    if len(hull) < 3:
        return None
    return np.array(hull)


def add_convex_hull_trace(fig, xvals, yvals, idx, color, name):
    if len(idx) < 3:
        return
    pts = np.column_stack([xvals[idx], yvals[idx]])
    hull = monotone_chain(pts)
    if hull is None or hull.shape[0] < 3:
        return

    hull_x = np.append(hull[:, 0], hull[0, 0])
    hull_y = np.append(hull[:, 1], hull[0, 1])

    fig.add_trace(
        go.Scatter(
            x=hull_x,
            y=hull_y,
            mode="lines",
            fill="toself",
            line=dict(color=color, width=1),
            fillcolor=rgba_with_alpha(color, 0.18),
            name=f"{name} hull",
            hoverinfo="skip",
            showlegend=False,
        )
    )


def make_param_scatter_plot(df, x_mode="center", y_mode="v", add_hulls=True):
    if x_mode == "center":
        xx = df["center_50"].to_numpy()
        xlab = "Center (c)"
    else:
        xx = df["u_50"].to_numpy()
        xlab = "u = logit(c) * v"

    if y_mode == "v":
        yy = df["gradient_50"].to_numpy()
        ylab = "Gradient (v)"
    else:
        yy = df["log10_gradient_50"].to_numpy()
        ylab = "log10(Gradient)"

    fig = go.Figure()

    class_order = [
        "Non-outlier",
        "Steeper",
        "Shallower",
        "P1-bias",
        "P0-bias",
        "P1-bias;Steeper",
        "P1-bias;Shallower",
        "P0-bias;Steeper",
        "P0-bias;Shallower",
    ]

    symbol_map = {
        "Non-outlier": "circle",
        "Steeper": "triangle-up",
        "Shallower": "triangle-down",
        "P1-bias": "square",
        "P0-bias": "diamond",
        "P1-bias;Steeper": "square",
        "P1-bias;Shallower": "square",
        "P0-bias;Steeper": "diamond",
        "P0-bias;Shallower": "diamond",
    }

    for cls in class_order:
        sub = df[df["plot_class"] == cls]
        if sub.empty:
            continue

        customdata = np.column_stack([
            sub["locus"],
            sub["row_id"],
            sub["outlier_type"],
            sub["plot_class"],
            sub["center_50"],
            sub["center_low"],
            sub["center_high"],
            sub["gradient_50"],
            sub["gradient_low"],
            sub["gradient_high"],
            sub["u_50"],
        ])

        fig.add_trace(
            go.Scatter(
                x=sub["center_50"] if x_mode == "center" else sub["u_50"],
                y=sub["gradient_50"] if y_mode == "v" else sub["log10_gradient_50"],
                mode="markers",
                name=cls,
                marker=dict(
                    size=8,
                    color=sub["color"].iloc[0],
                    symbol=symbol_map.get(cls, "circle"),
                    line=dict(width=0.5, color="rgba(0,0,0,0.5)")
                ),
                customdata=customdata,
                hovertemplate=param_hovertemplate(),
            )
        )

    if add_hulls:
        for cls in class_order:
            if cls == "Non-outlier":
                continue
            idx = np.where(df["plot_class"].to_numpy() == cls)[0]
            if len(idx) >= 3:
                add_convex_hull_trace(fig, xx, yy, idx, color_for_class(cls), cls)

    fig.update_layout(
        title="Genomic Cline Parameters",
        template="plotly_white",
        xaxis_title=xlab,
        yaxis_title=ylab,
        hovermode="closest",
        legend_title="Outlier class",
    )
    return fig


def make_single_locus_plot(df):
    h = np.linspace(0.001, 0.999, 400)
    fig = go.Figure()
    n = df.shape[0]

    for i, (_, row) in enumerate(df.iterrows()):
        phi_med = compute_phi(h, row["center_clipped"], row["gradient_50"])

        center_low = min(max(row["center_low"], 1e-12), 1 - 1e-12)
        center_high = min(max(row["center_high"], 1e-12), 1 - 1e-12)
        grad_low = max(row["gradient_low"], 1e-12)
        grad_high = max(row["gradient_high"], 1e-12)

        # Approximate uncertainty envelope from corner combinations
        phi_candidates = np.vstack([
            compute_phi(h, center_low, grad_low),
            compute_phi(h, center_low, grad_high),
            compute_phi(h, center_high, grad_low),
            compute_phi(h, center_high, grad_high),
        ])
        band_lower = phi_candidates.min(axis=0)
        band_upper = phi_candidates.max(axis=0)

        customdata = make_hover_customdata_line(row, len(h))

        fig.add_trace(
            go.Scatter(
                x=h,
                y=band_lower,
                mode="lines",
                line=dict(color="rgba(0,0,0,0)"),
                hoverinfo="skip",
                showlegend=False,
                visible=(i == 0),
            )
        )

        fig.add_trace(
            go.Scatter(
                x=h,
                y=band_upper,
                mode="lines",
                line=dict(color="rgba(0,0,0,0)"),
                fill="tonexty",
                fillcolor="rgba(150,150,150,0.25)",
                hoverinfo="skip",
                showlegend=False,
                visible=(i == 0),
            )
        )

        fig.add_trace(
            go.Scatter(
                x=h,
                y=phi_med,
                mode="lines",
                line=dict(color="black", width=2),
                name=row["locus"],
                visible=(i == 0),
                customdata=customdata,
                hovertemplate=cline_hovertemplate(),
            )
        )

    fig.add_trace(
        go.Scatter(
            x=[0, 1],
            y=[0, 1],
            mode="lines",
            line=dict(color="gray", width=1, dash="dash"),
            name="y=x",
            visible=True,
            hovertemplate="Hybrid index: %{x:.3f}<br>Ancestry probability: %{y:.3f}<extra></extra>",
        )
    )

    buttons = []
    for i, row in df.iterrows():
        visible = [False] * (3 * n) + [True]
        visible[3 * i] = True
        visible[3 * i + 1] = True
        visible[3 * i + 2] = True

        buttons.append(
            dict(
                label=row["locus"],
                method="update",
                args=[
                    {"visible": visible},
                    {"title": f"Genomic Cline: {row['locus']}"},
                ],
            )
        )

    fig.update_layout(
        title=f"Genomic Cline: {df.iloc[0]['locus']}",
        template="plotly_white",
        xaxis_title="Hybrid index",
        yaxis_title="Ancestry probability",
        xaxis=dict(range=[0, 1]),
        yaxis=dict(range=[0, 1]),
        hovermode="closest",
        updatemenus=[
            dict(
                buttons=buttons,
                direction="down",
                showactive=True,
                x=1.02,
                xanchor="left",
                y=1.0,
                yanchor="top",
            )
        ],
    )
    return fig


def main():
    parser = argparse.ArgumentParser(
        description="Generate three genomic cline Plotly HTML outputs with embedded MultiQC metadata."
    )
    parser.add_argument("--center", required=True, help="Center posterior summary TSV")
    parser.add_argument("--gradient", required=True, help="Gradient posterior summary TSV")
    parser.add_argument("--loci", required=True, help="Locus order file")

    parser.add_argument("--template-overlay", required=True, help="Template for cline overlay plot")
    parser.add_argument("--template-scatter", required=True, help="Template for parameter scatter plot")
    parser.add_argument("--template-single", required=True, help="Template for single-locus cline plot")

    parser.add_argument("--out-overlay", required=True, help="Output HTML for cline overlay plot")
    parser.add_argument("--out-scatter", required=True, help="Output HTML for parameter scatter plot")
    parser.add_argument("--out-single", required=True, help="Output HTML for single-locus cline plot")

    parser.add_argument("--out-table", required=True, help="Output TSV summary table for all loci")

    args = parser.parse_args()

    df = prepare_dataframe(args.center, args.gradient, args.loci)

    write_summary_table(df, args.out_table)

    overlay_meta = parse_template(args.template_overlay)
    scatter_meta = parse_template(args.template_scatter)
    single_meta = parse_template(args.template_single)

    overlay_fig = make_cline_overlay_plot(df)
    scatter_fig = make_param_scatter_plot(df, x_mode="center", y_mode="v", add_hulls=True)
    single_fig = make_single_locus_plot(df)

    write_html(overlay_fig, args.out_overlay, overlay_meta)
    write_html(scatter_fig, args.out_scatter, scatter_meta)
    write_html(single_fig, args.out_single, single_meta)


if __name__ == "__main__":
    main()
