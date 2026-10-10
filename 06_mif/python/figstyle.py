#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# mIF pipeline - shared plotting style
#
# Purpose   : The panel theme used by the mIF figures, the box/bar panel helper,
#             the compact-y-axis rule and the significance-bracket helpers.
#
# Inputs    : None (style definitions only).
#
# Outputs   : None.
#
# Usage     : import figstyle as fs
#             fs.prepare()   # apply the global rcParams once per script
#             fs.panel(ax, groups, present, style="box", title=..., ylabel=...)
# -----------------------------------------------------------------------------

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- Colours and groups ----------------
GROUP_ORDER = ["Stable", "Growth", "Solid"]
GROUP_LABELS = {"Stable": "Stable", "Growth": "GGO-growth", "Solid": "SOLID-growth"}
GROUP_COLORS = {"Stable": "#8DA0CB", "Growth": "#66C2A5", "Solid": "#FC8D62"}

# ---------------- Figure sizes (near-square; not tall and narrow) ----------------
FIG_PANEL = (2.35, 2.55)
FIG_PANEL_WIDE = (3.15, 2.55)
FIG_ROW3 = (7.3, 2.55)
FIG_ROW4 = (9.6, 2.55)
FIG_2X2 = (4.9, 5.3)
FIG_2X3 = (7.3, 5.3)
FIG_2X4 = (9.6, 5.3)
FIG_3X2 = (4.9, 7.2)

FS = {
    "font.family": "Arial",
    "font.size": 7.5,
    "axes.titlesize": 8.0,
    "axes.labelsize": 7.8,
    "xtick.labelsize": 7.2,
    "ytick.labelsize": 7.2,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.4,
    "ytick.major.size": 2.4,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "axes.titleweight": "normal",
    "axes.labelweight": "normal",
    "axes.titlepad": 3.5,
    "axes.labelpad": 2.0,
    "legend.frameon": False,
    "legend.handlelength": 1.2,
    "legend.borderaxespad": 0.2,
    "lines.linewidth": 0.8,
    "patch.linewidth": 0.5,
    "figure.dpi": 400,
    "savefig.dpi": 400,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
}


def prepare():
    plt.rcParams.update(FS)


def wrap_title(title, fig_w_in, fs_pt=7.6):
    """Wrap a long title to the figure width — otherwise the title stretches the
    bbox_inches='tight' output and the axes end up smaller than the canvas."""
    import textwrap
    per = max(20, int(fig_w_in / (0.52 * fs_pt / 72.0)))
    return "\n".join(textwrap.wrap(title, per)) or title


def classic(ax, grid=False, tick_fs=None):
    """theme_classic equivalent: keep only the left/bottom spines, thin lines,
    outward ticks; no grid by default.

    Font sizes are set explicitly here rather than inherited from rcParams, so a
    script-level font.size does not inflate the compact panel fonts.
    """
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.6)
    fs = FS["xtick.labelsize"] if tick_fs is None else tick_fs
    ax.tick_params(width=0.6, length=2.4, direction="out", pad=1.5, labelsize=fs)
    if grid:
        ax.yaxis.grid(True, linestyle=":", linewidth=0.4, alpha=0.45, color="gray")
        ax.set_axisbelow(True)
    else:
        ax.grid(False)
    return ax


# ---------------- Compact y axis ----------------
def _whisker_top(v):
    """Upper 1.5 x IQR whisker (clipped to the largest observed value)."""
    q1, q3 = np.percentile(v, [25, 75])
    return float(min(np.max(v), q3 + 1.5 * (q3 - q1)))


def tight_ylim(ax, groups, n_brackets=3, ymin=0.0, pad=0.05, step=0.115,
               hard_cap_ratio=1.75, show_clip_marker=False):
    """Compress the y-axis top to "range actually occupied by the data + room for
    the brackets".

    Rule: the top is scaled from each group's 1.5 x IQR whisker, not stretched by
    a single outlier or by the brackets. Points above the top are kept as filled
    triangles at the axis top (nothing is discarded); panel() places them at each
    group's real x. Returns dict(top, base, step_y, n_clipped, over) where `over`
    holds the out-of-range values per group.
    """
    groups = [np.asarray(g, dtype=float) for g in groups if len(g) > 0]
    if not groups:
        return None
    core_top = max(_whisker_top(g) for g in groups)
    data_max = max(float(np.max(g)) for g in groups)
    if core_top <= 0:
        core_top = data_max or 1.0
    ladder_h = core_top * (step * n_brackets + pad + 0.04)
    if data_max <= core_top * hard_cap_ratio:
        visible_top, clipped = data_max, 0
        base = data_max * 1.12
    else:
        visible_top = core_top
        base = core_top * 1.12
        clipped = int(sum(int((g > base + ladder_h).sum()) for g in groups))
    top = base + ladder_h
    ax.set_ylim(ymin, top)
    over = [g[g > top] for g in groups]          # out-of-range values per group (kept as-is)
    if clipped and show_clip_marker:
        # fallback marker: pin it to the top-right of the axes (axis-fraction
        # coordinates, so y must be 1.0, not the data-coordinate `top`)
        clip_marks(ax, [0.972], 1.0, transform=ax.transAxes)
    return dict(top=top, base=base, step_y=core_top * step, n_clipped=clipped, over=over)


def clip_marks(ax, xs, top, ms=3.2, color="black", transform=None):
    """Out-of-range points: draw a filled triangle at the matching x on the axis
    top (data are not lost).

    xs    — x of each out-of-range point (data coordinates; axis fractions when
            transform=ax.transAxes)
    top   — axis-top position (same coordinate system as xs); the triangle sits
            on this line
    """
    tr = transform or ax.transData
    for k, x in enumerate(xs):
        ax.plot([x], [top], marker="^", ms=ms, color=color, mec=color, mew=0,
                clip_on=False, zorder=8, transform=tr)


def clipped_marks(ax, vals, present, top, width=0.55, ms=3.3):
    """Draw one filled triangle per out-of-range point of each group, pinned to
    the axis top directly above that group's box.

    Matches the convention in 07_distance_shared_axis.py: points are neither
    dropped nor merged into a single symbol — the triangles are spread across the
    group's box width so the number of clipped points can be counted.
    """
    span = width * 0.40
    for i, v in enumerate(vals):
        v = np.asarray(v, dtype=float)
        k = int((v > top).sum())
        if k == 0:
            continue
        xs = (i + 1) + (np.linspace(-span, span, k) if k > 1 else np.array([0.0]))
        ax.plot(xs, np.full(k, top), marker="^", ms=ms, color="black", ls="none",
                mec="black", mew=0, clip_on=False, zorder=8)


# ---------------- Significance brackets ----------------
def _draw_bracket(ax, x1, x2, y, txt, transform, tip_frac, fs, lw):
    span = ax.get_ylim()[1] - ax.get_ylim()[0]
    if transform is None:                      # data coordinates
        t = span * tip_frac
        ys = [y - t, y, y, y - t]
    else:                                      # axis fractions (y in 0-1)
        t = tip_frac
        ys = [y - t, y, y, y - t]
    ax.plot([x1, x1, x2, x2], ys, lw=lw, color="black", solid_joinstyle="miter",
            clip_on=False, zorder=5, transform=(transform or ax.transData))
    off = span * 0.008 if transform is None else 0.012
    ax.text((x1 + x2) / 2.0, y + off, txt, ha="center", va="bottom", fontsize=fs,
            color="black", zorder=5, transform=(transform or ax.transData))


def _draw_vbracket(ax, x, y1, y2, txt, fs=5.4, lw=0.5, cap=0.045, txt_dx=0.06,
                   color="black", top_color=None, bot_color=None, gap=2.0):
    """Vertical significance marker (for line plots): a vertical line from y1 to
    y2 at x with short end caps, and a rotated P value at the top. The line
    length equals the actual difference between the two group values (the effect
    size is not exaggerated)."""
    ax.plot([x, x], [y1, y2], lw=lw, color=color, clip_on=False, zorder=6)
    for y, c in ((y1, bot_color), (y2, top_color)):
        ax.plot([x - cap, x], [y, y], lw=lw * 1.9 if c else lw,
                color=c or color, clip_on=False, zorder=6, solid_capstyle="butt")
    ax.annotate(txt, xy=(x + txt_dx, max(y1, y2)), xytext=(0, gap),
                textcoords="offset points", rotation=90, ha="center", va="bottom",
                fontsize=fs, color=color, zorder=7, annotation_clip=False)


def p_ladder(ax, present, pvals, base, step_y, fs=5.9, lw=0.5, tip_frac=0.030):
    """Linear axis: three pairwise brackets stacked from `base`, one step_y apart
    (data coordinates)."""
    order = [g for g in GROUP_ORDER if g in present]
    combos = [(order[i], order[j]) for i in range(len(order)) for j in range(i + 1, len(order))]
    for k, ((g1, g2), p) in enumerate(zip(combos, pvals)):
        x1, x2 = order.index(g1) + 1, order.index(g2) + 1
        txt = "P < 0.001" if p < 0.001 else f"P = {p:.3f}"
        _draw_bracket(ax, x1, x2, base + k * step_y, txt, None, tip_frac, fs, lw)


def p_ladder_axes(ax, present, pvals, ys=(0.80, 0.87, 0.94), fs=5.9, lw=0.5):
    """Log axis: brackets positioned in axis fractions (x data / y fraction) so
    they do not compete with the data for space."""
    tr = ax.get_xaxis_transform()
    order = [g for g in GROUP_ORDER if g in present]
    combos = [(order[i], order[j]) for i in range(len(order)) for j in range(i + 1, len(order))]
    for k, ((g1, g2), p) in enumerate(zip(combos, pvals)):
        x1, x2 = order.index(g1) + 1, order.index(g2) + 1
        txt = "P < 0.001" if p < 0.001 else f"P = {p:.3f}"
        _draw_bracket(ax, x1, x2, ys[k], txt, tr, 0.012, fs, lw)


# ---------------- Statistics (unchanged: KW + two-sided MWW + Bonferroni x3) ----------------
def group_stats(groups, present, bonferroni=True):
    from scipy.stats import kruskal, mannwhitneyu
    import itertools
    vals = [np.asarray(groups[g], dtype=float) for g in present]
    try:
        _, kw_p = kruskal(*vals)
    except ValueError:          # all values identical: zero variance, reported as no difference
        kw_p = 1.0
    pairs, praw = [], []
    for g1, g2 in itertools.combinations(present, 2):
        try:
            _, p = mannwhitneyu(groups[g1], groups[g2], alternative="two-sided")
        except ValueError:
            p = 1.0
        pairs.append((g1, g2))
        praw.append(p)
    padj = [min(p * len(praw), 1.0) for p in praw] if (bonferroni and praw) else praw
    return {"kw_p": kw_p, "pairs": list(zip(pairs, padj))}


# ---------------- Plot body ----------------
def box_jitter(ax, groups, present, width=0.55, jitter=0.075, alpha_box=0.85,
               pt_size=5.0, pt_alpha=0.55, seed=42, lw=0.5):
    rng = np.random.default_rng(seed)
    pos = [list(present).index(g) + 1 for g in present]
    bp = ax.boxplot([np.asarray(groups[g], dtype=float) for g in present],
                    positions=pos, widths=width, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", linewidth=0.9),
                    boxprops=dict(linewidth=lw, edgecolor="black"),
                    whiskerprops=dict(linewidth=lw, color="black"),
                    capprops=dict(linewidth=lw, color="black"))
    for i, g in enumerate(present):
        bp["boxes"][i].set_facecolor(GROUP_COLORS[g])
        bp["boxes"][i].set_alpha(alpha_box)
        y = np.asarray(groups[g], dtype=float)
        ax.scatter(pos[i] + rng.normal(0, jitter, size=len(y)), y, s=pt_size, c="black",
                   alpha=pt_alpha, linewidths=0, zorder=3, clip_on=True)
    _xticks(ax, present, pos)
    return bp


def bar_sem(ax, groups, present, width=0.55, show_points=True, seed=42,
            pt_size=4.5, pt_alpha=0.55):
    rng = np.random.default_rng(seed)
    pos = [list(present).index(g) + 1 for g in present]
    for i, g in enumerate(present):
        y = np.asarray(groups[g], dtype=float)
        m = float(np.mean(y))
        e = float(np.std(y, ddof=1) / np.sqrt(len(y))) if len(y) > 1 else 0.0
        ax.bar(pos[i], m, width=width, color=GROUP_COLORS[g], alpha=0.85,
               edgecolor="black", linewidth=0.5, zorder=2)
        ax.errorbar(pos[i], m, yerr=e, fmt="none", ecolor="black", elinewidth=0.5,
                    capsize=1.6, capthick=0.5, zorder=4)
        if show_points:
            ax.scatter(pos[i] + rng.normal(0, 0.075, size=len(y)), y, s=pt_size, c="black",
                       alpha=pt_alpha, linewidths=0, zorder=5)
    _xticks(ax, present, pos)


def _xticks(ax, present, pos):
    ax.set_xticks(pos)
    ax.set_xticklabels([GROUP_LABELS[g].replace("-", "-\n") if len(GROUP_LABELS[g]) > 10
                        else GROUP_LABELS[g] for g in present], linespacing=1.1,
                       fontsize=FS["xtick.labelsize"])
    ax.set_xlim(0.42, len(present) + 0.58)


def panel(ax, groups, present, style="box", title=None, ylabel=None, seed=42,
          log=False, kw_in_title=True, figsize=None, ymin=0.0, brackets=True, **kw):
    """One-stop panel: draw -> statistics -> brackets -> compact y axis -> axis
    style. Returns the stats dict.

    figsize  when given, resizes the whole figure to that size (single-panel plots).
    ymin     lower y limit (default 0; callers that passed an explicit ymin keep it).
    brackets False keeps only the KW p in the title (overview-style small panels,
             where brackets would be too cramped).
    """
    st = group_stats(groups, present)
    if style == "box":
        box_jitter(ax, groups, present, seed=seed)
    elif style == "bar":
        bar_sem(ax, groups, present, seed=seed)
    else:
        raise ValueError(style)

    vals = [groups[g] for g in present]
    pvals = [p for _, p in st["pairs"]] if brackets else []
    if log:
        geo = _log1p_axis(ax, vals, n_brackets=len(pvals))
        hi = max(float(np.max(np.asarray(v, dtype=float))) for v in vals)
        p_ladder_axes(ax, present, pvals, ys=_log_ladder_ys(ax, hi, len(pvals)))
    else:
        geo = tight_ylim(ax, vals, n_brackets=len(pvals), ymin=ymin)
        if geo["n_clipped"]:
            clipped_marks(ax, vals, present, geo["top"])
        p_ladder(ax, present, pvals, geo["base"], geo["step_y"])
    classic(ax)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=FS["axes.labelsize"])
    if title:
        if figsize is not None:
            title = wrap_title(title, figsize[0])
        if kw_in_title:
            kwp = "P < 0.001" if st["kw_p"] < 0.001 else f"P = {st['kw_p']:.3f}"
            ax.set_title(f"{title}\nKruskal-Wallis {kwp}", fontsize=7.6, linespacing=1.3)
        else:
            ax.set_title(title, fontsize=7.8)
    if figsize is not None:
        ax.figure.set_size_inches(*figsize)
    st["geo"] = geo
    return st


def _log1p_axis(ax, vals, data_top_frac=0.72, n_brackets=3):
    """log1p y axis for zero-inflated counts: ticks labelled with raw values and
    headroom reserved for brackets.

    On a log axis the data top must land at data_top_frac, which requires solving
    the axis top back through log1p; otherwise a small multiplicative margin is
    visually crushed to nothing.
    """
    import matplotlib.ticker as mticker
    allv = np.concatenate([np.asarray(v, dtype=float) for v in vals])
    hi = max(1.0, float(np.max(allv)))
    top = float(np.expm1(np.log1p(hi) / data_top_frac))
    ax.set_yscale("function",
                  functions=(lambda x: np.log1p(np.clip(x, 0, None)),
                             lambda x: np.expm1(x)))
    ax.set_ylim(0, top)
    ticks = [t for t in (0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)
             if t <= top * 1.02]
    if not ticks or ticks[-1] < 0.5 * hi:
        ticks.append(int(round(hi)))
    ax.set_yticks(sorted(set(ticks)))
    ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    return dict(top=top, mode="log", n_clipped=0)


def _log_ladder_ys(ax, hi, n=3):
    """Place bracket heights in axis fractions from where the data top actually
    sits on the axis."""
    top = ax.get_ylim()[1]
    frac_hi = np.log1p(hi) / np.log1p(top)
    fs = [0.12, 0.44, 0.76][:n] if n <= 3 else np.linspace(0.12, 0.76, n)
    return [frac_hi + (1 - frac_hi) * f for f in fs]
