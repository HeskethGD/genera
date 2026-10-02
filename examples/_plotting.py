"""Shared visual conventions for the optional example figures."""

TEAL = "#015758"
LINE_STYLES = (
    (0, (6, 3)),
    (0, (2, 2)),
    (0, (6, 2, 1, 2)),
    (0, (6, 2, 1, 2, 1, 2)),
)
MARKERS = ("o", "s", "^", "D")


def comparison_styles(index):
    """Pair a dashed numerical curve with distinctive analytic markers."""
    style = index % len(MARKERS)
    return (
        {"color": TEAL, "linestyle": LINE_STYLES[style]},
        {"color": TEAL, "linestyle": "none", "marker": MARKERS[style],
         "markersize": 4, "fillstyle": "none"},
    )
