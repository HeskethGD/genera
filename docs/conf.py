from pathlib import Path
import sys

from cycler import cycler

sys.path.insert(0, str(Path(__file__).parent.parent / "docs_ext"))
sys.path.insert(0, str(Path(__file__).parent.parent))

project = "Genera"
copyright = "2026, Genera contributors"
author = "Genera contributors"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "matplotlib.sphinxext.plot_directive",
    "genera_example",
]
nitpicky = True
exclude_patterns = ["_build"]
html_theme = "furo"
html_static_path = ["_static"]
html_logo = "_static/genera-logo.png"
html_title = "Genera documentation"
# Teal, mint, and gold sampled from the logo; keep backgrounds close to
# Furo's defaults so mathematical text and code remain easy to read.
html_theme_options = {
    "light_css_variables": {
        "color-brand-primary": "#015758",
        "color-brand-content": "#015758",
        "color-api-name": "#015758",
        "color-api-pre-name": "#015758",
        "color-link--visited": "#015758",
        "color-link--visited--hover": "#015758",
        "color-sidebar-background": "#f3f8f5",
        "color-background-hover": "#edf5ef",
        "color-background-border": "#d6e5dc",
    },
    "dark_css_variables": {
        "color-brand-primary": "#93c9a4",
        "color-brand-content": "#93c9a4",
        "color-api-name": "#93c9a4",
        "color-api-pre-name": "#93c9a4",
        "color-link--visited": "#93c9a4",
        "color-link--visited--hover": "#93c9a4",
        "color-sidebar-background": "#011d29",
        "color-background-hover": "#02333b",
        "color-background-border": "#235c4f",
    },
}
html_css_files = ["genera.css"]

plot_include_source = True
plot_formats = [("png", 96), "pdf"]
plot_html_show_formats = False
plot_html_show_source_link = False
# Distinguish curves by pattern and shape, keeping all strokes dark teal.
plot_rcparams = {
    "axes.prop_cycle": cycler(color=["#015758"] * 4)
        + cycler(linestyle=["-", "--", "-.", ":"]),
    "axes.facecolor": "white",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
}
plot_apply_rcparams = True


def track_example_dependencies(app, docname, source):
    """Re-read example pages when the code they exercise changes."""
    if not docname.startswith("examples/"):
        return
    root = Path(__file__).resolve().parent.parent
    for directory in ("src/genera", "examples", "docs_ext"):
        for path in (root / directory).rglob("*.py"):
            app.env.note_dependency(str(path))


def setup(app):
    app.connect("source-read", track_example_dependencies)
