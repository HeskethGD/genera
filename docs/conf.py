project = "Genera"
copyright = "2026, Genera contributors"
author = "Genera contributors"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "matplotlib.sphinxext.plot_directive",
]
nitpicky = True
exclude_patterns = ["_build"]
html_theme = "alabaster"

plot_include_source = True
plot_formats = [("png", 96), "pdf"]
plot_html_show_formats = False
plot_html_show_source_link = False
