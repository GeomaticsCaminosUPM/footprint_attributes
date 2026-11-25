"""Sphinx configuration for footprint_attributes."""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.abspath("../src"))

project = "footprint_attributes"
copyright = "2026, Ureña-Pliego M., Rodríguez-Saiz J., Núñez-Álvarez G., Marchamalo-Sacristán M., González-Rodrigo B."
author = "Ureña-Pliego M. et al."
release = "1.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "myst_parser",
    "nbsphinx",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = False
napoleon_use_rtype = False

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
autodoc_typehints = "description"
autodoc_mock_imports = ["scipy", "sklearn"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

html_theme = "furo"
html_title = "footprint_attributes"
html_static_path = ["_static"]
html_favicon = "_static/favicon.svg"
html_css_files = [
    "https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700;800&family=JetBrains+Mono:wght@400;600&display=swap",
    "custom.css",
]
html_theme_options = {
    "sidebar_hide_name": False,
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/GeomaticsCaminosUPM/footprint_attributes",
            "html": (
                '<svg stroke="currentColor" fill="currentColor" stroke-width="0" '
                'viewBox="0 0 16 16"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 '
                "0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01"
                "-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68"
                "-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28"
                "-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2"
                "-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36"
                ".09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 "
                "1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-"
                '.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-'
                '8-8-8z"/></svg>'
            ),
            "class": "",
        },
    ],
    "light_css_variables": {
        "font-stack": "Inter, 'Helvetica Neue', Arial, sans-serif",
        "font-stack--monospace": "'JetBrains Mono', ui-monospace, monospace",
        "color-brand-primary": "#111111",
        "color-brand-content": "#B24A2A",
        "color-background-primary": "#ffffff",
        "color-background-secondary": "#fafafa",
        "color-foreground-primary": "#111111",
        "color-background-border": "#e6e6e6",
        "color-sidebar-background": "#ffffff",
        "color-sidebar-background-border": "#e6e6e6",
        "color-link": "#B24A2A",
        "color-link--hover": "#111111",
        "color-link-underline": "transparent",
        "color-link-underline--hover": "#111111",
    },
}

# nbsphinx: notebooks live under examples/, referenced via ../examples/*.ipynb
nbsphinx_execute = "never"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "geopandas": ("https://geopandas.org/en/stable/", None),
    "shapely": ("https://shapely.readthedocs.io/en/stable/", None),
}
