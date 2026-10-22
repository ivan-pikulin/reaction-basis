"""Documentation for the installed package; no import-path overrides."""

from importlib.metadata import version as package_version

project = "Reaction Basis"
release = package_version("reaction-basis")
version = release
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.viewcode",
]
html_theme = "sphinx_rtd_theme"
html_theme_options = {"navigation_depth": 3, "collapse_navigation": False}
html_title = f"{project} {release}"
html_show_sphinx = False
html_show_copyright = False
exclude_patterns = ["_build"]
nitpicky = True
autodoc_member_order = "bysource"
autodoc_typehints = "none"
autodoc_class_signature = "separated"
autodoc_typehints_format = "short"
autodoc_default_options = {"show-inheritance": False}
