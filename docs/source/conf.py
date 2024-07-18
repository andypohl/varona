# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
from enum import Enum

from sphinx.application import Sphinx

sys.path.insert(0, os.path.abspath("../../src"))

project = "Varona"
copyright = "2024, Andy Pohl"
author = "Andy Pohl"
release = "0.0.2"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

templates_path = ["_templates"]
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
    "sphinx.ext.coverage",
    "sphinx.ext.viewcode",
    "sphinx_rtd_theme",
    "sphinx_multiversion",
]
exclude_patterns = []
add_module_names = False
html_theme = "sphinx_rtd_theme"
# html_logo = "_static/logo.svg"
html_theme_options = {
    "titles_only": False,
    "logo_only": True,
    "collapse_navigation": False,
    # "navigation_depth": 2,
}
html_static_path = ["_static"]
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "pysam": ("https://pysam.readthedocs.io/en/latest/", None),
}

# Add custom CSS files
html_css_files = [
    "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0-beta3/css/all.min.css",
]
smv_tag_whitelist = r"^v\d+\.\d+\.\d+.*$"
smv_branch_whitelist = r"^((sphinx|docs).*|main)$"
smv_released_pattern = r"^tags/v.*$"


def setup(app: Sphinx):
    def autodoc_skip_member(app, what, name, obj, skip, options):
        # Check if the object is an enum member by checking its type's base classes
        if isinstance(obj, Enum):
            return True  # Skip all enum members
        return skip  # Otherwise, keep the original decision

    app.connect("autodoc-skip-member", autodoc_skip_member)
    # Add custom JavaScript file
    app.add_js_file("custom.js")
    app.add_css_file("fa.css")
