# Configuration file for the Sphinx documentation builder.

import os
import sys

sys.path.insert(0, os.path.abspath("../../src"))

# -- Project information -----------------------------------------------------

project = "pybada"
author = "Your Name"
release = "0.1"  # Replace with your project's version

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx_gallery.gen_gallery",
]

sphinx_gallery_conf = {
    "examples_dirs": "../../examples",  # Path to your example scripts
    "gallery_dirs": "auto_examples",  # Path where to save generated output
    "filename_pattern": "^((?!skip_).)*$",  # Only include scripts not starting with 'skip_'
}

autosummary_generate = True
templates_path = ["_templates"]
exclude_patterns = []

# add_module_names = False
modindex_common_prefix = ["pyBADA."]

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
