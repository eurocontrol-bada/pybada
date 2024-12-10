# Configuration file for the Sphinx documentation builder.

import os
import sys

sys.path.append(os.path.abspath("../../src"))

# -- Project information -----------------------------------------------------

project = "pyBADA"
author = "Henrich Glaser - Opitz"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx_gallery.gen_gallery",
    "sphinx_search.extension",
]

sphinx_gallery_conf = {
    "examples_dirs": "../../examples",  # Path to your example scripts
    "gallery_dirs": "auto_examples",  # Path where to save generated output
    "filename_pattern": "^((?!skip_).)*$",  # Only include scripts not starting with 'skip_'
    "reset_modules_order": "after",  # Reset modules after each example
    "run_stale_examples": True,  # Force re-execution of examples
    "remove_config_comments": True,  # Clean up the code blocks
    "capture_repr": ("_repr_html_",),
    "default_thumb_file": os.path.abspath(
        os.path.join(
            os.path.dirname(__file__), "_static", "default_thumbnail.png"
        )
    ),
}

autosummary_generate = True
templates_path = ["_templates"]
exclude_patterns = []

html_css_files = [
    'css/custom.css',
]

# add_module_names = False
modindex_common_prefix = ["pyBADA."]

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

copyright = "2024, EUROCONTROL"