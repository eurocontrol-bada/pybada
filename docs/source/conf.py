# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath('../../src'))

# -- Project information -----------------------------------------------------

project = 'pybada'
author = 'EUROCONTROL'
copyright = "2024, EUROCONTROL"
release = '0.1.0'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinx_rtd_theme'
]

autosummary_generate = True
templates_path = ['_templates']
exclude_patterns = []

# add_module_names = False
modindex_common_prefix = ['pyBADA.']

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
# html_static_path = ['_static']