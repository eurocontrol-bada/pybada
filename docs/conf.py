# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../src'))
sys.path.insert(0, os.path.abspath('../src/pyBADA'))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "pyBADA"
copyright = "2024, EUROCONTROL"
author = "EUROCONTROL"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',     # For automatic documentation from docstrings
    'sphinx.ext.napoleon',    # For NumPy and Google style docstrings
    'sphinx.ext.viewcode',    # Optional: Add links to highlighted source code
    'sphinx.ext.autosummary', # Optional: Generate summary tables
]

templates_path = ["_templates"]
exclude_patterns = []


# Set the default docstring format for all modules
docstring_format = 'restructuredtext'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
