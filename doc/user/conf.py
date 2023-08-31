# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.append(os.path.abspath("C:/code/HydroChrono/doc/user/_ext"))

# -- Project information -----------------------------------------------------
project = 'HydroChrono'
copyright = ''
author = ''

# The full version, including alpha/beta/rc tags
release = '0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    #'hydrochrono_h5',
    #'matplotlib.sphinxext.mathmpl.math_directive',
    'matplotlib.sphinxext.mathmpl',
    'matplotlib.sphinxext.plot_directive',
    #'matplotlib.sphinxext.only_directives',
    'sphinxcontrib.bibtex',
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
    'sphinx.ext.intersphinx',
    'sphinx.ext.inheritance_diagram'
]

# Include TODO comments
todo_include_todos = True

# Bibtex file location
bibtex_bibfiles = ['C:/code/HydroChrono/doc/user/references.bib']

# matplotlib.sphinxext.plot_directive Configuration
# plot_basedir = 'C:/code/HydroChrono/doc/user'
# plot_working_directory = 'C:/code/HydroChrono_build/doc/user'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

ProjectName = f'HydroChrono'

# Breathe Configuration
#if 1 and 0:
#    extensions.append('breathe')
#    breathe_default_project = ProjectName
#    breathe_projects = { ProjectName : "..\\..\\..\\..\\HydroChrono\\xml" }


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'bizstyle'
#html_theme = 'haiku'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']
#html_logo = 'xxx.png'
#html_favicon = 'favicon.ico'

# html_theme_options = {
#    'full_logo': True
# }


rst_prolog = r"""
.. |UNIT_NONE| replace:: :math:`\rm  -`

.. |UNIT_TIME| replace:: :math:`\rm s`

.. |UNIT_LENGTH| replace:: :math:`\rm m`

.. |UNIT_TEMP| replace:: K

.. |PROP_TEMP| replace:: :math:`T`

"""
