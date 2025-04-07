# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
#
import os
import sys
import time
import datetime
import warnings

# autodoc required
sys.path.insert(
    0, os.path.abspath("../../src")
)  # Source code dir relative to this file
import easyclimate as ecl

# get year
localtime = time.localtime(time.time())
str_year = str(localtime[0])

project = "easyclimate"
copyright = f"2022-{datetime.datetime.now().year}, shenyulu（深雨露）"
author = "shenyulu"
release = "v" + ecl.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "recommonmark",
    "sphinx_markdown_tables",
    "sphinx.ext.mathjax",
    "sphinxcontrib.jupyter",
    "sphinx_inline_tabs",  # Add inline tabbed content to your Sphinx documentation
    "sphinx_gallery.gen_gallery",  # Add inline tabbed content to your Sphinx documentation
    # Sphinx AutoAPI Method
    "autoapi.extension",
    # Links to documentation for other projects
    "sphinx.ext.intersphinx",
    # copy button
    "sphinx_copybutton",
    "sphinx.ext.githubpages",
    "sphinx_design",
]

templates_path = ["_templates"]
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "changes/*.rst"]


# -- Options for AutoAPI extension -------------------------------------------
autoapi_type = "python"
autoapi_dirs = ["../../src"]
autoapi_add_toctree_entry = False
autoapi_root = "technical/api"

# autodoc_typehints = 'description'
# autosummary_generate = True  # Turn on sphinx.ext.autosummary

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'furo'
html_theme = "sphinx_book_theme"
html_static_path = ["_static"]

# Theme
html_theme_options = {
    "sidebar_hide_name": True,
    "top_of_page_button": "edit",
    "last-updated": True,
    "repository_url": "https://github.com/shenyulu/easyclimate",
    "use_repository_button": True,
    "icon_links": [
        {
            "name": "StackOverflow",
            "url": "https://stackoverflow.com/questions/tagged/easyclimate",
            "icon": "fa-brands fa-stack-overflow",
            "type": "fontawesome",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/easyclimate/",
            "icon": "fa-brands fa-python",
            "type": "fontawesome",
        },
    ],
}

# Logo
html_logo = "_static/easyclimate-logo.svg"


# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
# today = ''
# Else, today_fmt is used as the format for a strftime call.
today_fmt = "%Y-%m-%d"
# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = today_fmt

# Add edit button
html_theme_options = {
    "repository_url": "https://github.com/shenyulu/easyclimate",
    "use_repository_button": True,
    "repository_branch": "main",
    "use_repository_button": True,
    "use_issues_button": True,
    "use_download_button": True,
    "use_sidenotes": True,
}

# settings for sphinx-gallery
sphinx_gallery_conf = {
    "examples_dirs": "./dynamic_docs",  # path to your example scripts
    "gallery_dirs": "./auto_gallery",  # path to where to save gallery generated output
    "image_scrapers": ("matplotlib",),
    "compress_images": (
        "images",
        "thumbnails",
    ),  # require install `optipng`, download from http://optipng.sourceforge.net/
    "line_numbers": False,  # Line number
    "promote_jupyter_magic": True,
    #  Controlling what output is captured
    "capture_repr": ("_repr_html_", "__repr__", "__str__"),
    "run_stale_examples": True,
    "min_reported_time": False,
    "download_all_examples": False,
    #  'show_memory': True,
    "show_signature": False,
}
# supress warnings in gallery output
# https://sphinx-gallery.github.io/stable/configuration.html
warnings.filterwarnings("ignore", category=UserWarning,
                        message='Matplotlib is currently using agg, which is a'
                                ' non-GUI backend, so cannot show the figure.')

# Linked document
# https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html#confval-intersphinx_mapping
# https://pydoctor.readthedocs.io/en/latest/sphinx-integration.html
# e.g., https://docs.python.org/3/objects.inv
intersphinx_mapping = {
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "xarray": ("https://docs.xarray.dev/en/stable", None),
    "geocat-viz": ("https://geocat-viz.readthedocs.io/en/latest", None),
    "geocat-comp": ("https://geocat-comp.readthedocs.io/en/latest/", None),
    "dask": ("https://docs.dask.org/en/latest", None),
    "python": ("https://docs.python.org/3", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
    "statsmodels": ("https://www.statsmodels.org/stable", None),
    "xeofs": ("https://xeofs.readthedocs.io/en/latest/", None),
    "metpy": ("https://unidata.github.io/MetPy/latest/", None),
}
