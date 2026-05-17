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
import os
from dotenv import load_dotenv

# autodoc required
sys.path.insert(
    0, os.path.abspath("../../src")
)  # Source code dir relative to this file
import easyclimate as ecl

# copyright
localtime = time.localtime(time.time())
str_year = str(localtime[0])

project = "easyclimate"
copyright = f"2022-{datetime.datetime.now().year}, Shenyulu（深雨露） and easyclimate developers"
author = "shenyulu and easyclimate developers"
release = "v" + ecl.__version__


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "recommonmark",
    "sphinx_markdown_tables",
    "sphinx.ext.mathjax",
    # Add inline tabbed content to your Sphinx documentation
    "sphinx_inline_tabs",
    "sphinx_gallery.gen_gallery",
    # Sphinx AutoAPI Method
    "autoapi.extension",
    # Links to documentation for other projects
    "sphinx.ext.intersphinx",
    # copy button
    "sphinx_copybutton",
    "sphinx.ext.githubpages",
    "sphinx_design",
    # Embedding icons from over 200,000 open-source vector icons (https://icon-sets.iconify.design/)
    "sphinx_iconify",
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

html_theme = 'shibuya'
html_static_path = ["_static"]

# Logo
html_logo = "_static/easyclimate_logo_mini.png"

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
# today = ''
# Else, today_fmt is used as the format for a strftime call.
today_fmt = "%Y-%m-%d"
# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = today_fmt

# Shibuya html options
html_context = {
    # Source files copy for `Copy page`
    "source_type": "github",
    "source_user": "shenyulu",
    "source_repo": "easyclimate",
}

# Shibuya theme options
html_theme_options = {
    # Main color
    "accent_color": "bronze",
    # Copy page
    "show_ai_links": True,
    "open_in_chatgpt": True,
    "open_in_claude": True,
    "open_in_perplexity": True,

    "github_url": "https://github.com/shenyulu/easyclimate",
    "repository_url": "https://github.com/shenyulu/easyclimate",
    "use_repository_button": True,
    "repository_branch": "main",
    "use_repository_button": True,
    "use_issues_button": True,
    "use_download_button": True,
    "use_sidenotes": True,
    "nav_links": [
        {
            "title": "Getting Started",
            "children": [
                {"title": "Overview", "url": "overview"},
                {"title": "Installation", "url": "install"},

            ],
        },
        {
            "title": "Gallery", "url": "auto_gallery/index"
        },
        {
            "title": "API Index", "url": "api_index/index"
        },
        {
            "title": "Reference",
            "children": [

                {"title": "Release Notes", "url": "changes"},
                {"title": "Cite", "url": "static_docs/cite"},
                {
                    "title": "Open Source Licenses",
                    "url": "https://easyclimate-backend.readthedocs.io/en/latest/src/softlist.html",
                    "external": True,
                },
            ],
        },
        {
            "title": "Community",
            "children": [
                {"title": "Contributing", "url": "contributing"},
                {
                    "title": "How to Contribute",
                    "url": "https://github.com/shenyulu/easyclimate/blob/main/CONTRIBUTING.md",
                    "external": True,
                },
                {
                    "title": "GitHub",
                    "url": "https://github.com/shenyulu/easyclimate",
                    "external": True,
                },
                {"title": "Sponsor", "url": "sponsor"},
            ],
        },
    ],
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

    # Rerunning stale examples
    "run_stale_examples": False,
    # "run_stale_examples": True,    # rebuild examples only

    "min_reported_time": False,
    "download_all_examples": False,
    #  'show_memory': True,
    "show_signature": False,
    'remove_config_comments': True,
    # Modules for which function level galleries are created.  In
    # this case sphinx_gallery and numpy in a tuple of strings.
    "doc_module": "easyclimate",
    # Insert links to documentation of objects in the examples
    "reference_url": {"easyclimate": None},
    'parallel': 4,
    # mini-galleries
    ## directory where function/class granular galleries are stored
    'backreferences_dir'  : 'gen_modules/backreferences',
    ## Modules for which function/class level galleries are created. In
    ## this case sphinx_gallery and numpy in a tuple of strings.
    'doc_module'          : ('sphinx_gallery', 'easyclimate'),
    ## Regexes to match objects to exclude from implicit backreferences.
    ## The default option is an empty set, i.e. exclude nothing.
    ## To exclude everything, use: '.*'
    'exclude_implicit_doc': {r'pyplot\.show'},
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
    "geopandas": ("https://geopandas.org/en/stable", None),
}
