.. _contributing:

************************
Contributing to easyclimate
************************

.. figure:: _static/fig6.jpg
    :scale: 40%
    :align: center

    Photo by Thomas Despeyroux on Unsplash

Thank you for your interest in contributing to *easyclimate*.
This page summarizes the contribution workflow that matches the current repository layout,
tooling, and maintenance practices.

*easyclimate* welcomes many kinds of contributions, including:

* bug reports
* feature suggestions
* bug fixes
* new tests
* documentation improvements
* examples and tutorials
* code cleanup and refactoring

All contributors are expected to follow the
`Code of Conduct <https://github.com/shenyulu/easyclimate/blob/main/CODE_OF_CONDUCT.md>`_.


Overview
========

The main development resources for this project are:

* `Source repository <https://github.com/shenyulu/easyclimate>`_
* `Issue tracker <https://github.com/shenyulu/easyclimate/issues>`_
* `Documentation site <https://easyclimate.readthedocs.io/en/latest/>`_
* `Contribution guide in Markdown <https://github.com/shenyulu/easyclimate/blob/main/CONTRIBUTING.md>`_

If you are new to the project, a good way to start is:

1. Read the user documentation.
2. Browse open issues.
3. Pick a small documentation fix, test improvement, or isolated bug.
4. Open a pull request from your fork.


How to ask questions
====================

If you need clarification before contributing:

1. Check whether the answer already exists in the documentation.
2. Search existing `GitHub issues <https://github.com/shenyulu/easyclimate/issues>`_.
3. If needed, open a new issue and describe your question clearly.

When asking, include relevant details such as:

* your operating system
* Python version
* *easyclimate* version
* a minimal code example, if applicable


Reporting bugs
==============

Good bug reports make it much easier to reproduce and fix problems.

Before opening an issue:

* verify that you are using a recent version of *easyclimate*
* search for existing reports describing the same problem
* check whether the problem can be reproduced with a minimal example

When filing a bug report, include:

* what you expected to happen
* what actually happened
* a minimal reproducible example
* traceback or error message, if available
* system information and package versions

You can print version information with:

.. code:: python

    import easyclimate as ecl

    ecl.show_versions()

For security-sensitive issues, do not disclose details publicly in the issue tracker.
Instead, contact the maintainer privately using the contact information in the project metadata.


Suggesting enhancements
=======================

Enhancement proposals are also handled through
`GitHub issues <https://github.com/shenyulu/easyclimate/issues>`_.

An effective enhancement request should explain:

* the current limitation
* the proposed improvement
* why the change is useful for multiple users
* possible API or behavior changes
* example usage when relevant

If you already plan to implement the enhancement, mentioning that in the issue or in the pull
request description is helpful.


Development workflow
====================

The recommended workflow is based on GitHub forks.

1. Fork the repository on GitHub.
2. Clone your fork locally.
3. Add the main repository as ``upstream``.
4. Create a dedicated branch for your work.
5. Commit focused changes with clear messages.
6. Push the branch to your fork.
7. Open a pull request.

Example setup:

.. code:: bash

    git clone https://github.com/your-user-name/easyclimate.git
    cd easyclimate
    git remote add upstream https://github.com/shenyulu/easyclimate.git
    git fetch upstream --tags

Create a branch before editing:

.. code:: bash

    git checkout -b docs/update-contributing-guide


Repository layout
=================

Understanding the repository structure helps contributors make targeted changes.

Important directories and files include:

* ``src/`` — package source code
* ``test/`` — automated tests and test data
* ``docs/source/`` — Sphinx documentation sources
* ``docs/requirements.txt`` — documentation build dependencies
* ``pyproject.toml`` — package metadata and build system configuration
* ``test_requirements.txt`` — test dependencies
* ``.pre-commit-config.yaml`` — configured pre-commit hooks


Setting up a development environment
====================================

The project supports Python ``>=3.10`` and currently declares classifiers for Python 3.10 through 3.14.
For local development, Python 3.12 is a practical default because it is also used in the pre-commit configuration.

A conda-based workflow fits the existing documentation and scripts well.

Windows PowerShell
------------------

.. code:: powershell

    conda create -n easyclimate-dev python=3.12
    conda activate easyclimate-dev
    pip install uv numpy numba
    uv pip install -r .\test_requirements.txt -r .\docs\requirements.txt build setuptools wheel setuptools-scm
    python -m build --wheel --no-isolation
    pip install .\dist\easyclimate-*.whl

Linux Bash
----------

.. code:: bash

    conda create -n easyclimate-dev python=3.12
    conda activate easyclimate-dev
    pip install uv numpy numba
    uv pip install -r ./test_requirements.txt -r ./docs/requirements.txt build setuptools wheel setuptools-scm
    python -m build --wheel --no-isolation
    pip install ./dist/easyclimate-*.whl

Notes:

* Documentation-only changes do not always require a full wheel build, but using the same environment is convenient.
* The project includes helper scripts under ``docs/`` such as ``build_docs_windows.ps1`` and ``build_docs_linux.sh``.


Running tests
=============

The test suite is located in ``test/`` and uses ``pytest``.
For most contributions, running the most relevant subset of tests is sufficient before opening a pull request.

Examples:

.. code:: bash

    pytest
    pytest test/test_core_stat.py
    pytest test/test_plot_bar.py

If your change affects calculations, interpolation, plotting, or physics utilities, add or update tests whenever possible.


Documentation contributions
===========================

Documentation sources live in ``docs/source/`` and are written primarily in reStructuredText.

Typical documentation contributions include:

* fixing wording or typos
* improving API explanations
* adding examples
* aligning docs with current behavior
* updating contribution and installation instructions

To build documentation locally, install dependencies from ``docs/requirements.txt`` and use one of the helper scripts in ``docs/``.

Examples:

* ``docs/build_docs_windows.ps1``
* ``docs/build_docs_linux.sh``

When editing docs:

* keep section titles and cross-references stable where possible
* prefer concise examples that users can reproduce
* update related pages if behavior or file paths changed


Code style and pre-commit
=========================

This repository includes a ``pre-commit`` configuration in ``.pre-commit-config.yaml``.
Configured hooks currently include basic file checks and ``black`` formatting for Python files,
plus a local deprecated-function check.

Install and run hooks locally:

.. code:: bash

    pre-commit install
    pre-commit run --all-files

Current style expectations:

* keep changes focused and easy to review
* use descriptive names
* preserve compatibility with the supported Python versions
* format Python code consistently
* avoid unrelated refactoring in the same pull request


Pull requests
=============

When your changes are ready:

1. Rebase or merge the latest upstream changes if needed.
2. Push your branch to your fork.
3. Open a pull request against the main repository.
4. Clearly explain what changed and why.
5. Link the relevant issue when applicable.

Helpful pull request descriptions usually include:

* scope of the change
* motivation
* testing performed
* screenshots for documentation or plotting changes, if useful

Keep pull requests as focused as possible. Smaller PRs are easier to review and merge.


Review checklist
================

Before submitting a pull request, verify the following:

* the change addresses a specific problem or documentation gap
* the modified code or docs are limited to the intended scope
* relevant tests were added or updated when needed
* local checks or pre-commit hooks pass
* documentation was updated for user-visible changes
* commit history is readable


Additional references
=====================

* `Project README <https://github.com/shenyulu/easyclimate/blob/main/README.md>`_
* `Issue tracker <https://github.com/shenyulu/easyclimate/issues>`_
* `Code of Conduct <https://github.com/shenyulu/easyclimate/blob/main/CODE_OF_CONDUCT.md>`_
* `pyproject.toml <https://github.com/shenyulu/easyclimate/blob/main/pyproject.toml>`_
