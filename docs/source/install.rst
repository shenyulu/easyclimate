.. _install:

Installation Guide
====================================

Welcome to the **easyclimate** installation guide! 🚀 We're excited to help you get started with our powerful climate analysis tool.
Follow these simple steps to install **easyclimate** on your system.

The easyclimate package is currently built and tested for specific platforms due to compatibility and dependency constraints.
Below are the supported platforms and notes for users on other systems.

- **Windows x86-64/AMD64** (Windows 10+)
- **Linux x86-64/AMD64** (glibc 2.28 or later, including: Debian 10+, Ubuntu 18.10+, Fedora 29+, CentOS/RHEL 8+)

These platforms are fully tested, and pre-built wheels (``.whl``) are available on PyPI for easy installation via following methods:

.. tab-set::

    .. tab-item:: :iconify:`devicon:pypi` PyPI

        Using the `PyPI <https://pypi.org/project/pip/>`__ package manager:

        .. code:: bash

            python -m pip install easyclimate

        If you don't have ``pip`` installed, this `Python installation guide <https://docs.python-guide.org/starting/installation/>`__ can guide you through the process.

    .. tab-item:: :iconify:`material-icon-theme:uv` Astral uv

        1. Install `uv <https://docs.astral.sh/uv/>`__

        .. code:: bash

            pip install uv

        2. Install ``easyclimate`` by `uv <https://docs.astral.sh/uv/>`__

        .. code:: bash

            uv pip install easyclimate

    .. tab-item:: :iconify:`devicon:anaconda` conda/mamba

        🛠️ Support is coming soon! Stay tuned for updates—we're working on it!

    .. tab-item:: :iconify:`fluent-emoji-flat:hammer-and-wrench` Development version

        You can use ``PyPI`` to install the latest **unreleased** version from
        GitHub (⚠️ **NOT recommended** in most situations):

        .. code:: bash

            python -m pip install --upgrade git+https://github.com/shenyulu/easyclimate@dev

        .. note::

            The commands above should be executed in a terminal. On Windows, use the
            ``cmd.exe`` or the "Anaconda Prompt" app if you're using Anaconda.

        .. tip::

            For users who have difficulty accessing GitHub, we set up a GIT official mirror to access

            .. code:: bash

                python -m pip install --upgrade git+https://gitee.com/shenyulu/easyclimate@dev

.. warning::

    Unfortunately, *easyclimate currently does NOT officially support macOS*, including both Intel-based Macs and Apple Silicon (M-series) Macs.

Python Version Requirement
------------------------------------

**easyclimate** requires **Python 3.10 or higher**. To check your Python version, run:

.. code:: bash

    python --version

Make sure you're up to date! 🐍

.. tip::

    See more `Status of Python versions <https://devguide.python.org/versions/>`__.

.. _dependencies:

Dependencies
------------------------------------

**easyclimate** comes with all the necessary dependencies for a smooth experience. Here's what gets installed:

.. tab-set::

    .. tab-item:: Base requirements

        Essential packages for core functionality.

        .. literalinclude:: ../../release_requirements.txt

    .. tab-item:: Test requirements

        Packages needed for running tests.

        .. literalinclude:: ../../test_requirements.txt

    .. tab-item:: Docs build requirements

        Tools for building the documentation.

        .. literalinclude:: ../requirements.txt

Building the Documentation
------------------------------------

Want to build the documentation yourself? 📚 Follow these steps:

- Install the docs build requirements listed above.
    .. code:: bash

        pip install -r docs/requirements.txt
- Go to the ``docs`` directory.
- Run the build script:
    .. tab-set::

        .. tab-item:: Windows Powershell

            .. code:: powershell

                .\build_docs_windows.ps1

            .. hint::

                On Windows, we've included ``optipng.exe`` for you! 😉 You might **NOT** need to install `optipng <https://optipng.sourceforge.net/>`__ for image optimization.

        .. tab-item:: Linux Bash

            .. code:: bash

                ./build_docs_linux.sh

            .. hint::

                On Linux, you might need to install `optipng <https://optipng.sourceforge.net/>`__ for image optimization.

                .. code:: bash

                    sudo apt-get install optipng

.. tip::

    For more control, you need to clean the build directory, build the HTML documentation, and copy example notebooks.



We hope this guide makes installing **easyclimate** a breeze! If you have any questions or run into issues,
feel free to reach out. Happy climate analyzing! 🌍
