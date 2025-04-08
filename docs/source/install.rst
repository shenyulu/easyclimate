.. _install:

Installation Guide
====================================

Welcome to the **easyclimate** installation guide! 🚀 We're excited to help you get started with our powerful climate analysis tool.
Follow these simple steps to install **easyclimate** on your system.

The easyclimate package is currently built and tested for specific platforms due to compatibility and dependency constraints.
Below are the supported platforms and notes for users on other systems.

- Windows x86-64/AMD64
- Linux x86-64/AMD64

These platforms are fully tested, and pre-built wheels (``.whl``) are available on PyPI for easy installation via following methods:

.. tab-set::

    .. tab-item:: PyPI

        Using the `PyPI <https://pypi.org/project/pip/>`__ package manager:

        .. code:: bash

            python -m pip install easyclimate

        If you don't have ``pip`` installed, this `Python installation guide <https://docs.python-guide.org/starting/installation/>`__ can guide you through the process.

    .. tab-item:: conda/mamba

        🛠️ Support is coming soon! Stay tuned for updates—we're working on it!

    .. tab-item:: Development version

        You can use ``PyPI`` to install the latest **unreleased** version from
        GitHub (⚠️ **NOT recommended** in most situations):

        .. code:: bash

            python -m pip install --upgrade git+https://github.com/shenyulu/easyclimate@dev

        .. note::

            The commands above should be executed in a terminal. On Windows, use the
            ``cmd.exe`` or the "Anaconda Prompt" app if you're using Anaconda.

Python Version Requirement
------------------------------------

**easyclimate** requires **Python 3.10 or higher**. To check your Python version, run:

.. code:: bash

    python --version

Make sure you're up to date! 🐍


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

        .. tab-item:: Windows

            .. code:: powershell

                .\build_docs_windows.ps1

            .. hint::

                On Windows, we've included ``optipng.exe`` for you! 😉 You might **NOT** need to install `optipng <https://optipng.sourceforge.net/>`__ for image optimization.

        .. tab-item:: Linux

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

About easyclimate-backend
------------------------------------
`Easyclimate-backend <https://easyclimate-backend.readthedocs.io/>`__ is the *core* powerhouse behind the easyclimate front-end package,
providing a suite of high-performance,
low-level functions for climate data analysis. Implemented in languages like ``Fortran`` and ``C``,
these functions ensure that your climate data processing is both efficient and accurate.

Because of this, you may also need to install a pre-compiled package or compile it yourself on
`Windows <https://easyclimate-backend.readthedocs.io/en/latest/src/building_windows.html>`__,
`Linux <https://easyclimate-backend.readthedocs.io/en/latest/src/building_linux.html>`__, or
`manylinux package <https://easyclimate-backend.readthedocs.io/en/latest/src/building_manylinux.html>`__.
