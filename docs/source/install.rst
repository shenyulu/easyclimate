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

    .. tab-item:: uv

        1. Install `uv <https://docs.astral.sh/uv/>`__

        .. code:: bash

            pip install uv

        2. Install `numba <https://numba.pydata.org/>`__ and `numpy <https://numpy.org/>`__ solely.

        .. code:: bash

            pip install numba numpy

        3. Install ``easyclimate`` by `uv <https://docs.astral.sh/uv/>`__

        .. code:: bash

            uv pip install easyclimate

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

        .. tip::

            For users who have difficulty accessing GitHub, we set up a GIT official mirror to access

            .. code:: bash

                python -m pip install --upgrade git+https://gitee.com/shenyulu/easyclimate@dev

.. tip::

    For *Linux* users who use *older systems*, special attention is required! It is necessary to install a relatively recent version of the `gcc <https://gcc.gnu.org/projects/c-status.html>`__ and `g++ <https://gcc.gnu.org/projects/cxx-status.html>`__ compiler suites to meet the **C++17** requirements. The latest versions of these suites can also be installed as follows:

    .. tab-set::

        .. tab-item:: conda (Recommended)

            .. code:: bash

                conda install -c conda-forge gcc gxx

        .. tab-item:: Ubuntu/Debian

            1. Update package list

            .. code:: bash

                sudo apt update

            2. Install the ``build-essential`` package by typing the following command:

            .. code:: bash

                sudo apt install build-essential

            This command installs a bunch of new software packages, including gcc, g++, and make.

            3. You may also need to install the manual pages for development using GNU/Linux:

            .. code:: bash

                sudo apt-get install manpages-dev

            4. To verify whether the GCC compiler has been successfully installed, use the following command to print the GCC version:

            .. code:: bash

                gcc --version

        .. tab-item:: CentOS/RHEL/Fedora

            1. Start a login shell as the root user.

            .. code:: bash

                sudo -i

            2. Install the ``gcc``.

            .. code:: bash

                yum install gcc

            Press Y and Enter to confirm. Or use the command without confirmation:

            .. code:: bash

                yum -y install gcc

            3. Install the ``g++``.

            .. code:: bash

                yum install gcc-c++

            Press Y and Enter to confirm. Or use the command without confirmation:

            .. code:: bash

                yum -y install gcc-c++


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

About easyclimate-backend
------------------------------------
The `easyclimate-backend <https://easyclimate-backend.readthedocs.io/>`__ is the **core** powerhouse behind the easyclimate front-end package,
providing a suite of high-performance,
low-level functions for climate data analysis. Implemented in languages like ``Fortran`` and ``C``,
these functions ensure that your climate data processing is both efficient and accurate.

Because of this, you may also need to install a pre-compiled package or compile it yourself on
`Windows <https://easyclimate-backend.readthedocs.io/en/latest/src/building_windows.html>`__,
`Linux <https://easyclimate-backend.readthedocs.io/en/latest/src/building_linux.html>`__, or
`manylinux package <https://easyclimate-backend.readthedocs.io/en/latest/src/building_manylinux.html>`__.
