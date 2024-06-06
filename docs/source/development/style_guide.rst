Style Guide
===========

Description
-----------
Python is the main programming language of the TEA project.
This document describes the rules that must be followed if the
source code gets extended.

Pre-commit hook
---------------
Before committing any changes, a custom commit-hook will be run **automatically**.
You **must solve** any issues before committing any changes!

You can choose to run the pre-commit hook configuration by yourself beforehand by
running the command:

.. code-block:: PowerShell

    pre-commit run --all-files

All available pre-commit hooks are listed here: https://pre-commit.com/hooks.html.
The addition of a pre-commit hook **needs an approval** of an author.

To be able to run the pre-commit hooks on Windows, it could be
possible to update the SSL lib. To do this download this `OpenSSL Installer`_.

.. _OpenSSL Installer: https://w-hs.sciebo.de/s/8F0rm9ChFAYRZu9/download

Code formatting
---------------
The overall code formatting is done with the auto-formatter `pyink`_.
This will be done if the pre-commit hooks are ran.
Alternatively you could also run *pyink* during development.

.. _pyink: https://github.com/google/pyink

Linting
-------
You have to run `ruff`_ over your code to check any static errors.
The configuration to use is defined in the `pyproject.toml`.

.. _ruff: https://docs.astral.sh/ruff/

Type annotation
---------------
Python is a dynamically typed language, but in this project
Python is used as a **statically typed** language.
The decision emphasizes robust and less error-prone code.
Therefore, you have to use Python's type annotation feature.

Annotations of python built-ins
*******************************
Annotating variables using python builtins where it is possible.

.. code-block:: Python

    i: int = 0

Annotations of TEA built-ins
****************************
Annotating variables using TEA builtins where class of TEA are used

.. code-block:: Python

    tmp_actions: list['action.Action'] = []

Annotations of library built-ins
********************************
Annotating variables using library builtins where data types of
libraries are used.

.. code-block:: Python

    import numpy as np

    tmp_tasks: np.ndarray = np.ndarray([])

Naming conventions
------------------
* Package: snake_case
* Module: snake_case
* Class: PascalCase
* Method: snake_case
    * private: _ prefix (single underscore)

    .. code-block:: Python

        def _create_directory_structure(self) -> None:

* Function: snake_case
* Variable: snake_case
    * argument: a/an_var_name, if no specific variable is meant.

    .. code-block:: Python

        def run_action(an_action: "action.Action") -> None:

    * argument: the_var_name, if a specific variable is meant.

    .. code-block:: Python

        def save_task_result(the_task_result: 'task_result.TaskResult') -> None:

    * method/function scope: **tmp\_** prefix

    .. code-block:: Python

      ...
      tmp_destination_filepath: str = "/home/rhel_user/scratch/log.txt"
      ...

* Global variable: g\_ prefix + snake_case

Imports
-------
Never use wildcard imports. Always import the module **not** the class itself.

.. code-block:: Python

    from tea.util import tea_logging # Correct: Module is imported

    from tea.util import * # Wrong! Wildcard import
    from os.path import exists # Wrong! Function/Class import

Use official abbreviations for common python libraries.

.. code-block:: Python

    import numpy as np
    import pandas as pd

Terminology
-----------
Path, dir, file & filepath
**************************
* Always use ``path`` if a directory path is meant.
* Always use ``dir`` if a directory name is meant.
* Always use ``filepath`` if an absolute path to a file is meant.
* Always use ``file`` if a name of a file is meant.

Difference between TODO and fixme
*********************************
* Add a ``# TODO`` if there is a task which needs to be done.
* Add a ``# fixme`` if there is an important note which needs to be quickly found.

Editor folds
------------
Always wrap argument checks into an editor-fold (Ctrl+Alt+T) and
insert a line break before **and** after the ending of the editor-fold.
Example:

.. code-block:: Python

    # <editor-fold desc="Checks">
    if the_threadpool is None:
        logger.error("the_threadpool is None.")
        raise tea_exception.IllegalArgumentError("the_threadpool is None.")

    # </editor-fold>
