Setup Dev Environment
=====================

Step 1
------
Clone the GitHub repository of TEA.

.. code-block:: powershell

    git clone https://github.com/urban233/tea.git

Step 2
------
Create an empty mamba environment.

.. code-block:: powershell

    mamba create -n tea-dev python=3.9

Step 3
------
Install the following packages: "PyQt5", "sphinx", "pydata-sphinx-theme", "ruff", "pyink", "pytest", "pytest-qt"

.. code-block:: powershell

    mamba activate tea-dev
    pip install PyQt5, sphinx, pydata-sphinx-theme, ruff, pyink, pytest, pytest-qt
