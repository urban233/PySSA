# PySSA Code Rules
Authors: Hannah Kullik & Martin Urban

## Contents of this document
* [Description](#Description)
* [Linting](#Linting)
* [Type annotation](#Type-annotation)
* [Imports](#Imports)
* [Exception handling](#Exception-handling)
* [Terminology](#Terminology)
* [Code formatting](#Code-formatting)

## Description
Python is the main programming language of the PySSA project. 
This document describes the rules which must be followed if the 
source code gets extended.

## Linting
You have to run `ruff` over your code to check any static errors.
The configuration to use is defined in the `pyproject.toml`.

## Type annotation
Python is a dynamically typed language, but in this project 
Python is used as a **statically typed** language.
The decision emphasizes robust and less error-prone code.
Therefore, you have to use Python's type annotation feature.

### Annotations of python builtins
Annotating variables using python builtins where it is possible.
```python
i: int = 0
```
### Annotations of pyssa builtins
Annotating variables using pyssa builtins where data structures of 
pyssa are used.
```python
protein_pairs_for_analysis: list['protein_pair.ProteinPair'] = []
```
### Annotations of library builtins
Annotating variables using library builtins where data types of 
libraries are used.
```python
import numpy as np

distances_of_amino_acid_pairs: np.ndarray = np.ndarray([])
```

### Annotations of return values
If a function/ method has a return value that will not be used, that 
function call needs to be wrapped inside the `void` function.
The `void` function is the only function which gets imported as function 
and not as module:
```python
from pyssa.util.void import void
from pyssa.util import main_window_util

void(main_window_util.setup_app_settings(self.app_settings))  # void indicates that there is a 
                                                              # return value but it is not used
```

## Naming conventions
* Package: snake_case
* Module: snake_case
* Class: PascalCase
* Method: snake_case
  * private: _ prefix (single underscore)
  ```python
  def _create_directory_structure(self) -> None:
  ```
* Function: snake_case
* Variable: snake_case
  * argument: a/an_var_name, if no specific variable is meant.
  ```python
  def export_protein_as_pdb_file(a_filepath: str) -> None:
  ```
  * argument: the_var_name, if a specific variable is meant.
  ```python
  def load_settings(the_app_settings: 'settings.Settings') -> None:
  ```
  * method/function scope: tmp_ prefix
  ```python
  ...
  tmp_destination_filepath: str = "/home/rhel_user/scratch/log.txt"
  ...
  ```
* Global variable: g_ prefix + snake_case

## Imports
Never use wildcard imports. Always import the module **not** the class itself.
```python
from pymol import cmd # Correct: Module is imported

from pymol import * # Wrong! Wildcard import
from os.path import exists # Wrong! Function/Class import
```
Use official abbreviations for common python libraries.
```python
import numpy as np
import pandas as pd
```

## Exception handling

### Argument checks
Always check for None:
```python
def copy_fasta_file(a_source_filepath, a_destination_filepath):
    if a_source_filepath is None:
        logger.error(f"The argument 'a_source_filepath' is illegal: {a_source_filepath}!")
        raise exceptions.IllegalArgumentError("An argument is illegal.")
    if a_destination_filepath is None:
        logger.error(f"The argument 'a_destination_filepath' is illegal: {a_destination_filepath}!")
        raise exceptions.IllegalArgumentError("An argument is illegal.")
```

Raise **IllegalArgumentError** if *unmodified* argument 
is **not** usable for the function/method:
```python
import os


def copy_fasta_file(a_source_filepath: pathlib.Path, a_destination_filepath: pathlib.Path):
  ...
  if not os.path.exists(a_source_filepath):  # argument is unmodified
    raise FileNotFoundError()
```

Raise custom exception if argument
is *modified* **and** is **not** usable for the function/method
```python
import os


def copy_fasta_file(a_source_filepath: pathlib.Path, a_destination_filepath: pathlib.Path):
  ...
  if not os.path.exists(a_source_filepath.parent):  # .parent is a modified version of the argument
    raise exceptions.DirectoryNotFoundError("")
```

### try-except blocks
Always wrap `cmd` commands of the PyMOL api into a try-except block.

```python
import pymol

try:
  cmd.scene(f"{tmp_protein_pair.protein_1.get_molecule_object()}"
            f"{tmp_protein_pair.protein_2.get_molecule_object()}",
            action="recall")
except pymol.CmdException:
    logger.error("...")
    raise ...
```


## Terminology
### Path, dir, file & filepath
* Always use `path` if a directory path is meant.
* Always use `dir` if a directory name is meant.
* Always use `filepath` if an absolute path to a file is meant.
* Always use `file` if a name of a file is meant.

### Difference between TODO and fixme
* Add a `# TODO` if there is a task which needs to be done.
* Add a `# fixme` if there is an important note which needs to be quickly found.

## Code formatting
Always wrap argument checks into an editor-fold (Ctrl+Alt+T) and 
insert a line break before **and** after the ending of the editor-fold.
Example:
```python
# <editor-fold desc="Checks">
if the_fasta_path is None:
    logger.error("The argument filename is illegal.")
    raise exception.IllegalArgumentError("")

# </editor-fold>

```
