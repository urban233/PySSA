# PySSA Code Rules
## Description
Python is the main programming language of the PySSA project.
This document describes the rules which must be followed if the
source code gets extended.
The rule of thumb is to follow the [Google style guide for Python](https://google.github.io/styleguide/pyguide.html), but anything 
mentioned in this document **overrides** the Google style guide.
These rules can be outdated or changed at any time.

## Linting
You have to run [ruff](https://docs.astral.sh/ruff/) over your code to check any static errors.
The configuration to use is defined in the `pyproject.toml`.
Do **NOT** commit any code which is not linted with `ruff` and all ruff errors 
are fixed.

## Type annotation
Python is a dynamically typed language, but this project
makes heavy use of type annotations.
The decision emphasizes robust and less error-prone code.
Therefore, you have to use Python's type annotation feature.

### Annotations of Python built-ins
Annotating variables using python built-ins where it is possible.
```python
i: int = 0
program_name: str = "PySSA"
```
### Annotations of PySSA built-ins
Annotating variables using PySSA built-ins where data structures of
PySSA are used.
```python
protein_pairs_for_analysis: list['protein_pair.ProteinPair'] = []
```
### Annotations of library built-ins
Annotating variables using library built-ins where data types of
libraries are used.
```python
import numpy as np

distances_of_amino_acid_pairs: np.ndarray = np.ndarray([])
```

### Annotations of return values
If a function/ method has a return value that will not be used, that
function call needs to be wrapped inside the `rvoid` function.
The `rvoid` function is the only function which gets imported as function
and not as module:

```python
from pyssa.util.void import rvoid
from pyssa.util import main_window_util

rvoid(main_window_util.setup_app_settings(self.app_settings))  
# void indicates that there is a return value, but it is not used
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
  * slots: __ prefix (double underscore)
  ```python
  def __slot_open_project(self) -> None:
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
        raise exception.IllegalArgumentError("An argument is illegal.")
    if a_destination_filepath is None:
        logger.error(f"The argument 'a_destination_filepath' is illegal: {a_destination_filepath}!")
        raise exception.IllegalArgumentError("An argument is illegal.")
```

Raise **IllegalArgumentError** if *unmodified* argument
is **not** usable for the function/method:
```python
import os


def copy_fasta_file(a_source_filepath: pathlib.Path, a_destination_filepath: pathlib.Path):
  ...
  if not os.path.exists(a_source_filepath):  # argument is unmodified
    raise exception.IllegalArgumentError("An argument is illegal.")
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
Always wrap `cmd` commands of the PyMOL API into a try-except block.

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

## Communication
### QMainWindow & QDialogs
The communication between any QMainWindow and QDialog is done with
signals and slots. This ensures that no unauthorized memory access violations occur.
#### How-to
1. Define a custom pyqtsignal in the QDialog class:
```python
...

class DialogAddModel(Qt.QtWidgets.QDialog):
    """Class for a dialog to add proteins to a project."""

    """
    A pyqtsignal that is used to hand-over the protein structure information.
    """
    return_value = pyqtSignal(tuple)  # this is a custom PyQt signal

    ...
```
2. Emit the signal where communication should occur.
```python
...

def add_model(self) -> None:
    """Emits a custom pyqtsignal and closes the dialog."""
    self.return_value.emit((self.ui.txt_add_protein.text(), True))
    self.close()

...
```
3. Connect the signal in the QMainWindow with the QDialog object and the slot function

```python
...


def add_existing_protein(self) -> None:
  """Opens a dialog to add an existing protein structure to the project."""
  self.tmp_dialog = dialog_add_model.AddProteinView()
  self.tmp_dialog.return_value.connect(self.post_add_existing_protein)  # here is the connection
  self.tmp_dialog.show()


...
```
4. Be sure that the slot function has the value of the signal as an function argument
```python
...

def post_add_existing_protein(self, return_value: tuple):  # in this case the value is a tuple
    ...
```

## Threading
Within PySSA the custom `tea` package will be used if multithreading is necessary.
Do **NOT** use the `Action` class directly only use the `Task` class!

### Usage

```python
...


def opens_project(self):
    """Initiates the task to open an existing project."""
    self._active_task = tasks.LegacyTask(self.__async_open_project, post_func=self.__await_open_project)
    self._active_task.start()


def __async_open_project(self) -> tuple:
    """Runs in the separate QThread and does CPU-bound work."""
    tmp_project_path = pathlib.Path(f"{self._workspace_path}/{self._view.ui.txt_open_selected_project.text()}")
    return ("result", project.Project.deserialize_project(tmp_project_path, self._application_settings))


def __await_post_project(self, a_result: tuple):
    """Runs after the QThread finished."""
    ...
```
The `Task` class gets an "async" function and optionally an "await" function.
The function that runs in the QThread must have the signature `__async`
(double underscore). The function that runs after the QThread finished must
have the signature `__await`.
This design decision is based on intuition because the `__async` function
runs **asynchronous** in the QThread and the `__await` function **waits**
for the QThread (`__async` function) to finish,

## Database
PySSA uses a SQLite database for every single project. 

### Basic interaction with DatabaseX classes
The interaction is managed through the `DatabaseManager` class.
The interaction with the manager from a controller is done through the 
`DatabaseThread` class. The `DatabaseThread` has a queue which 
accepts objects of the type `DatabaseOperation`. 

To run an INSERT statement from a controller, you have to create
a `DatabaseOperation` object with the `SQLQueryType` (in this case INSERT_...)
and put it into the queue of the `DatabaseThread`.

```python
def _delete_protein(self):
    """Deletes an existing protein from the project."""
    tmp_protein: "protein.Protein" = self._view.ui.proteins_tree_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE)
    
    # Below is the creation of the DatabaseOperation object
    tmp_database_operation = database_operation.DatabaseOperation(enums.SQLQueryType.DELETE_EXISTING_PROTEIN, 
                                                                  (0, tmp_protein.get_id()))
    # Here the DatabaseOperation object will be put into the queue of the DatabaseThread
    self._database_thread.put_database_operation_into_queue(tmp_database_operation)
    
    # -- The rest of the function
    self._interface_manager.get_current_project().delete_specific_protein(tmp_protein.get_molecule_object())
    self._interface_manager.refresh_protein_model()
    self._interface_manager.refresh_main_view()
```

### Adding new SQL queries to the DatabaseThread
Every SQL statement has to be implemented in the `DatabaseManager` class!
For proper functionality of the `DatabaseThread` class 
it is necessary to add the SQL statements from the 
database manager into a _wrapper_ function and map this function against 
an appropriate SQLQueryType enum.

An example for a wrapper function.
```python
@staticmethod
def __wrapper_delete_existing_protein(the_db_manager, the_buffered_data: tuple):
    # It is import to unpack the first element of the tuple with an _ !   
    _, tmp_protein_id = the_buffered_data  
    the_db_manager.delete_existing_protein(tmp_protein_id)
```

An example for the mapping process
```python
def _setup_operations_mapping(self):
    self._operations_mapping = {
        enums.SQLQueryType.INSERT_NEW_PROTEIN: self.__wrapper_insert_new_protein,
        enums.SQLQueryType.DELETE_EXISTING_PROTEIN: self.__wrapper_delete_existing_protein
    }
```
and the SQLQueryType enum class
```python
class SQLQueryType(enum.Enum):
    """An enum for all possible sql queries for the database thread."""
    INSERT_NEW_PROTEIN = 'insert_new_protein'
    DELETE_EXISTING_PROTEIN = 'delete_existing_protein'
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
The overall code formatting is done with the auto-formatter black.
This will be done if the pre-commit hooks are ran.

### Editor folds
Always wrap argument checks into an editor-fold (Ctrl+Alt+T).
Example:
```python
# <editor-fold desc="Checks">
if the_fasta_path is None:
    logger.error("The argument filename is illegal.")
    raise exception.IllegalArgumentError("")
# </editor-fold>

```
