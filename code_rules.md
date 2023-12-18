## PySSA Code Rules

### Type annotation
#### Annotations of python builtins
```python
i: int = 0
```

#### Annotations of pyssa builtins
```python
analysis_list: list['protein_pair.ProteinPair'] = []
```

### Exception handling

#### Argument checks
Raise **IllegalArgumentError** if 
* *unmodified* argument is not usable for the function/method

Raise custom exception if
* argument value is *modified* **and** is not usable for the function/method
  * e.g.: a_filepath.parent -> DirectoryNotFoundError

##### Styling
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

#### try-except blocks


### TODO's
#### Difference between TODO and fixme
Add a # TODO if there is a task which needs to be done.

Add a # fixme if there is an important note which needs to be quickly found.
