## PySSA Code Rules

### Exception handling

#### Argument checks
Raise **IllegalArgumentError** if 
* *unmodified* argument is not usable for the function/method

Raise custom exception if
* argument value is *modified* **and** is not usable for the function/method
  * e.g.: a_filepath.parent -> DirectoryNotFoundError

#### try-except blocks


### TODO's
#### Difference between TODO and fixme
Add a # TODO if there is a task which needs to be done.

Add a # fixme if there is an important note which needs to be quickly found.
