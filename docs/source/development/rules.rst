Rules
=====

Versioning
----------
TEA uses a loose variant of semantic versioning (`SemVer`_) to govern deprecations,
API compatibility, and version numbering.

A TEA release number is made up of ``MAJOR.MINOR.PATCH``.

API breaking changes should only occur in **major** releases.
These changes will be documented, with clear guidance on what is changing,
why it’s changing, and how to migrate existing code to the new behavior.

Whenever possible, a deprecation path will be provided
rather than an outright breaking change.

TEA will introduce deprecations in **minor** releases. These deprecations will preserve the existing behavior while emitting a warning that provide guidance on:

* How to achieve similar behavior if an alternative is available
* The TEA version in which the deprecation will be enforced.

We will not introduce new deprecations in patch releases.

Deprecations will only be enforced in **major** releases.
For example, if a behavior is deprecated in TEA 1.1.0,
it will continue to work, with a warning, for all releases in the 1.x series.
The behavior will change and the deprecation removed in the next major release (2.0.0).

.. _SemVer: https://semver.org

Code style
----------
TEA's code style is nearly completely based on the `Google Python Style Guide`_
with a few exceptions.

.. _Google Python Style Guide: https://google.github.io/styleguide/pyguide.html

The additions and changes from the google style guide are described in this
document :doc:`style_guide`.


.. toctree::
    :hidden:

    style_guide