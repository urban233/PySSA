

class PyMOL_LOCK:
    _state: bool

    def __init__(self):
        self._state: bool = False

    def lock(self):
        """Sets the state to True."""
        self._state = True

    def unlock(self):
        """Sets the state to False."""
        self._state = False

    def is_locked(self):
        """Checks if the lock is locked.

        Notes:
            The lock is locked if the state is True.
        """
        return self._state
