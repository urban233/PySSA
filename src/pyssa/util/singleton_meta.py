#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/urban233/PySSA>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""Module for a singleton meta class type."""
import threading

class SingletonMeta(type):
  """
  A thread-safe implementation of Singleton using a Metaclass.
  """
  _instances = {}
  _lock: threading.Lock = threading.Lock()

  def __call__(cls, *args, **kwargs):
    # First check (no lock) - fast path for when the instance already exists
    if cls not in cls._instances:
      # Acquire lock only if the instance doesn't exist yet
      with cls._lock:
        # Second check (with lock) - ensures another thread didn't create
        # the instance while we were waiting for the lock
        if cls not in cls._instances:
          instance = super().__call__(*args, **kwargs)
          cls._instances[cls] = instance
    return cls._instances[cls]
