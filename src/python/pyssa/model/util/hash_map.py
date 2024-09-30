from typing import Optional

import numpy as np
import collections


class HashMap:
  """Implementation of a non automatically resizable hash map."""

  def __init__(self, a_size: int):
    """Constructor.

    Args:
      a_size: The size of the hash map

    Notes:
      Be aware that the hash map is 25% larger than the given size!
    """
    # <editor-fold desc="Checks">
    if a_size == 0 or a_size < 0:
      raise ValueError("'a_size' must have a positive integer value!")
    # </editor-fold>
    self.size = int(a_size + a_size * 0.25)
    self.table: np.ndarray[collections.deque] = np.empty(self.size, dtype=object)
    for i in range(self.size):
      self.table[i] = collections.deque()  # Each index will store a linked list (for chaining)

  def __repr__(self) -> str:
    """Defines the string representation of the hash map."""
    return repr(self.table)

  def _hash_index(self, key) -> int:
    """Compute the hash index for a given key"""
    return hash(key) % self.size

  def insert(self, a_key: object, a_value: object) -> bool:
    """Insert a key-value pair into the hash table"""
    bucket = self.table[self._hash_index(a_key)]
    for i, (tmp_key, tmp_value) in enumerate(bucket):
      if tmp_key == a_key:
        return False
    # If key is not found, append the new key-value pair
    bucket.append((a_key, a_value))
    return True

  def update(self, a_key: object, a_value: object) -> bool:
    """Updates a key-value pair into the hash table """
    bucket = self.table[self._hash_index(a_key)]
    # Check if the key exists in the bucket and update it
    for i, (tmp_key, tmp_value) in enumerate(bucket):
      if tmp_key == a_key:
        bucket[i] = (a_key, a_value)  # Update key-value
        return True
    return False  # If key is not found, return False

  def get(self, key) -> Optional[object]:
    """ Search for a key in the hash table and return its value """
    index = self._hash_index(key)
    bucket = self.table[index]
    for k, v in bucket:
      if k == key:
        return v  # Key found, return value
    return None  # Key not found

  def remove(self, key) -> bool:
    """ Remove a key-value pair from the hash table """
    index = self._hash_index(key)
    bucket = self.table[index]
    for i, (k, v) in enumerate(bucket):
      if k == key:
        del bucket[i]  # Remove the key-value pair
        return True
    return False  # Key not found

  def keys(self):
    """Returns a list of all keys in the hash table"""
    tmp_keys = []
    for tmp_bucket in self.table:
      for tmp_key_value_pair in tmp_bucket:
        tmp_keys.append(tmp_key_value_pair[0])
    return tmp_keys


if __name__ == "__main__":
  """Test main function"""
  tmp_hash_map = HashMap(10)
  tmp_hash_map.insert("test", 11)
  tmp_hash_map.insert("test2", 12)
  tmp_hash_map.keys()

  # print(tmp_hash_map.size)
  # print(tmp_hash_map)
  # print(tmp_hash_map.get("test2"))
  # tmp_hash_map.remove("test")
  # print(tmp_hash_map)
