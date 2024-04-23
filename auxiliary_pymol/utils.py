import base64
import os


def create_base64_string_from_file(filepath: str) -> str:
    """Creates a base64 string from a binary file.

    Args:
        filepath: a filepath to a binary file.

    Returns:
        a base64 encoded string.
        an empty string if filepath could not be found.
    """
    if os.path.exists(filepath):
        with open(filepath, "rb") as binary_file:
            binary_data = binary_file.read()
            binary_file.close()
        # Encode binary data to base64 and decode to utf-8 string
        return base64.b64encode(binary_data).decode("utf-8")
    else:
        return ""
