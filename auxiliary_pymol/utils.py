import base64
import os
import pathlib


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


def write_binary_file_from_base64_string(filepath: pathlib.Path, base64_data) -> None:  # noqa: ANN001
    """Writes base64 data to a binary file.

    Args:
        filepath: a filepath to a binary file.
        base64_data: a base64 string to write.
    """
    # Decode base64 string to binary data
    binary_data_export = base64.b64decode(base64_data)
    # Ensure that the directory exists
    directory = os.path.dirname(filepath)
    if not os.path.exists(directory):
        os.makedirs(directory)

    with open(filepath, "wb") as f:
        f.write(binary_data_export)
