from PyQt5 import QtWidgets
from PyQt5 import QtCore


class CustomLineEdit(QtWidgets.QLineEdit):
    def keyPressEvent(self, event):
        # Get the key code
        key_text = event.text()

        # Check if the key is allowed, Backspace is allowed, or if it's an empty string (allowing empty input)
        allowed_chars = set("0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_-")
        if key_text not in allowed_chars and key_text != "" and event.key() != QtCore.Qt.Key_Backspace:
            # Ignore the key event
            return

        # Call the base class implementation to handle other keys
        super(CustomLineEdit, self).keyPressEvent(event)
