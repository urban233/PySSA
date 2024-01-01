import pyautogui


class TestPyssa:
    def __init__(self):
        # moves to (519,1060) in 1 sec
        pyautogui.moveTo(1536, 42, duration=1)

        # simulates a click at the present
        # mouse position
        pyautogui.click()
