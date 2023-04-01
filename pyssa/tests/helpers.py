from PyQt5.QtTest import QTest


def check_if_label_changed_in_500_ms(label, text: str):
    for _ in range(50):  # Check 50 times (total of 500 ms wait time)
        QTest.qWait(10)  # Wait for 10 ms
        if label.text() == text:
            break


def change_page(button_for_page, label_page_title, page_title_text):
    # Click the button
    button_for_page.animateClick()
    # Wait for the label text to update
    check_if_label_changed_in_500_ms(label_page_title, page_title_text)
