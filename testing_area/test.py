from selenium import webdriver

def open_webpage(url, driver=None):
    if not driver:
        driver = webdriver.Chrome()  # You can use other browsers by specifying the appropriate webdriver

    driver.get(url)

if __name__ == '__main__':
    # Example usage:
    url1 = "https://example.com/page1"
    url2 = "https://example.com/page2"

    driver = None

    open_webpage(url1, driver)
    # Do something with the first tab

    open_webpage(url2, driver)
    # Continue working with the same tab, now on a different subpage

    # Close the browser when done
    driver.quit()

