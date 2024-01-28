import os
import requests
from concurrent.futures import ThreadPoolExecutor

def download_file(file_url, file_path):
    with requests.get(file_url, stream=True) as response:
        with open(file_path, 'wb') as file:
            for chunk in response.iter_content(chunk_size=128):
                file.write(chunk)

def download_directory(url, local_path, num_threads=5):
    response = requests.get(url)

    if response.status_code == 200:
        if not os.path.exists(local_path):
            os.makedirs(local_path)

        files = response.json()

        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            futures = []

            for file_name in files:
                file_url = f"{url}/{file_name}"
                file_path = os.path.join(local_path, file_name)
                future = executor.submit(download_file, file_url, file_path)
                futures.append(future)

            for future in futures:
                future.result()

        print("Download complete.")
    else:
        print(f"Failed to download directory. Status code: {response.status_code}")
if __name__ == "__main__":
    # Example usage
    server_url = "http://192.168.40.67:8000/local_music/"
    local_directory = "local_directory"

    download_directory(server_url, local_directory)
