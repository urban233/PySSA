import os
import requests


def download_1pt4_pdb_file() -> str:
    local_filepath = os.path.join(os.path.join(os.path.expanduser("~"), "Downloads"),
                                  'pytest protein kalata family with a lot of spaces.pdb')
    url = 'https://files.rcsb.org/download/1PT4.pdb'
    response = requests.get(url)
    if response.status_code == 200:
        with open(local_filepath, 'wb') as file:
            file.write(response.content)
    else:
        print(f"Failed to download file. Status code: {response.status_code}")
        assert False
    return local_filepath


def download_pdb_file_by_id(a_pdb_id: str, a_basename: str) -> str:
    local_filepath = os.path.join(os.path.join(os.path.expanduser("~"), "Downloads"),
                                  f'{a_basename}.pdb')
    url = f'https://files.rcsb.org/download/{a_pdb_id.upper()}.pdb'
    response = requests.get(url)
    if response.status_code == 200:
        with open(local_filepath, 'wb') as file:
            file.write(response.content)
    else:
        print(f"Failed to download file. Status code: {response.status_code}")
        assert False
    return local_filepath


def download_file_from_url(url: str, filename: str) -> str:
    local_filepath = os.path.join(os.path.join(os.path.expanduser("~"), "Downloads"),
                                  filename)
    response = requests.get(url)
    if response.status_code == 200:
        with open(local_filepath, 'wb') as file:
            file.write(response.content)
    else:
        print(f"Failed to download file. Status code: {response.status_code}")
        assert False
    return local_filepath
