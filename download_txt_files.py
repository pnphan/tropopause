import os
import requests
from bs4 import BeautifulSoup
import zipfile
from io import BytesIO

# Directory URL
BASE_URL = "https://www.ncei.noaa.gov/data/integrated-global-radiosonde-archive/access/data-por/"

# Directory to save downloaded files
SAVE_DIR = "downloaded_txt_files"
os.makedirs(SAVE_DIR, exist_ok=True)

def get_zip_file_links(base_url):
    """Fetch all .zip file links from the given URL."""
    response = requests.get(base_url)
    if response.status_code != 200:
        raise Exception(f"Failed to access {base_url}. HTTP Status Code: {response.status_code}")
    
    soup = BeautifulSoup(response.text, "html.parser")
    # Find all links ending with .zip
    zip_links = [base_url + link['href'] for link in soup.find_all('a', href=True) if link['href'].endswith('.zip')]
    return zip_links

def download_and_extract_zip(url, save_dir):
    """Download a .zip file from a URL, extract .txt files, and save them to the specified directory."""
    try:
        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with zipfile.ZipFile(BytesIO(response.content)) as z:
                for file_name in z.namelist():
                    if file_name.endswith('.txt'):
                        extracted_path = os.path.join(save_dir, file_name)
                        os.makedirs(os.path.dirname(extracted_path), exist_ok=True)
                        with open(extracted_path, 'wb') as f:
                            f.write(z.read(file_name))
                        print(f"Extracted: {file_name}")
        else:
            print(f"Failed to download: {url}. HTTP Status Code: {response.status_code}")
    except zipfile.BadZipFile:
        print(f"Failed to extract: {url}. Not a valid zip file.")

def main():
    try:
        print("Fetching .zip file links...")
        zip_links = get_zip_file_links(BASE_URL)
        print(f"Found {len(zip_links)} .zip files. Starting download and extraction...")
        
        for link in zip_links:
            download_and_extract_zip(link, SAVE_DIR)
        
        print("All .txt files extracted successfully.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
