import os
import requests
from concurrent.futures import ThreadPoolExecutor

def download_file(url):
    filename = url.split('/')[-1]
    filepath = os.path.join(snakemake.output[0], filename)
    # create the filepath if it doesn't exist
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        with open(filepath, 'wb') as file:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    file.write(chunk)
    except Exception as e:
        print(f"Error downloading {url}: {e}")

def main():
    urls_file = snakemake.input[0]  # Path to the file containing the list of URLs
    nthreads = snakemake.threads  # Number of threads for parallelization

    # Create the output directory if it doesn't exist
    os.makedirs(snakemake.output[0], exist_ok=True)

    # Read the list of URLs from the file
    with open(urls_file, 'r') as file:
        urls = file.read().splitlines()

    # Download files using parallel processing
    with ThreadPoolExecutor(max_workers=nthreads) as executor:
        executor.map(download_file, urls)

    print("Download completed!")

if __name__ == '__main__':
    main()