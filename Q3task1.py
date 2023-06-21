import requests
import re

class miRBaseVersionFetcher:
    def __init__(self):
        self.url = "https://www.mirbase.org/ftp/CURRENT/README"

    def get_mirbase_version(self):
        response = requests.get(self.url)
        readme_text = response.text
        # Search for the version information in the README file
        version_match = re.search(r"Release (\d+\.\d+)", readme_text)
        if version_match:
            miRBase_version = version_match.group(1)
            return miRBase_version
        else:
            return None

    def print_mirbase_version(self):
        mirbase_version = self.get_mirbase_version()
        if mirbase_version:
            print("Current version of miRBase:", mirbase_version)
        else:
            print("Unable to retrieve the current version of miRBase.")

version_fetcher = miRBaseVersionFetcher()

version_fetcher.print_mirbase_version()
