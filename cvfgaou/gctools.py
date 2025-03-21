"""Some convenience functions for Google Cloud"""

from google.cloud import storage
from google.cloud.storage.blob import Blob

def blob_exists(uri):
    return Blob.from_string(uri, storage.Client()).exists()
