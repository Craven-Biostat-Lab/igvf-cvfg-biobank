"""Some convenience functions for Google Cloud"""

from google.cloud import storage
from google.cloud.storage.blob import Blob

def blob_exists(uri):
    """Check if a blob at a URI exists."""
    return Blob.from_string(uri, storage.Client()).exists()

def list_blobs(uri: str):
    """List blobs given a gs:// URI."""

    if not uri.startswith('gs://'):
        raise ValueError(f'Expected URI "{uri}" to start with, "gs://".')
    
    # The "5" below is the length of "gs://"
    parts = uri[5:].split('/', maxsplit=1)

    # Everything before the first slash should be the bucket
    bucket = storage.Client().get_bucket(parts[0])

    # Everything after the bucket will be treated as blob prefix
    blobs = bucket.list_blobs() if len(parts) == 1 else bucket.list_blobs(prefix=parts[1])

    # Results will be gs:// URIs
    return [f'gs://{b.bucket.name}/{b.name}' for b in blobs]
