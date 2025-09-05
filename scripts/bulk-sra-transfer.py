import os
import subprocess
import logging
import boto3
from concurrent.futures import ThreadPoolExecutor, as_completed
from botocore.exceptions import ClientError
import hashlib

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

METADATA_BUCKET = "sra-pub-metadata-us-east-1"
METADATA_PREFIX = "sra/metadata/"
DESTINATION_BUCKET = "seqtoid-public-references"
DESTINATION_PREFIX = "test_files/"
S3_CLIENT = boto3.client("s3", region_name="us-east-1")
AWS_PROFILE = "AWSAdministratorAccess-941377154785"


try:
    SESSION = boto3.Session(profile_name=AWS_PROFILE)
    S3_CLIENT = SESSION.client("s3", region_name="us-east-1")
except Exception as e:
    logger.error(f"Failed to initialize boto3 session with profile {AWS_PROFILE}: {e}")
    exit(1)

def refresh_sso_token():
    """Refresh AWS SSO token by running aws sso login."""
    try:
        logger.info(f"Refreshing SSO token for profile {AWS_PROFILE}")
        result = subprocess.run(
            ["aws", "--profile", AWS_PROFILE, "sso", "login"],
            capture_output=True,
            text=True,
            check=True
        )
        logger.info(f"SSO login successful: {result.stdout}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to refresh SSO token: {e.stderr}")
        raise
    except Exception as e:
        logger.error(f"Error during SSO login: {e}")
        raise

def get_metadata_files():
    """List metadata files in s3://sra-pub-metadata-us-east-1/sra/metadata/."""
    try:
        response = S3_CLIENT.list_objects_v2(Bucket=METADATA_BUCKET, Prefix=METADATA_PREFIX)
        return [obj["Key"] for obj in response.get("Contents", []) if obj["Key"].endswith(".csv")]
    except ClientError as e:
        logger.error(f"Failed to list metadata files: {e}")
        raise

def main(sra_names_file):
    with open("SRAnames.txt") as f:
        sra_numbers = [line.strip() for line in f if line.strip()]

    print(sra_numbers)

    metadata_files = get_metadata_files()
    logger.info(f"Found {len(metadata_files)} metadata files")
    print(f"Found {len(metadata_files)} metadata files")

if __name__ == "__main__":

    main("SRAnames.txt")