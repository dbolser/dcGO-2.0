"""Unit tests for the downloader's checksum verification.

Pinning a dataset to an immutable release URL only buys reproducibility if the
bytes are actually checked, so this covers the verification path rather than the
(network-bound) download itself.
"""

import hashlib
import importlib.util
from pathlib import Path

import pytest

DOWNLOADER_PATH = Path(__file__).resolve().parents[2] / "scripts" / "download_data.py"


@pytest.fixture(scope="module")
def downloader():
    """Import the CLI script as a module so its helpers can be exercised."""
    spec = importlib.util.spec_from_file_location("download_data", DOWNLOADER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def payload(tmp_path):
    path = tmp_path / "doid.obo"
    path.write_bytes(b"format-version: 1.2\n")
    return path


class TestVerifyChecksum:
    def test_matching_digest_passes(self, downloader, payload):
        digest = hashlib.sha256(payload.read_bytes()).hexdigest()
        downloader.verify_checksum(payload, ("sha256", digest))  # does not raise

    def test_mismatch_names_both_digests(self, downloader, payload):
        with pytest.raises(downloader.ChecksumError, match="sha256 mismatch"):
            downloader.verify_checksum(payload, ("sha256", "00" * 32))

    def test_streams_in_chunks(self, downloader, tmp_path):
        # Larger than one read block, so the incremental update path is used.
        path = tmp_path / "big.bin"
        blob = b"x" * (downloader.CHUNK_SIZE * 2 + 7)
        path.write_bytes(blob)
        assert (
            downloader.file_digest(path, "sha256") == hashlib.sha256(blob).hexdigest()
        )

    def test_other_algorithms_supported(self, downloader, payload):
        digest = hashlib.md5(payload.read_bytes()).hexdigest()
        downloader.verify_checksum(payload, ("md5", digest))
