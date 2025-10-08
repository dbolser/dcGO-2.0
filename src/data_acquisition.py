"""Data acquisition module for the dcGO pipeline.

This module handles downloading and setting up all required datasets including:
- UniProt sequences (Swiss-Prot and TrEMBL)
- Gene Ontology annotations (GOA)
- GO ontology files
- Pfam HMM profiles
- InterProScan software

The module provides robust downloading with progress tracking, error handling,
and proper setup of external tools like InterProScan. Supports both synchronous
and asynchronous operations using modern Python 3.12 features.
"""

import asyncio
import ftplib
import gzip
import hashlib
import os
import stat
import tarfile
import urllib.parse
from pathlib import Path
from typing import Dict, Optional, List, Tuple, Any
from dataclasses import dataclass
from contextlib import asynccontextmanager, contextmanager
from concurrent.futures import ThreadPoolExecutor
import time

import aiohttp
import aiofiles
import requests
from loguru import logger
from tqdm import tqdm

from config.settings import Config


class DataAcquisitionError(Exception):
    """Exception raised when data acquisition fails."""


@dataclass
class DownloadProgress:
    """Tracks download progress with modern dataclass."""

    source_name: str
    url: str
    filepath: Path
    total_size: Optional[int] = None
    downloaded_size: int = 0
    status: str = "pending"  # pending, downloading, completed, failed
    error: Optional[str] = None
    start_time: Optional[float] = None
    end_time: Optional[float] = None

    @property
    def progress_percentage(self) -> float:
        """Calculate download progress percentage."""
        if self.total_size and self.total_size > 0:
            return min(100.0, (self.downloaded_size / self.total_size) * 100.0)
        return 0.0

    @property
    def download_speed(self) -> float:
        """Calculate download speed in bytes per second."""
        if self.start_time and self.downloaded_size > 0:
            elapsed = time.time() - self.start_time
            return self.downloaded_size / elapsed if elapsed > 0 else 0.0
        return 0.0


class DataAcquisition:
    """Handles acquisition and setup of all required datasets for dcGO pipeline.

    This class provides methods to download large datasets from various sources
    with proper progress tracking, error handling, and verification. It also
    handles the setup of external tools like InterProScan.

    Supports both synchronous and asynchronous operations for efficient parallel
    downloading using Python 3.12 features including modern type hints and
    async context managers.

    Attributes:
        config: Configuration object containing data source URLs and paths
        data_dir: Directory where raw data files are stored
        chunk_size: Size of chunks for streaming downloads (default 64KB)
        timeout: Request timeout in seconds (default 30)
        max_concurrent: Maximum concurrent downloads for async operations
        progress_callbacks: List of progress callback functions
    """

    def __init__(
        self,
        config: Config,
        chunk_size: int = 65536,
        timeout: int = 30,
        max_concurrent: int = 3,
    ) -> None:
        """Initialize data acquisition system.

        Args:
            config: Configuration object with data sources and paths
            chunk_size: Size of chunks for streaming downloads in bytes
            timeout: Request timeout in seconds
            max_concurrent: Maximum concurrent downloads for async operations
        """
        self.config = config
        self.data_dir = config.DATA_DIR / "raw"
        self.chunk_size = chunk_size
        self.timeout = timeout
        self.max_concurrent = max_concurrent
        self.progress_callbacks: List[callable] = []
        self._download_progress: Dict[str, DownloadProgress] = {}

        # Ensure data directory exists
        self.data_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Data acquisition initialized. Data directory: {self.data_dir}")

    def add_progress_callback(self, callback: callable) -> None:
        """Add a progress callback function.

        Args:
            callback: Function that receives DownloadProgress objects
        """
        self.progress_callbacks.append(callback)

    def _notify_progress(self, progress: DownloadProgress) -> None:
        """Notify all registered progress callbacks."""
        for callback in self.progress_callbacks:
            try:
                callback(progress)
            except Exception as e:
                logger.warning(f"Progress callback failed: {e}")

    async def download_all_datasets_async(self) -> Dict[str, Path]:
        """Asynchronously download all required datasets for the dcGO pipeline.

        Downloads datasets concurrently with improved performance and progress tracking.
        Uses async context managers and modern Python 3.12 async features.

        Returns:
            Dictionary mapping source names to local file paths

        Raises:
            DataAcquisitionError: If any critical download fails
        """
        logger.info("Starting async download of all required datasets")
        downloaded_files = {}

        # Create download tasks for all datasets
        download_tasks = []
        for source_name, url in self.config.DATASOURCES.items():
            # Determine local filename
            parsed_url = urllib.parse.urlparse(url)
            filename = Path(parsed_url.path).name
            if not filename:
                filename = f"{source_name}.data"

            filepath = self.data_dir / filename

            # Skip if file already exists and is non-empty
            if filepath.exists() and filepath.stat().st_size > 0:
                logger.info(
                    f"File {filepath} already exists ({filepath.stat().st_size:,} bytes), skipping"
                )
                downloaded_files[source_name] = filepath
                continue

            # Create download task
            task = self._download_async(source_name, url, filepath)
            download_tasks.append(task)

        # Execute downloads with concurrency limit
        if download_tasks:
            semaphore = asyncio.Semaphore(self.max_concurrent)

            async def limited_download(task):
                async with semaphore:
                    return await task

            # Run all downloads concurrently
            results = await asyncio.gather(
                *[limited_download(task) for task in download_tasks],
                return_exceptions=True,
            )

            # Process results
            for i, result in enumerate(results):
                task = download_tasks[i]
                source_name = (
                    task.get_name() if hasattr(task, "get_name") else f"task_{i}"
                )

                if isinstance(result, Exception):
                    logger.error(f"Failed to download {source_name}: {result}")
                    if source_name in [
                        "uniprot_sprot",
                        "goa_annotations",
                        "go_ontology",
                    ]:
                        raise DataAcquisitionError(
                            f"Critical dataset {source_name} failed: {result}"
                        )
                elif isinstance(result, tuple):
                    name, path = result
                    downloaded_files[name] = path

        logger.info(
            f"Async dataset download completed. Downloaded {len(downloaded_files)} files"
        )
        return downloaded_files

    async def _download_async(
        self, source_name: str, url: str, filepath: Path
    ) -> Tuple[str, Path]:
        """Asynchronously download a single file.

        Args:
            source_name: Name of the data source
            url: URL to download from
            filepath: Local path to save the file

        Returns:
            Tuple of source name and downloaded file path

        Raises:
            DataAcquisitionError: If download fails
        """
        progress = DownloadProgress(source_name, url, filepath)
        progress.status = "downloading"
        progress.start_time = time.time()
        self._download_progress[source_name] = progress

        try:
            if url.startswith("ftp://"):
                # FTP downloads still use sync approach due to library limitations
                await asyncio.get_event_loop().run_in_executor(
                    None, self._download_ftp_sync, url, filepath, None
                )
            else:
                await self._download_http_async(url, filepath, progress)

            progress.status = "completed"
            progress.end_time = time.time()
            self._notify_progress(progress)

            logger.info(
                f"Successfully downloaded {filepath} ({filepath.stat().st_size:,} bytes)"
            )
            return source_name, filepath

        except Exception as e:
            progress.status = "failed"
            progress.error = str(e)
            progress.end_time = time.time()
            self._notify_progress(progress)
            raise DataAcquisitionError(f"Failed to download {source_name}: {e}") from e

    async def _download_http_async(
        self, url: str, filepath: Path, progress: DownloadProgress
    ) -> None:
        """Asynchronously download file using HTTP/HTTPS.

        Args:
            url: HTTP(S) URL to download
            filepath: Local path to save file
            progress: Progress tracking object

        Raises:
            DataAcquisitionError: If download fails
        """
        timeout = aiohttp.ClientTimeout(total=None, connect=self.timeout)

        try:
            async with aiohttp.ClientSession(
                timeout=timeout, headers={"User-Agent": "dcGO-Pipeline/1.0"}
            ) as session:
                async with session.get(url) as response:
                    response.raise_for_status()

                    # Get content length for progress tracking
                    content_length = response.headers.get("content-length")
                    if content_length:
                        progress.total_size = int(content_length)

                    # Create progress bar
                    progress_bar = tqdm(
                        total=progress.total_size,
                        unit="B",
                        unit_scale=True,
                        desc=filepath.name,
                        miniters=1,
                    )

                    # Download with async file operations
                    async with aiofiles.open(filepath, "wb") as f:
                        async for chunk in response.content.iter_chunked(
                            self.chunk_size
                        ):
                            await f.write(chunk)
                            progress.downloaded_size += len(chunk)
                            progress_bar.update(len(chunk))
                            self._notify_progress(progress)

                    progress_bar.close()

        except aiohttp.ClientError as e:
            raise DataAcquisitionError(f"Failed to download {url}: {e}") from e
        except IOError as e:
            raise DataAcquisitionError(f"Failed to write file {filepath}: {e}") from e

    def download_with_progress(
        self,
        url: str,
        filepath: Path,
        expected_size: Optional[int] = None,
        verify_checksum: Optional[str] = None,
    ) -> None:
        """Download files with progress tracking and optional verification.

        Args:
            url: URL to download from
            filepath: Local path where file should be saved
            expected_size: Expected file size in bytes for progress tracking
            verify_checksum: Expected MD5 checksum for verification

        Raises:
            DataAcquisitionError: If download fails or verification fails
        """
        logger.info(f"Downloading {url} to {filepath}")

        # Handle different protocols
        if url.startswith("ftp://"):
            self._download_ftp_sync(url, filepath, expected_size)
        else:
            self._download_http(url, filepath, expected_size)

        # Verify checksum if provided
        if verify_checksum:
            if not self._verify_checksum(filepath, verify_checksum):
                filepath.unlink()  # Remove corrupted file
                raise DataAcquisitionError(
                    f"Checksum verification failed for {filepath}"
                )

        logger.info(
            f"Successfully downloaded {filepath} ({filepath.stat().st_size:,} bytes)"
        )

    def _download_http(
        self, url: str, filepath: Path, expected_size: Optional[int] = None
    ) -> None:
        """Download file using HTTP/HTTPS with streaming and progress bar.

        Args:
            url: HTTP(S) URL to download
            filepath: Local path to save file
            expected_size: Expected file size for progress tracking

        Raises:
            DataAcquisitionError: If HTTP download fails
        """
        try:
            response = requests.get(
                url,
                stream=True,
                timeout=self.timeout,
                headers={"User-Agent": "dcGO-Pipeline/1.0"},
            )
            response.raise_for_status()

            # Get file size from headers if not provided
            if expected_size is None:
                content_length = response.headers.get("content-length")
                if content_length:
                    expected_size = int(content_length)

            # Create progress bar
            progress_bar = tqdm(
                total=expected_size,
                unit="B",
                unit_scale=True,
                desc=filepath.name,
                miniters=1,
            )

            # Download with progress tracking
            with open(filepath, "wb") as f:
                for chunk in response.iter_content(chunk_size=self.chunk_size):
                    if chunk:  # Filter out keep-alive chunks
                        f.write(chunk)
                        progress_bar.update(len(chunk))

            progress_bar.close()

        except requests.RequestException as e:
            raise DataAcquisitionError(f"Failed to download {url}: {e}") from e
        except IOError as e:
            raise DataAcquisitionError(f"Failed to write file {filepath}: {e}") from e

    def _download_ftp_sync(
        self, url: str, filepath: Path, expected_size: Optional[int] = None
    ) -> None:
        """Download file using FTP with progress tracking.

        Args:
            url: FTP URL to download
            filepath: Local path to save file
            expected_size: Expected file size for progress tracking

        Raises:
            DataAcquisitionError: If FTP download fails
        """
        try:
            parsed_url = urllib.parse.urlparse(url)
            hostname = parsed_url.hostname
            username = parsed_url.username or "anonymous"
            password = parsed_url.password or "anonymous@"
            remote_path = parsed_url.path

            # Connect to FTP server
            with ftplib.FTP(hostname, timeout=self.timeout) as ftp:
                ftp.login(username, password)
                ftp.set_pasv(True)  # Use passive mode

                # Get file size if not provided
                if expected_size is None:
                    try:
                        expected_size = ftp.size(remote_path)
                    except ftplib.error_perm:
                        # SIZE command not supported, continue without size
                        pass

                # Create progress bar
                progress_bar = tqdm(
                    total=expected_size,
                    unit="B",
                    unit_scale=True,
                    desc=filepath.name,
                    miniters=1,
                )

                # Download with progress tracking
                with open(filepath, "wb") as f:

                    def callback(data: bytes) -> None:
                        f.write(data)
                        progress_bar.update(len(data))

                    ftp.retrbinary(
                        f"RETR {remote_path}", callback, blocksize=self.chunk_size
                    )

                progress_bar.close()

        except ftplib.all_errors as e:
            raise DataAcquisitionError(f"Failed to download {url} via FTP: {e}") from e
        except IOError as e:
            raise DataAcquisitionError(f"Failed to write file {filepath}: {e}") from e

    def _verify_checksum(self, filepath: Path, expected_md5: str) -> bool:
        """Verify file integrity using MD5 checksum.

        Args:
            filepath: Path to file to verify
            expected_md5: Expected MD5 hash in hexadecimal format

        Returns:
            True if checksum matches, False otherwise
        """
        logger.info(f"Verifying checksum for {filepath}")

        hash_md5 = hashlib.md5()
        try:
            with open(filepath, "rb") as f:
                for chunk in iter(lambda: f.read(self.chunk_size), b""):
                    hash_md5.update(chunk)

            calculated_md5 = hash_md5.hexdigest()
            matches = calculated_md5.lower() == expected_md5.lower()

            if matches:
                logger.info(f"Checksum verification passed for {filepath}")
            else:
                logger.error(
                    f"Checksum mismatch for {filepath}: expected {expected_md5}, got {calculated_md5}"
                )

            return matches

        except IOError as e:
            logger.error(f"Failed to read file for checksum verification: {e}")
            return False

    def download_all_datasets(self) -> Dict[str, Path]:
        """Download all required datasets for the dcGO pipeline.

        Downloads UniProt sequences, GO annotations, GO ontology,
        Pfam HMMs, and InterProScan software. Skips files that already exist
        unless they appear corrupted or incomplete.

        Returns:
            Dictionary mapping source names to local file paths

        Raises:
            DataAcquisitionError: If any critical download fails
        """
        logger.info("Starting download of all required datasets")
        downloaded_files = {}

        for source_name, url in self.config.DATASOURCES.items():
            try:
                # Determine local filename
                parsed_url = urllib.parse.urlparse(url)
                filename = Path(parsed_url.path).name
                if not filename:
                    # Fallback for URLs without clear filenames
                    filename = f"{source_name}.data"

                filepath = self.data_dir / filename

                # Check if file already exists and is non-empty
                if filepath.exists() and filepath.stat().st_size > 0:
                    logger.info(
                        f"File {filepath} already exists ({filepath.stat().st_size:,} bytes), skipping download"
                    )
                else:
                    # Download the file
                    self.download_with_progress(url, filepath)

                downloaded_files[source_name] = filepath

            except DataAcquisitionError as e:
                logger.error(f"Failed to download {source_name} from {url}: {e}")
                # For critical datasets, re-raise the error
                if source_name in ["uniprot_sprot", "goa_annotations", "go_ontology"]:
                    raise
                # For optional datasets, continue with warning
                logger.warning(f"Continuing without {source_name} dataset")

        logger.info(
            f"Dataset download completed. Downloaded {len(downloaded_files)} files"
        )
        return downloaded_files

    async def setup_interproscan_async(self, interpro_archive: Path) -> Path:
        """Asynchronously extract and configure InterProScan for domain scanning.

        Uses async file operations and threading for CPU-intensive extraction.
        Implements modern Python 3.12 patterns with proper error handling.

        Args:
            interpro_archive: Path to InterProScan tar.gz archive

        Returns:
            Path to the main InterProScan executable script

        Raises:
            DataAcquisitionError: If extraction or setup fails
        """
        logger.info(f"Setting up InterProScan from {interpro_archive}")

        extract_dir = self.data_dir / "interproscan"

        # Check if already extracted
        if extract_dir.exists():
            potential_paths = list(extract_dir.glob("interproscan-*/interproscan.sh"))
            if potential_paths:
                interpro_script = potential_paths[0]
                if interpro_script.exists():
                    logger.info(f"InterProScan already extracted at {interpro_script}")
                    await self._configure_interproscan_async(interpro_script)
                    return interpro_script

        # Extract the archive in a thread pool
        try:
            logger.info("Extracting InterProScan archive...")
            extract_dir.mkdir(exist_ok=True)

            # Run extraction in thread pool to avoid blocking
            loop = asyncio.get_event_loop()
            with ThreadPoolExecutor(max_workers=1) as executor:
                interpro_script = await loop.run_in_executor(
                    executor,
                    self._extract_interproscan_sync,
                    interpro_archive,
                    extract_dir,
                )

            # Configure InterProScan asynchronously
            await self._configure_interproscan_async(interpro_script)

            return interpro_script

        except Exception as e:
            raise DataAcquisitionError(f"Failed to setup InterProScan: {e}") from e

    def _extract_interproscan_sync(
        self, interpro_archive: Path, extract_dir: Path
    ) -> Path:
        """Synchronous InterProScan extraction for thread pool execution."""
        with tarfile.open(interpro_archive, "r:gz") as tar:
            members = tar.getmembers()
            total_members = len(members)

            with tqdm(total=total_members, desc="Extracting", unit="files") as pbar:
                for member in members:
                    tar.extract(member, extract_dir)
                    pbar.update(1)

        logger.info("InterProScan archive extracted successfully")

        # Find the main script
        interpro_scripts = list(extract_dir.glob("interproscan-*/interproscan.sh"))
        if not interpro_scripts:
            raise DataAcquisitionError("InterProScan script not found after extraction")

        return interpro_scripts[0]

    async def _configure_interproscan_async(self, interpro_script: Path) -> None:
        """Asynchronously configure InterProScan installation."""
        try:
            # Make scripts executable
            await asyncio.get_event_loop().run_in_executor(
                None, self._make_executable, interpro_script
            )

            # Validate installation structure
            interpro_dir = interpro_script.parent
            required_dirs = ["bin", "data", "lib"]
            missing_dirs = [d for d in required_dirs if not (interpro_dir / d).exists()]

            if missing_dirs:
                logger.warning(f"Some expected directories are missing: {missing_dirs}")
            else:
                logger.info("InterProScan directory structure looks complete")

            # Set up temporary directory
            temp_dir = interpro_dir / "temp"
            temp_dir.mkdir(exist_ok=True)

            # Set environment variables
            os.environ["JAVA_TOOL_OPTIONS"] = f"-Djava.io.tmpdir={temp_dir}"

            logger.info("InterProScan async configuration completed successfully")

        except Exception as e:
            raise DataAcquisitionError(f"Failed to configure InterProScan: {e}") from e

    def _make_executable(self, script_path: Path) -> None:
        """Make script and related files executable."""
        # Make main script executable
        current_permissions = script_path.stat().st_mode
        script_path.chmod(current_permissions | stat.S_IEXEC)
        logger.info(f"Made {script_path} executable")

        # Make other scripts executable
        interpro_dir = script_path.parent
        for script_pattern in ["*.sh", "bin/*"]:
            for script_file in interpro_dir.glob(script_pattern):
                if script_file.is_file():
                    current_perms = script_file.stat().st_mode
                    script_file.chmod(current_perms | stat.S_IEXEC)

    def setup_interproscan(self, interpro_archive: Path) -> Path:
        """Extract and configure InterProScan for domain scanning.

        Extracts the InterProScan archive, sets appropriate permissions,
        and prepares the installation for use. This includes making the
        main script executable and setting up any required configuration.

        Args:
            interpro_archive: Path to InterProScan tar.gz archive

        Returns:
            Path to the main InterProScan executable script

        Raises:
            DataAcquisitionError: If extraction or setup fails
        """
        logger.info(f"Setting up InterProScan from {interpro_archive}")

        extract_dir = self.data_dir / "interproscan"

        # Check if already extracted
        if extract_dir.exists():
            # Look for existing installation
            potential_paths = list(extract_dir.glob("interproscan-*/interproscan.sh"))
            if potential_paths:
                interpro_script = potential_paths[0]
                if interpro_script.exists():
                    logger.info(f"InterProScan already extracted at {interpro_script}")
                    self._configure_interproscan(interpro_script)
                    return interpro_script

        # Extract the archive
        try:
            interpro_script = self._extract_interproscan_sync(
                interpro_archive, extract_dir
            )
            self._configure_interproscan(interpro_script)
            return interpro_script

        except Exception as e:
            raise DataAcquisitionError(
                f"Failed to extract InterProScan archive: {e}"
            ) from e

    def _configure_interproscan(self, interpro_script: Path) -> None:
        """Configure InterProScan installation.

        Sets proper file permissions, validates the installation,
        and performs any necessary configuration steps.

        Args:
            interpro_script: Path to the main InterProScan script

        Raises:
            DataAcquisitionError: If configuration fails
        """
        try:
            # Make main script executable
            current_permissions = interpro_script.stat().st_mode
            interpro_script.chmod(current_permissions | stat.S_IEXEC)
            logger.info(f"Made {interpro_script} executable")

            # Make other scripts executable if they exist
            interpro_dir = interpro_script.parent
            for script_pattern in ["*.sh", "bin/*"]:
                for script_path in interpro_dir.glob(script_pattern):
                    if script_path.is_file():
                        current_perms = script_path.stat().st_mode
                        script_path.chmod(current_perms | stat.S_IEXEC)

            # Check for required subdirectories
            required_dirs = ["bin", "data", "lib"]
            missing_dirs = [d for d in required_dirs if not (interpro_dir / d).exists()]

            if missing_dirs:
                logger.warning(f"Some expected directories are missing: {missing_dirs}")
            else:
                logger.info("InterProScan directory structure looks complete")

            # Set up temporary directory if needed
            temp_dir = interpro_dir / "temp"
            temp_dir.mkdir(exist_ok=True)

            # Set environment variable for Java temp directory
            os.environ["JAVA_TOOL_OPTIONS"] = f"-Djava.io.tmpdir={temp_dir}"

            logger.info("InterProScan configuration completed successfully")

        except (OSError, IOError) as e:
            raise DataAcquisitionError(f"Failed to configure InterProScan: {e}") from e

    def download_specific_dataset(
        self,
        source_name: str,
        custom_url: Optional[str] = None,
        force_redownload: bool = False,
    ) -> Path:
        """Download a specific dataset by name.

        Args:
            source_name: Name of the data source to download
            custom_url: Custom URL to use instead of default
            force_redownload: Whether to redownload even if file exists

        Returns:
            Path to the downloaded file

        Raises:
            DataAcquisitionError: If the source is unknown or download fails
            ValueError: If source_name is not recognized
        """
        if source_name not in self.config.DATASOURCES and not custom_url:
            raise ValueError(f"Unknown data source: {source_name}")

        url = custom_url or self.config.DATASOURCES[source_name]

        # Preserve original filename and organize by source subdirectory
        # Extract the original filename from URL
        parsed_url = urllib.parse.urlparse(url)
        original_filename = Path(parsed_url.path).name

        # Create source-specific subdirectory
        source_dir = self.data_dir / source_name
        source_dir.mkdir(parents=True, exist_ok=True)

        filepath = source_dir / original_filename

        # Check if redownload is needed
        if filepath.exists() and not force_redownload:
            logger.info(
                f"File {filepath} already exists, use force_redownload=True to override"
            )
            return filepath

        # Download the file
        logger.info(f"Downloading {source_name} dataset")
        logger.info(f"  Source: {url}")
        logger.info(f"  Destination: {filepath}")
        self.download_with_progress(url, filepath)

        return filepath

    def verify_all_downloads(self) -> Dict[str, bool]:
        """Verify integrity of all downloaded files.

        Performs basic integrity checks on downloaded files including
        size validation and format verification where applicable.

        Returns:
            Dictionary mapping source names to verification status
        """
        logger.info("Verifying downloaded datasets")
        verification_results = {}

        for source_name in self.config.DATASOURCES:
            # Find the downloaded file
            potential_files = list(self.data_dir.glob(f"{source_name}*"))

            if not potential_files:
                verification_results[source_name] = False
                logger.warning(f"No file found for {source_name}")
                continue

            filepath = potential_files[0]

            # Basic checks
            if not filepath.exists():
                verification_results[source_name] = False
                continue

            if filepath.stat().st_size == 0:
                verification_results[source_name] = False
                logger.warning(f"File {filepath} is empty")
                continue

            # Format-specific checks
            try:
                if filepath.suffix == ".gz":
                    # Test gzip files can be opened
                    with gzip.open(filepath, "rt") as f:
                        f.readline()  # Try to read first line
                elif filepath.suffix == ".obo":
                    # Basic OBO file validation
                    with open(filepath, "r") as f:
                        first_line = f.readline().strip()
                        if not first_line.startswith("format-version:"):
                            logger.warning(f"OBO file {filepath} may be corrupted")
                            verification_results[source_name] = False
                            continue

                verification_results[source_name] = True
                logger.info(f"Verification passed for {source_name}")

            except Exception as e:
                logger.warning(f"Verification failed for {source_name}: {e}")
                verification_results[source_name] = False

        passed = sum(verification_results.values())
        total = len(verification_results)
        logger.info(f"Dataset verification completed: {passed}/{total} files passed")

        return verification_results

    def get_download_summary(self) -> Dict[str, dict]:
        """Get summary information about downloaded datasets.

        Returns:
            Dictionary with information about each downloaded file
        """
        summary = {}

        for source_name in self.config.DATASOURCES:
            potential_files = list(self.data_dir.glob(f"{source_name}*"))

            if potential_files:
                filepath = potential_files[0]
                file_stats = filepath.stat()
                summary[source_name] = {
                    "path": str(filepath),
                    "size_bytes": file_stats.st_size,
                    "size_mb": round(file_stats.st_size / 1024 / 1024, 2),
                    "modified": file_stats.st_mtime,
                    "exists": True,
                }
            else:
                summary[source_name] = {
                    "path": None,
                    "size_bytes": 0,
                    "size_mb": 0,
                    "modified": None,
                    "exists": False,
                }

        return summary

    @contextmanager
    def download_session(self):
        """Context manager for download sessions with proper cleanup."""
        session_start = time.time()
        logger.info("Starting download session")

        try:
            yield self
        except Exception as e:
            logger.error(f"Download session failed: {e}")
            raise
        finally:
            session_duration = time.time() - session_start
            logger.info(f"Download session completed in {session_duration:.2f} seconds")

    @asynccontextmanager
    async def async_download_session(self):
        """Async context manager for download sessions with proper cleanup."""
        session_start = time.time()
        logger.info("Starting async download session")

        try:
            yield self
        except Exception as e:
            logger.error(f"Async download session failed: {e}")
            raise
        finally:
            session_duration = time.time() - session_start
            logger.info(
                f"Async download session completed in {session_duration:.2f} seconds"
            )

    async def verify_all_downloads_async(self) -> Dict[str, bool]:
        """Asynchronously verify integrity of all downloaded files.

        Uses thread pool for I/O operations to avoid blocking the event loop.
        Implements modern Python 3.12 async patterns.

        Returns:
            Dictionary mapping source names to verification status
        """
        logger.info("Verifying downloaded datasets asynchronously")
        verification_tasks = []

        # Create verification tasks
        for source_name in self.config.DATASOURCES:
            potential_files = list(self.data_dir.glob(f"{source_name}*"))
            if potential_files:
                filepath = potential_files[0]
                task = self._verify_file_async(source_name, filepath)
                verification_tasks.append(task)

        # Run verifications concurrently
        if verification_tasks:
            results = await asyncio.gather(*verification_tasks, return_exceptions=True)
            verification_results = {}

            for i, result in enumerate(results):
                if isinstance(result, tuple):
                    source_name, status = result
                    verification_results[source_name] = status
                elif isinstance(result, Exception):
                    logger.error(f"Verification task failed: {result}")
                    # Handle failure case
                    source_name = list(self.config.DATASOURCES.keys())[i]
                    verification_results[source_name] = False
        else:
            verification_results = {name: False for name in self.config.DATASOURCES}

        passed = sum(verification_results.values())
        total = len(verification_results)
        logger.info(
            f"Async dataset verification completed: {passed}/{total} files passed"
        )

        return verification_results

    async def _verify_file_async(
        self, source_name: str, filepath: Path
    ) -> Tuple[str, bool]:
        """Asynchronously verify a single file."""
        try:
            # Run verification in thread pool
            loop = asyncio.get_event_loop()
            with ThreadPoolExecutor(max_workers=1) as executor:
                is_valid = await loop.run_in_executor(
                    executor, self._verify_single_file, filepath
                )

            if is_valid:
                logger.info(f"Verification passed for {source_name}")
            else:
                logger.warning(f"Verification failed for {source_name}")

            return source_name, is_valid

        except Exception as e:
            logger.warning(f"Verification failed for {source_name}: {e}")
            return source_name, False

    def _verify_single_file(self, filepath: Path) -> bool:
        """Verify a single file integrity."""
        # Basic checks
        if not filepath.exists() or filepath.stat().st_size == 0:
            return False

        # Format-specific checks
        try:
            if filepath.suffix == ".gz":
                with gzip.open(filepath, "rt") as f:
                    f.readline()  # Try to read first line
            elif filepath.suffix == ".obo":
                with open(filepath, "r") as f:
                    first_line = f.readline().strip()
                    if not first_line.startswith("format-version:"):
                        return False

            return True

        except Exception:
            return False

    async def download_with_retry_async(
        self,
        url: str,
        filepath: Path,
        max_retries: int = 3,
        retry_delay: float = 1.0,
        backoff_factor: float = 2.0,
    ) -> None:
        """Download with automatic retry and exponential backoff.

        Uses modern Python 3.12 async features with comprehensive error handling.

        Args:
            url: URL to download from
            filepath: Local path to save file
            max_retries: Maximum number of retry attempts
            retry_delay: Initial delay between retries in seconds
            backoff_factor: Multiplier for delay after each retry

        Raises:
            DataAcquisitionError: If all retry attempts fail
        """
        last_exception = None

        for attempt in range(max_retries + 1):
            try:
                if url.startswith("ftp://"):
                    await asyncio.get_event_loop().run_in_executor(
                        None, self._download_ftp_sync, url, filepath, None
                    )
                else:
                    progress = DownloadProgress("retry_download", url, filepath)
                    await self._download_http_async(url, filepath, progress)

                logger.info(
                    f"Successfully downloaded {filepath} on attempt {attempt + 1}"
                )
                return

            except Exception as e:
                last_exception = e

                if attempt < max_retries:
                    delay = retry_delay * (backoff_factor**attempt)
                    logger.warning(
                        f"Download attempt {attempt + 1} failed: {e}. Retrying in {delay:.1f}s"
                    )
                    await asyncio.sleep(delay)
                else:
                    logger.error(f"All {max_retries + 1} download attempts failed")

        raise DataAcquisitionError(
            f"Failed to download {url} after {max_retries + 1} attempts: {last_exception}"
        )

    def get_progress_summary(self) -> Dict[str, Dict[str, Any]]:
        """Get current progress summary for all ongoing downloads.

        Returns:
            Dictionary with progress information for each download
        """
        summary = {}

        for source_name, progress in self._download_progress.items():
            summary[source_name] = {
                "status": progress.status,
                "progress_percentage": progress.progress_percentage,
                "download_speed_mbps": progress.download_speed / (1024 * 1024),
                "total_size_mb": (progress.total_size / (1024 * 1024))
                if progress.total_size
                else None,
                "downloaded_mb": progress.downloaded_size / (1024 * 1024),
                "error": progress.error,
                "url": progress.url,
                "filepath": str(progress.filepath),
            }

        return summary

    def cleanup_old_downloads(self, keep_latest: int = 3) -> None:
        """Clean up old downloaded files to save space.

        Args:
            keep_latest: Number of latest versions to keep for each dataset
        """
        logger.info(f"Cleaning up old downloads, keeping latest {keep_latest} versions")

        for source_name in self.config.DATASOURCES:
            pattern = f"{source_name}*"
            matching_files = list(self.data_dir.glob(pattern))

            if len(matching_files) > keep_latest:
                # Sort by modification time, newest first
                matching_files.sort(key=lambda p: p.stat().st_mtime, reverse=True)

                # Remove old files
                for old_file in matching_files[keep_latest:]:
                    try:
                        old_file.unlink()
                        logger.info(f"Removed old file: {old_file}")
                    except OSError as e:
                        logger.warning(f"Failed to remove {old_file}: {e}")
