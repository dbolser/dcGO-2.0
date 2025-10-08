"""
dcGO Database Management System

This module provides comprehensive database management functionality for the
domain-centric Gene Ontology (dcGO) pipeline, including SQLite integration,
data storage, export capabilities, and statistical reporting.

Author: dcGO Pipeline Team
License: MIT
"""

import sqlite3
import pandas as pd
from pathlib import Path
from typing import Dict, List, Any, Optional, Union
from contextlib import contextmanager
from loguru import logger
import json
from datetime import datetime
import threading


class DatabaseError(Exception):
    """Custom exception for database operations"""

    pass


class dcGODatabaseManager:
    """
    Production-ready database management system for dcGO pipeline.

    Provides comprehensive functionality for:
    - Database schema initialization with performance indices
    - Atomic transaction handling with rollback support
    - Annotation storage with metadata tracking
    - TSV export with customizable filtering
    - Statistical reporting and summary generation
    - Connection pooling and thread safety
    """

    def __init__(self, db_path: Union[str, Path]) -> None:
        """
        Initialize database manager with SQLite backend.

        Args:
            db_path: Path to SQLite database file

        Raises:
            DatabaseError: If database initialization fails
        """
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()

        logger.info(f"Initializing dcGO database at {self.db_path}")

        try:
            self._initialize_database()
            self._validate_schema()
        except Exception as e:
            logger.error(f"Failed to initialize database: {e}")
            raise DatabaseError(f"Database initialization failed: {e}")

    @contextmanager
    def get_connection(self):
        """
        Thread-safe database connection context manager.

        Yields:
            sqlite3.Connection: Database connection with proper error handling
        """
        conn = None
        try:
            with self._lock:
                conn = sqlite3.connect(
                    self.db_path,
                    timeout=30.0,  # 30 second timeout
                    isolation_level=None,  # Autocommit mode
                    check_same_thread=False,
                )
                conn.execute("PRAGMA journal_mode=WAL")  # Write-Ahead Logging
                conn.execute("PRAGMA synchronous=NORMAL")  # Performance optimization
                conn.execute("PRAGMA temp_store=MEMORY")  # Use memory for temp tables
                conn.execute("PRAGMA mmap_size=268435456")  # 256MB memory map

            yield conn

        except sqlite3.Error as e:
            logger.error(f"Database connection error: {e}")
            if conn:
                conn.rollback()
            raise DatabaseError(f"Database connection failed: {e}")
        finally:
            if conn:
                conn.close()

    def _initialize_database(self) -> None:
        """
        Initialize SQLite database with optimized schema and indices.

        Creates tables:
        - domain_annotations: Core annotation data with performance indices
        - association_statistics: Statistical test results and contingency data
        - pipeline_metadata: Pipeline execution metadata and configuration
        """
        logger.info("Creating database schema with performance optimizations")

        with self.get_connection() as conn:
            cursor = conn.cursor()

            # Main annotations table with comprehensive schema
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS domain_annotations (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                domain_id TEXT NOT NULL CHECK(length(domain_id) > 0),
                go_id TEXT NOT NULL CHECK(go_id GLOB 'GO:*'),
                fdr_q_value REAL NOT NULL CHECK(fdr_q_value >= 0 AND fdr_q_value <= 1),
                association_score REAL NOT NULL CHECK(association_score >= 0 AND association_score <= 100),
                annotation_type TEXT NOT NULL CHECK(annotation_type IN ('direct', 'propagated')),
                direct_source_term TEXT NOT NULL CHECK(direct_source_term GLOB 'GO:*'),
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                UNIQUE(domain_id, go_id, annotation_type)
            )
            """)

            # Performance indices for fast querying
            cursor.execute(
                "CREATE INDEX IF NOT EXISTS idx_domain_id ON domain_annotations(domain_id)"
            )
            cursor.execute(
                "CREATE INDEX IF NOT EXISTS idx_go_id ON domain_annotations(go_id)"
            )
            cursor.execute(
                "CREATE INDEX IF NOT EXISTS idx_annotation_type ON domain_annotations(annotation_type)"
            )
            cursor.execute(
                "CREATE INDEX IF NOT EXISTS idx_fdr_q_value ON domain_annotations(fdr_q_value)"
            )
            cursor.execute(
                "CREATE INDEX IF NOT EXISTS idx_association_score ON domain_annotations(association_score)"
            )
            cursor.execute(
                "CREATE INDEX IF NOT EXISTS idx_direct_source ON domain_annotations(direct_source_term)"
            )
            cursor.execute(
                "CREATE INDEX IF NOT EXISTS idx_domain_type ON domain_annotations(domain_id, annotation_type)"
            )

            # Statistical data table for detailed analysis
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS association_statistics (
                domain_id TEXT NOT NULL,
                go_id TEXT NOT NULL,
                a_count INTEGER NOT NULL CHECK(a_count >= 0),
                b_count INTEGER NOT NULL CHECK(b_count >= 0),
                c_count INTEGER NOT NULL CHECK(c_count >= 0),
                d_count INTEGER NOT NULL CHECK(d_count >= 0),
                odds_ratio REAL NOT NULL CHECK(odds_ratio >= 0),
                p_value REAL NOT NULL CHECK(p_value >= 0 AND p_value <= 1),
                hypergeometric_score REAL CHECK(hypergeometric_score >= 0 AND hypergeometric_score <= 100),
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                PRIMARY KEY (domain_id, go_id),
                FOREIGN KEY (domain_id, go_id) REFERENCES domain_annotations(domain_id, go_id)
            )
            """)

            # Pipeline metadata and configuration tracking
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS pipeline_metadata (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL,
                data_type TEXT NOT NULL DEFAULT 'string',
                description TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
            """)

            # Summary statistics view for quick access
            cursor.execute("""
            CREATE VIEW IF NOT EXISTS annotation_summary AS
            SELECT 
                annotation_type,
                COUNT(*) as total_count,
                COUNT(DISTINCT domain_id) as unique_domains,
                COUNT(DISTINCT go_id) as unique_go_terms,
                AVG(association_score) as avg_score,
                MIN(fdr_q_value) as min_fdr,
                MAX(fdr_q_value) as max_fdr
            FROM domain_annotations 
            GROUP BY annotation_type
            """)

            conn.commit()
            logger.info("Database schema initialized successfully")

    def _validate_schema(self) -> None:
        """Validate database schema integrity."""
        with self.get_connection() as conn:
            cursor = conn.cursor()

            # Check if required tables exist
            cursor.execute("""
            SELECT name FROM sqlite_master 
            WHERE type='table' AND name IN ('domain_annotations', 'association_statistics', 'pipeline_metadata')
            """)

            tables = {row[0] for row in cursor.fetchall()}
            required_tables = {
                "domain_annotations",
                "association_statistics",
                "pipeline_metadata",
            }

            if not required_tables.issubset(tables):
                missing = required_tables - tables
                raise DatabaseError(f"Missing required tables: {missing}")

            logger.info("Database schema validation passed")

    def store_annotations(
        self,
        annotations: List[Any],
        metadata: Optional[Dict[str, Any]] = None,
        statistics: Optional[List[Dict[str, Any]]] = None,
    ) -> None:
        """
        Store domain-GO annotations with full transaction support.

        Args:
            annotations: List of annotation objects with domain, go_term, q_value, etc.
            metadata: Optional pipeline metadata dictionary
            statistics: Optional statistical test results for detailed analysis

        Raises:
            DatabaseError: If storage operation fails
        """
        if not annotations:
            logger.warning("No annotations provided for storage")
            return

        logger.info(f"Storing {len(annotations)} annotations in database")

        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()

                # Begin explicit transaction
                cursor.execute("BEGIN IMMEDIATE")

                try:
                    # Prepare annotation data
                    annotation_data = []
                    for ann in annotations:
                        annotation_data.append(
                            {
                                "domain_id": str(ann.domain),
                                "go_id": str(ann.go_term),
                                "fdr_q_value": float(ann.q_value),
                                "association_score": float(ann.association_score),
                                "annotation_type": str(ann.annotation_type),
                                "direct_source_term": str(ann.direct_source_term),
                            }
                        )

                    # Clear existing annotations for this run
                    cursor.execute("DELETE FROM domain_annotations")
                    logger.info("Cleared existing annotations")

                    # Insert new annotations using executemany for performance
                    insert_query = """
                    INSERT INTO domain_annotations (
                        domain_id, go_id, fdr_q_value, association_score, 
                        annotation_type, direct_source_term
                    ) VALUES (?, ?, ?, ?, ?, ?)
                    """

                    annotation_tuples = [
                        (
                            ann["domain_id"],
                            ann["go_id"],
                            ann["fdr_q_value"],
                            ann["association_score"],
                            ann["annotation_type"],
                            ann["direct_source_term"],
                        )
                        for ann in annotation_data
                    ]

                    cursor.executemany(insert_query, annotation_tuples)
                    logger.info(f"Inserted {len(annotation_tuples)} annotations")

                    # Store statistical data if provided
                    if statistics:
                        cursor.execute("DELETE FROM association_statistics")

                        stats_query = """
                        INSERT INTO association_statistics (
                            domain_id, go_id, a_count, b_count, c_count, d_count,
                            odds_ratio, p_value, hypergeometric_score
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """

                        stats_tuples = [
                            (
                                stat["domain"],
                                stat["go_term"],
                                stat["a"],
                                stat["b"],
                                stat["c"],
                                stat["d"],
                                stat["odds_ratio"],
                                stat["p_value"],
                                stat.get("hyper_score", 0),
                            )
                            for stat in statistics
                        ]

                        cursor.executemany(stats_query, stats_tuples)
                        logger.info(f"Inserted {len(stats_tuples)} statistical records")

                    # Store metadata if provided
                    if metadata:
                        self._store_metadata(cursor, metadata)

                    # Commit transaction
                    cursor.execute("COMMIT")
                    logger.info("All data committed successfully")

                except Exception as e:
                    cursor.execute("ROLLBACK")
                    logger.error(f"Transaction rolled back due to error: {e}")
                    raise

        except Exception as e:
            logger.error(f"Failed to store annotations: {e}")
            raise DatabaseError(f"Annotation storage failed: {e}")

    def _store_metadata(self, cursor: sqlite3.Cursor, metadata: Dict[str, Any]) -> None:
        """Store pipeline metadata with type information."""
        # Clear existing metadata
        cursor.execute("DELETE FROM pipeline_metadata")

        for key, value in metadata.items():
            # Determine data type
            if isinstance(value, (int, float)):
                data_type = "numeric"
                value_str = str(value)
            elif isinstance(value, bool):
                data_type = "boolean"
                value_str = str(value)
            elif isinstance(value, (list, dict)):
                data_type = "json"
                value_str = json.dumps(value)
            else:
                data_type = "string"
                value_str = str(value)

            cursor.execute(
                """
            INSERT INTO pipeline_metadata (key, value, data_type, updated_at)
            VALUES (?, ?, ?, ?)
            """,
                (key, value_str, data_type, datetime.now().isoformat()),
            )

        logger.info(f"Stored {len(metadata)} metadata entries")

    def export_tsv(
        self,
        output_path: Union[str, Path],
        annotation_type: Optional[str] = None,
        min_score: Optional[float] = None,
        max_fdr: Optional[float] = None,
        domains: Optional[List[str]] = None,
    ) -> Path:
        """
        Export annotations to TSV format with flexible filtering.

        Args:
            output_path: Output file path
            annotation_type: Filter by annotation type ('direct' or 'propagated')
            min_score: Minimum association score threshold
            max_fdr: Maximum FDR q-value threshold
            domains: List of specific domains to export

        Returns:
            Path to exported TSV file

        Raises:
            DatabaseError: If export operation fails
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        logger.info(f"Exporting annotations to {output_path}")

        try:
            with self.get_connection() as conn:
                # Build dynamic query with filters
                query_parts = [
                    "SELECT domain_id, go_id, fdr_q_value, association_score,",
                    "       annotation_type, direct_source_term, created_at",
                    "FROM domain_annotations",
                    "WHERE 1=1",
                ]
                params = []

                if annotation_type:
                    query_parts.append("AND annotation_type = ?")
                    params.append(annotation_type)

                if min_score is not None:
                    query_parts.append("AND association_score >= ?")
                    params.append(min_score)

                if max_fdr is not None:
                    query_parts.append("AND fdr_q_value <= ?")
                    params.append(max_fdr)

                if domains:
                    placeholders = ",".join(["?" for _ in domains])
                    query_parts.append(f"AND domain_id IN ({placeholders})")
                    params.extend(domains)

                query_parts.append(
                    "ORDER BY domain_id, annotation_type, association_score DESC"
                )

                query = " ".join(query_parts)

                # Execute query and export
                df = pd.read_sql_query(query, conn, params=params)

                if df.empty:
                    logger.warning("No annotations match the specified filters")
                    return output_path

                # Export to TSV with proper formatting
                df.to_csv(output_path, sep="\t", index=False, float_format="%.6f")

                logger.info(f"Exported {len(df)} annotations to {output_path}")
                return output_path

        except Exception as e:
            logger.error(f"Failed to export TSV: {e}")
            raise DatabaseError(f"TSV export failed: {e}")

    def get_summary_statistics(self) -> Dict[str, Any]:
        """
        Generate comprehensive summary statistics from the database.

        Returns:
            Dictionary containing detailed statistics about annotations
        """
        logger.info("Generating summary statistics")

        try:
            with self.get_connection() as conn:
                stats = {}

                # Basic counts
                cursor = conn.execute("SELECT COUNT(*) FROM domain_annotations")
                stats["total_annotations"] = cursor.fetchone()[0]

                if stats["total_annotations"] == 0:
                    logger.warning("No annotations found in database")
                    return stats

                # Annotation type breakdown
                cursor = conn.execute("""
                SELECT annotation_type, COUNT(*) 
                FROM domain_annotations 
                GROUP BY annotation_type
                """)
                for ann_type, count in cursor.fetchall():
                    stats[f"{ann_type}_annotations"] = count

                # Unique counts
                cursor = conn.execute(
                    "SELECT COUNT(DISTINCT domain_id) FROM domain_annotations"
                )
                stats["unique_domains"] = cursor.fetchone()[0]

                cursor = conn.execute(
                    "SELECT COUNT(DISTINCT go_id) FROM domain_annotations"
                )
                stats["unique_go_terms"] = cursor.fetchone()[0]

                cursor = conn.execute(
                    "SELECT COUNT(DISTINCT direct_source_term) FROM domain_annotations"
                )
                stats["unique_direct_terms"] = cursor.fetchone()[0]

                # Score statistics
                cursor = conn.execute("""
                SELECT 
                    MIN(association_score) as min_score,
                    MAX(association_score) as max_score,
                    AVG(association_score) as avg_score,
                    MIN(fdr_q_value) as min_fdr,
                    MAX(fdr_q_value) as max_fdr,
                    AVG(fdr_q_value) as avg_fdr
                FROM domain_annotations
                """)
                score_stats = cursor.fetchone()
                stats.update(
                    {
                        "min_association_score": score_stats[0],
                        "max_association_score": score_stats[1],
                        "avg_association_score": score_stats[2],
                        "min_fdr_q_value": score_stats[3],
                        "max_fdr_q_value": score_stats[4],
                        "avg_fdr_q_value": score_stats[5],
                    }
                )

                # Top scoring direct associations
                cursor = conn.execute("""
                SELECT domain_id, go_id, association_score, fdr_q_value
                FROM domain_annotations 
                WHERE annotation_type = 'direct'
                ORDER BY association_score DESC 
                LIMIT 10
                """)
                stats["top_direct_associations"] = [
                    {
                        "domain": row[0],
                        "go_term": row[1],
                        "score": row[2],
                        "fdr": row[3],
                    }
                    for row in cursor.fetchall()
                ]

                # Domain statistics
                cursor = conn.execute("""
                SELECT domain_id, COUNT(*) as annotation_count
                FROM domain_annotations
                GROUP BY domain_id
                ORDER BY annotation_count DESC
                LIMIT 10
                """)
                stats["most_annotated_domains"] = [
                    {"domain": row[0], "count": row[1]} for row in cursor.fetchall()
                ]

                # GO term statistics
                cursor = conn.execute("""
                SELECT go_id, COUNT(*) as domain_count
                FROM domain_annotations
                GROUP BY go_id
                ORDER BY domain_count DESC
                LIMIT 10
                """)
                stats["most_common_go_terms"] = [
                    {"go_term": row[0], "domain_count": row[1]}
                    for row in cursor.fetchall()
                ]

                # Metadata if available
                cursor = conn.execute(
                    "SELECT key, value, data_type FROM pipeline_metadata"
                )
                metadata = {}
                for key, value, data_type in cursor.fetchall():
                    if data_type == "numeric":
                        try:
                            metadata[key] = float(value)
                        except ValueError:
                            metadata[key] = value
                    elif data_type == "boolean":
                        metadata[key] = value.lower() == "true"
                    elif data_type == "json":
                        try:
                            metadata[key] = json.loads(value)
                        except json.JSONDecodeError:
                            metadata[key] = value
                    else:
                        metadata[key] = value

                if metadata:
                    stats["pipeline_metadata"] = metadata

                logger.info("Summary statistics generated successfully")
                return stats

        except Exception as e:
            logger.error(f"Failed to generate statistics: {e}")
            raise DatabaseError(f"Statistics generation failed: {e}")

    def query_annotations(
        self,
        domain_id: Optional[str] = None,
        go_id: Optional[str] = None,
        annotation_type: Optional[str] = None,
        limit: Optional[int] = None,
    ) -> pd.DataFrame:
        """
        Query annotations with flexible filtering.

        Args:
            domain_id: Specific domain ID to query
            go_id: Specific GO term ID to query
            annotation_type: Filter by annotation type
            limit: Maximum number of results to return

        Returns:
            DataFrame containing matching annotations
        """
        try:
            with self.get_connection() as conn:
                query_parts = ["SELECT * FROM domain_annotations WHERE 1=1"]
                params = []

                if domain_id:
                    query_parts.append("AND domain_id = ?")
                    params.append(domain_id)

                if go_id:
                    query_parts.append("AND go_id = ?")
                    params.append(go_id)

                if annotation_type:
                    query_parts.append("AND annotation_type = ?")
                    params.append(annotation_type)

                query_parts.append("ORDER BY association_score DESC")

                if limit:
                    query_parts.append(f"LIMIT {limit}")

                query = " ".join(query_parts)
                return pd.read_sql_query(query, conn, params=params)

        except Exception as e:
            logger.error(f"Query failed: {e}")
            raise DatabaseError(f"Annotation query failed: {e}")

    def get_database_info(self) -> Dict[str, Any]:
        """
        Get comprehensive database information including schema and performance metrics.

        Returns:
            Dictionary with database information
        """
        try:
            with self.get_connection() as conn:
                info = {}

                # Database file size
                info["database_path"] = str(self.db_path)
                info["database_size_mb"] = self.db_path.stat().st_size / (1024 * 1024)

                # Table information
                cursor = conn.execute("""
                SELECT name, sql FROM sqlite_master 
                WHERE type='table' 
                ORDER BY name
                """)

                tables = {}
                for name, sql in cursor.fetchall():
                    cursor.execute(f"SELECT COUNT(*) FROM {name}")
                    row_count = cursor.fetchone()[0]
                    tables[name] = {"row_count": row_count, "schema": sql}

                info["tables"] = tables

                # Index information
                cursor = conn.execute("""
                SELECT name, sql FROM sqlite_master 
                WHERE type='index' AND sql IS NOT NULL
                ORDER BY name
                """)

                info["indices"] = {name: sql for name, sql in cursor.fetchall()}

                # Database settings
                settings = {}
                for pragma in [
                    "journal_mode",
                    "synchronous",
                    "temp_store",
                    "mmap_size",
                ]:
                    cursor = conn.execute(f"PRAGMA {pragma}")
                    settings[pragma] = cursor.fetchone()[0]

                info["settings"] = settings

                return info

        except Exception as e:
            logger.error(f"Failed to get database info: {e}")
            raise DatabaseError(f"Database info retrieval failed: {e}")

    def backup_database(self, backup_path: Union[str, Path]) -> Path:
        """
        Create a backup of the database.

        Args:
            backup_path: Path for the backup file

        Returns:
            Path to the backup file
        """
        backup_path = Path(backup_path)
        backup_path.parent.mkdir(parents=True, exist_ok=True)

        logger.info(f"Creating database backup at {backup_path}")

        try:
            with self.get_connection() as conn:
                backup = sqlite3.connect(backup_path)
                conn.backup(backup)
                backup.close()

            logger.info("Database backup completed successfully")
            return backup_path

        except Exception as e:
            logger.error(f"Backup failed: {e}")
            raise DatabaseError(f"Database backup failed: {e}")

    def close(self) -> None:
        """Close database manager and clean up resources."""
        logger.info("Closing database manager")
        # SQLite connections are closed in context manager, no persistent connections to close

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()


# Example usage and testing functions
if __name__ == "__main__":
    # Example usage
    from dataclasses import dataclass
    from loguru import logger

    @dataclass
    class MockAnnotation:
        domain: str
        go_term: str
        q_value: float
        association_score: float
        annotation_type: str
        direct_source_term: str

    # Configure logger
    logger.add("dcgo_database.log", rotation="10 MB", level="INFO")

    # Test database manager
    db_path = "/tmp/test_dcgo.db"

    with dcGODatabaseManager(db_path) as db:
        # Create test annotations
        test_annotations = [
            MockAnnotation(
                "PF00001", "GO:0003677", 0.001, 85.5, "direct", "GO:0003677"
            ),
            MockAnnotation(
                "PF00001", "GO:0003674", 0.001, 85.5, "propagated", "GO:0003677"
            ),
            MockAnnotation(
                "PF00002", "GO:0008150", 0.005, 72.3, "direct", "GO:0008150"
            ),
        ]

        # Test metadata
        test_metadata = {
            "total_proteins": 50000,
            "total_domains": 15000,
            "fdr_threshold": 0.01,
            "pipeline_version": "1.0.0",
        }

        # Store data
        db.store_annotations(test_annotations, test_metadata)

        # Export TSV
        tsv_path = db.export_tsv("/tmp/test_export.tsv")
        print(f"Exported to: {tsv_path}")

        # Get statistics
        stats = db.get_summary_statistics()
        print("Summary Statistics:")
        for key, value in stats.items():
            if key not in ["top_direct_associations", "most_annotated_domains"]:
                print(f"  {key}: {value}")
