import os
from pathlib import Path

import psycopg2
import pytest

from conftest import load_app_module

TEST_DB_REQUIRED = all(
    name in os.environ
    for name in ["DB_HOST", "DB_PORT", "DB_NAME", "DB_USER", "DB_PASSWORD"]
)


def _exec_sql_file(cur, path):
    cur.execute(Path(path).read_text())


@pytest.fixture(scope="session")
def initialized_test_database():
    if not TEST_DB_REQUIRED:
        pytest.skip("Database environment variables are not configured for integration tests")

    conn = psycopg2.connect(
        host=os.environ["DB_HOST"],
        port=os.environ["DB_PORT"],
        dbname=os.environ["DB_NAME"],
        user=os.environ["DB_USER"],
        password=os.environ["DB_PASSWORD"],
    )
    conn.autocommit = True
    try:
        with conn.cursor() as cur:
            _exec_sql_file(cur, Path(__file__).with_name("data") / "test_schema.sql")
            _exec_sql_file(cur, Path(__file__).with_name("data") / "test_seed.sql")
    finally:
        conn.close()


@pytest.fixture
def db_client(monkeypatch, initialized_test_database):
    monkeypatch.setenv("DB_HOST", os.environ["DB_HOST"])
    monkeypatch.setenv("DB_PORT", os.environ["DB_PORT"])
    monkeypatch.setenv("DB_NAME", os.environ["DB_NAME"])
    monkeypatch.setenv("DB_USER", os.environ["DB_USER"])
    monkeypatch.setenv("DB_PASSWORD", os.environ["DB_PASSWORD"])
    module = load_app_module()
    module.app.config.update(TESTING=True)
    with module.app.test_client() as client:
        yield client


def test_search_loci_returns_result(db_client):
    response = db_client.post("/search_loci", data={"species": "human", "loci": "chr1:100-220"})
    assert response.status_code == 200
    assert b"chr1:100-200" in response.data
    assert b"101" in response.data


def test_search_gene_uses_configured_gene_table(db_client):
    response = db_client.post("/search_gene", data={"species": "human", "gene": "MYC", "window": "5000"})
    assert response.status_code == 200
    assert b"MYC" in response.data
    assert b"101" in response.data


def test_search_gsm_uses_summary_table(db_client):
    response = db_client.post("/search_gsm", data={"species": "human", "gsm": "GSM100001"})
    assert response.status_code == 200
    assert b"GSM100001" in response.data
    assert b"Bcell" in response.data


def test_search_celltype_returns_gain_and_metadata(db_client):
    response = db_client.post("/search_celltype", data={"species": "human", "celltype": "Bcell"})
    assert response.status_code == 200
    assert b"Bcell" in response.data
    assert b"101" in response.data
