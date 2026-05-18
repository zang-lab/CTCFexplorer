import importlib
import os
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def load_app_module():
    if "app" in sys.modules:
        return importlib.reload(sys.modules["app"])
    import app  # noqa: F401
    return sys.modules["app"]


@pytest.fixture
def app_client(monkeypatch):
    monkeypatch.setenv("DB_NAME", os.environ.get("DB_NAME", "CTCFDB_PostgreSQL"))
    monkeypatch.setenv("DB_USER", os.environ.get("DB_USER", "postgres"))
    monkeypatch.setenv("DB_PASSWORD", os.environ.get("DB_PASSWORD", ""))
    monkeypatch.setenv("DB_HOST", os.environ.get("DB_HOST", "localhost"))
    monkeypatch.setenv("DB_PORT", os.environ.get("DB_PORT", "5432"))
    module = load_app_module()
    module.app.config.update(TESTING=True)
    with module.app.test_client() as client:
        yield client
