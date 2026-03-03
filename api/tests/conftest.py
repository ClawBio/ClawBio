"""Shared fixtures for API tests."""

import pytest
from fastapi.testclient import TestClient

from api.app import create_app
from api.jobs import JobStore, store as _default_store


@pytest.fixture()
def app():
    return create_app()


@pytest.fixture()
def client(app):
    return TestClient(app)


@pytest.fixture(autouse=True)
def _clear_store():
    """Reset the in-memory store between tests."""
    _default_store._jobs.clear()
    yield
    _default_store._jobs.clear()
