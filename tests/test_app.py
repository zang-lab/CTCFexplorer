def test_homepage_renders(app_client):
    response = app_client.get("/")
    assert response.status_code == 200
    assert b"CTCF" in response.data


def test_invalid_loci_returns_helpful_message(app_client):
    response = app_client.post("/search_loci", data={"species": "human", "loci": "bad-input"})
    assert response.status_code == 200
    assert b"Use format chr#:start-end" in response.data


def test_invalid_union_rejected_without_db_query(app_client):
    response = app_client.post("/search_union", data={"species": "mouse", "union": "abc"})
    assert response.status_code == 200
    assert b"Union ID must be numeric" in response.data
