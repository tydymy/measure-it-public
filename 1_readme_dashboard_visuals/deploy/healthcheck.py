"""Container health check: GET $MEASURE_IT_HEALTH_URL (default the API's /health); exit 0 only on HTTP 200.

The API's /health also reports `status: degraded` (HTTP 200) when data/processed or SOURCE_REGISTRY.yaml is missing:
that counts as unhealthy here, because a container without its data bundle cannot answer anything."""
import json
import os
import sys
import urllib.request

url = os.environ.get("MEASURE_IT_HEALTH_URL", "http://127.0.0.1:8000/health")
try:
    with urllib.request.urlopen(url, timeout=8) as r:
        body = r.read()
        if r.status != 200:
            sys.exit(1)
except Exception as e:  # noqa: BLE001
    print(f"unhealthy: {url}: {e}", file=sys.stderr)
    sys.exit(1)
try:
    doc = json.loads(body)
except ValueError:
    sys.exit(0)   # streamlit's /_stcore/health answers "ok" as plain text
if isinstance(doc, dict) and doc.get("status") not in (None, "ok"):
    print(f"unhealthy: {url}: status {doc.get('status')!r}: {doc.get('reason')}", file=sys.stderr)
    sys.exit(1)
sys.exit(0)
