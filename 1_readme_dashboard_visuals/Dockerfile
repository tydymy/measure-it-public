# syntax=docker/dockerfile:1.7
#
# measure-it-public: one image, three services (API, dashboard, MCP server). Read-only, offline.
#
#   docker build -t measure-it:local .                              # code only; data mounted at run time
#   docker build -t measure-it:baked --target baked \
#       --build-context bundle=dist/bundle .                        # code + an unpacked data bundle baked in
#
# The image never holds raw data. At run time it expects an unpacked data bundle (`measure-it bundle`,
# docs/DEPLOYMENT.md) at /app: data/processed, the raw MANIFEST / DATA_AUDIT sidecars, results/, configs/ and
# SOURCE_REGISTRY.yaml. MEASURE_IT_OFFLINE=1 is set, so no tool reaches the network.

ARG PYTHON_IMAGE=python:3.12-slim-bookworm
ARG UV_VERSION=0.12.2

# ------------------------------------------------------------------------------------------------ uv binary
FROM ghcr.io/astral-sh/uv:${UV_VERSION} AS uv

# ------------------------------------------------------------------------------------------------ builder
FROM ${PYTHON_IMAGE} AS builder
COPY --from=uv /uv /uvx /bin/
# a C compiler for the one dependency without an aarch64 wheel (mat-io builds from its sdist); builder stage only
RUN apt-get update \
    && apt-get install -y --no-install-recommends gcc libc6-dev \
    && rm -rf /var/lib/apt/lists/*
ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PYTHON_DOWNLOADS=never \
    UV_PYTHON=/usr/local/bin/python3.12 \
    UV_PROJECT_ENVIRONMENT=/app/.venv
WORKDIR /app
# 1. third-party dependencies only (cached layer; changes only when uv.lock changes)
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --frozen --no-dev --no-install-project
# 2. the project itself (editable, so measure_it.config.PROJECT_ROOT resolves to /app)
COPY pyproject.toml uv.lock README.md ./
COPY src ./src
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-dev

# ------------------------------------------------------------------------------------------------ runtime
FROM ${PYTHON_IMAGE} AS runtime
LABEL org.opencontainers.image.title="measure-it-public" \
      org.opencontainers.image.description="Measure It to Cure It: public-data engine (API, dashboard, MCP), read-only" \
      org.opencontainers.image.licenses="MIT"
RUN groupadd --gid 10001 app && useradd --uid 10001 --gid app --create-home --home-dir /home/app app
WORKDIR /app
COPY --from=builder /app/.venv /app/.venv
COPY pyproject.toml uv.lock README.md SOURCE_REGISTRY.yaml ./
COPY src ./src
COPY dashboard ./dashboard
COPY configs ./configs
COPY docs ./docs
COPY deploy/healthcheck.py /usr/local/bin/measure-it-healthcheck
# mount points of the data bundle (empty in the code-only image)
RUN mkdir -p data/processed data/raw data/_http_cache data/interim results/tables \
    && chown -R app:app /app/data /app/results
ENV PATH=/app/.venv/bin:$PATH \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MEASURE_IT_OFFLINE=1 \
    HOME=/home/app \
    MPLCONFIGDIR=/tmp/matplotlib \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_BROWSER_SERVER_ADDRESS=localhost \
    MEASURE_IT_HEALTH_URL=http://127.0.0.1:8000/health
USER app
EXPOSE 8000 8501 8765
# The URL is per service (compose sets MEASURE_IT_HEALTH_URL); the default is the API.
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD ["python", "/usr/local/bin/measure-it-healthcheck"]
CMD ["measure-it", "serve-api", "--host", "0.0.0.0", "--port", "8000", "--offline"]

# ------------------------------------------------------------------------------------------------ baked (optional)
# `--build-context bundle=<unpacked bundle dir>`: copies the data bundle into the image (about 2.2 GB).
FROM runtime AS baked
COPY --from=bundle --chown=app:app . /app/

# ------------------------------------------------------------------------------------------------ default target
# `docker build .` without --target builds this (= runtime, no data); the baked stage above is skipped.
FROM runtime AS default
