FROM python:3.12-slim

ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app

# Package metadata + sources first so dependency installation is cached.
COPY pyproject.toml README.md ./
COPY src ./src
RUN pip install --upgrade pip \
    && pip install ".[test]" \
    && pip install setuptools wheel build

# Tests and the one-shot verifier (not needed to serve the API).
COPY tests ./tests
COPY verify ./verify

EXPOSE 8000
CMD ["uvicorn", "nanopore_aligner.main:app", "--host", "0.0.0.0", "--port", "8000"]
