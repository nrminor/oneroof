@default:
    just --list

alias readme := make-readme
alias zip := compress_html
alias doc := docs
alias env := setup-env
alias install := setup-env
alias local := setup-env
alias dev := setup-env
alias all_docker := docker
alias container_prep := docker
alias py := python
alias doit := all

# Render the main quarto document `docs/index.qmd`
render:
    quarto render docs/index.qmd

# Convert the resulting markdown file from `docs/index.qmd` to `README.md` in the project root.
make-readme:
    @mv docs/index.md ./README.md

# Render the developer documentation in docs/developer.qmd.
render-dev:
    quarto render docs/developer.qmd

# Render the developer documentation in docs/developer.qmd.
render-data-mgmt:
    quarto render docs/data_management.qmd

# Render the pipeline architecture documentation.
render-architecture:
    quarto render docs/pipeline_architecture.qmd

# Render the whats-that-file documentation.
render-whats-that-file:
    quarto render docs/whats-that-file.qmd

# Compress HTML files rendered from both the main and the developer documentation.
compress_html:
    @gzip -f docs/index.html
    @gzip -f docs/developer.html
    @gzip -f docs/pipeline_architecture.html
    @gzip -f docs/whats-that-file.html
    # @gzip -f docs/_management.html

# Render the main docs and the developer docs.
qmd: render render-dev render-data-mgmt render-architecture render-whats-that-file

# Run all quarto recipes in sequence.
docs: render make-readme render-dev render-architecture render-whats-that-file compress_html

# Set up the Python environment with Pixi.
setup-env:
    pixi install --frozen

# Build the Docker image for the pipeline locally.
docker-build:
    docker build -t nrminor/dorado-and-friends:v0.2.3 .

# Push the docker image to Docker Hub (requires N.R. Minor's login credentials).
docker-push:
    docker push nrminor/dorado-and-friends:v0.2.3

# Run both docker recipes in sequence
docker: docker-build docker-push

# Export the project's PyPI dependencies to requirements format.
py-freeze:
    uv export --format requirements-txt > requirements.txt

# Run all Ruff Python lints
py-lints:
    ruff check . --exit-zero --fix --unsafe-fixes

# Format all Python scripts in the current working directory.
py-format:
    ruff format .

# Sort Python imports in all Python scripts in the project.
py-sort-imports:
    ruff check . -n --select=I --fix

# Run all Python recipes in sequence.
python: py-lints py-format py-sort-imports

# === Python Testing Commands ===

# Run all Python tests with pytest
py-test:
    pytest bin/ -v

# Run Python tests with coverage report
py-test-cov:
    pytest bin/ --cov=bin --cov-report=html --cov-report=term-missing

# Run Python tests for a specific module
py-test-module MODULE:
    pytest bin/test_{{MODULE}}.py -v

# Run Python tests across all supported Python versions
py-test-tox:
    tox

# Run Python tests in parallel
py-test-parallel:
    pytest bin/ -n auto -v

# Clean Python test artifacts
py-test-clean:
    rm -rf .coverage
    rm -rf htmlcov/
    rm -rf .pytest_cache/
    rm -rf .tox/
    find bin -name "__pycache__" -type d -exec rm -rf {} +

# Run all recipes in sequence with one another.
all: docs setup-env python docker

# === Testing Commands ===

# Run all nf-test tests
test:
    nf-test test --profile test,docker

# Run tests with verbose output
test-verbose:
    nf-test test --verbose --profile test,docker

# Run only module tests
test-modules:
    nf-test test --tag modules --profile test,docker

# Run only workflow tests
test-workflows:
    nf-test test --tag workflows --profile test,docker

# Run only pipeline tests
test-pipeline:
    nf-test test --tag pipeline --profile test,docker

# Run a specific test file
test-file FILE:
    nf-test test {{FILE}} --profile test,docker

# Update test snapshots after intentional changes
test-update:
    nf-test test --update-snapshot --profile test,docker

# Clean test outputs and working directories
test-clean:
    rm -rf .nf-test/
    rm -rf tests/output/
    rm -rf tests/work/
    rm -rf tests/.nextflow/

# === Globus Integration Commands ===

# Set up Globus environment (copies template to .env if not exists)
globus-setup:
    @if [ ! -f globus/config/.env ]; then \
        echo "Creating Globus configuration from template..."; \
        cp globus/config/.env.template globus/config/.env; \
        echo "Please edit globus/config/.env with your values"; \
    else \
        echo "Globus config already exists at globus/config/.env"; \
    fi

# Deploy the Globus action provider
globus-deploy: globus-setup
    cd globus/scripts && ./deploy.sh

# Deploy action provider as systemd service (requires sudo)
globus-deploy-systemd: globus-setup
    cd globus/scripts && sudo ./deploy.sh --systemd

# Deploy action provider with Docker
globus-deploy-docker: globus-setup
    cd globus/scripts && ./deploy.sh --docker

# Register the OneRoof flow with Globus
globus-register: globus-setup
    cd globus/scripts && ./register_flow.py

# Test the Globus flow with sample Nanopore data
globus-test-nanopore flow_id="" input="/test/pod5" primers="/test/primers.bed" ref="/test/ref.fasta" gbk="/test/ref.gbk" output="/test/output":
    #!/usr/bin/env bash
    cd globus/scripts
    if [ -n "{{flow_id}}" ]; then
        ./test_flow.py --flow-id {{flow_id}} --platform nanopore --input-data {{input}} --primer-bed {{primers}} --refseq {{ref}} --ref-gbk {{gbk}} --output-path {{output}}
    else
        ./test_flow.py --platform nanopore --input-data {{input}} --primer-bed {{primers}} --refseq {{ref}} --ref-gbk {{gbk}} --output-path {{output}}
    fi

# Test the Globus flow with sample Illumina data
globus-test-illumina flow_id="" input="/test/fastq" primers="/test/primers.bed" ref="/test/ref.fasta" gbk="/test/ref.gbk" output="/test/output":
    #!/usr/bin/env bash
    cd globus/scripts
    if [ -n "{{flow_id}}" ]; then
        ./test_flow.py --flow-id {{flow_id}} --platform illumina --input-data {{input}} --primer-bed {{primers}} --refseq {{ref}} --ref-gbk {{gbk}} --output-path {{output}}
    else
        ./test_flow.py --platform illumina --input-data {{input}} --primer-bed {{primers}} --refseq {{ref}} --ref-gbk {{gbk}} --output-path {{output}}
    fi

# Check status of Globus action provider
globus-status:
    @if systemctl is-active --quiet oneroof-action-provider; then \
        echo "Action provider is running as systemd service"; \
        systemctl status oneroof-action-provider; \
    elif docker ps | grep -q oneroof-action-provider; then \
        echo "Action provider is running in Docker"; \
        docker logs --tail 20 oneroof-action-provider; \
    else \
        echo "Action provider is not running"; \
    fi

# View Globus action provider logs
globus-logs:
    @if systemctl is-active --quiet oneroof-action-provider; then \
        journalctl -u oneroof-action-provider -f; \
    elif docker ps | grep -q oneroof-action-provider; then \
        docker logs -f oneroof-action-provider; \
    else \
        echo "Action provider is not running"; \
    fi

# Stop Globus action provider
globus-stop:
    @if systemctl is-active --quiet oneroof-action-provider; then \
        sudo systemctl stop oneroof-action-provider; \
    elif docker ps | grep -q oneroof-action-provider; then \
        docker stop oneroof-action-provider; \
    else \
        echo "Action provider is not running"; \
    fi

# Clean up Globus temporary files and logs
globus-clean:
    rm -f globus/scripts/temp_*.json
    rm -f globus/config/flow_id.txt

# Full Globus setup: configure, deploy, and register
globus-init: globus-setup globus-deploy globus-register
    @echo "Globus integration setup complete!"
    @echo "Test with: just globus-test-nanopore or just globus-test-illumina"
