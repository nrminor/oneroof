@default:
    just --list

# Render only docs/index.qmd to generate README
[private]
[group('docs')]
render-readme:
    quarto render docs/index.qmd --to gfm --output-dir .

# Convert the resulting markdown file from `docs/index.qmd` to `README.md` in the project root and fix badge placement.
[group('docs')]
make-readme: render-readme
    #!/usr/bin/env -S uv run
    import re
    import shutil

    # First move the file from root (where it was rendered) to README.md
    shutil.move('index.md', 'README.md')

    # Read the README
    with open('README.md', 'r') as f:
        content = f.read()

    # Find the badge line (starts with [![Nextflow])
    badge_pattern = r'^(\[\!\[Nextflow\].*?)$'
    match = re.search(badge_pattern, content, re.MULTILINE)

    if match:
        badge_line = match.group(1)
        # Remove badge from current position (and any trailing newlines)
        content = re.sub(badge_pattern + r'\n+', '', content, flags=re.MULTILINE)
        # Split into lines
        lines = content.split('\n')
        # Insert badge after title (line 1)
        lines.insert(1, badge_line)
        # Write back
        with open('README.md', 'w') as f:
            f.write('\n'.join(lines))

# Render the developer documentation in docs/developer.qmd.
[private]
[group('docs')]
render-dev:
    quarto render docs/developer.qmd

# Render the data management documentation in docs/data_management.qmd.
[private]
[group('docs')]
render-data-mgmt:
    quarto render docs/data_management.qmd

# Render the pipeline architecture documentation.
[private]
[group('docs')]
render-architecture:
    quarto render docs/pipeline_architecture.qmd

# Render the whats-that-file documentation.
[private]
[group('docs')]
render-whats-that-file:
    quarto render docs/whats-that-file.qmd

# Compress HTML files rendered from both the main and the developer documentation.
[private]
[group('docs')]
compress_html:
    @gzip -f docs/index.html
    @gzip -f docs/developer.html
    @gzip -f docs/pipeline_architecture.html
    @gzip -f docs/whats-that-file.html
    # @gzip -f docs/_management.html

# Render the main docs and the developer docs.
[private]
[group('docs')]
qmd: render-readme render-dev render-data-mgmt render-architecture render-whats-that-file

# Run all quarto recipes in sequence.
[group('docs')]
docs: make-readme render-dev render-architecture render-whats-that-file compress_html

# Render the Quarto documentation website
[group('docs')]
render-site:
    quarto render
    @echo "✓ Website rendered to _site/"
    @# Copy PDFs back to docs directory if desired
    @if [ -d "_site/docs" ]; then \
        cp _site/docs/*.pdf docs/ 2>/dev/null && echo "✓ PDFs copied to docs/" || true; \
        cp _site/docs/*.md docs/ 2>/dev/null && echo "✓ Markdown files copied to docs/" || true; \
    fi

# Preview the Quarto documentation website locally
[group('docs')]
preview-site:
    quarto preview

# Publish documentation to GitHub Pages (requires gh-pages branch)
[group('docs')]
publish-docs:
    quarto publish gh-pages --no-prompt

# Clean Quarto build artifacts
[group('docs')]
clean-site:
    rm -rf _site .quarto

# Copy docs/index.html to site root (internal recipe)
[private]
[group('docs')]
copy-index:
    @if [ -f "_site/docs/index.html" ]; then \
        cp _site/docs/index.html _site/index.html && echo "✓ Set docs/index.html as homepage"; \
    else \
        echo "⚠️  No _site directory found - using default homepage"; \
    fi

# Set up the Python environment with Pixi.
[group('setup')]
setup-env:
    pixi install --frozen

# Build the Docker image for the pipeline locally.
[group('docker')]
docker-build:
    docker build -t nrminor/dorado-and-friends:v0.2.3 .

# Push the docker image to Docker Hub (requires N.R. Minor's login credentials).
[group('docker')]
docker-push:
    docker push nrminor/dorado-and-friends:v0.2.3

# Run both docker recipes in sequence
[group('docker')]
docker: docker-build docker-push

# Export the project's PyPI dependencies to requirements format.
[group('python')]
py-freeze:
    uv export --format requirements-txt > requirements.txt

# Run all Ruff Python lints
[group('python')]
py-lints:
    ruff check . --exit-zero --fix --unsafe-fixes

# Format all Python scripts in the current working directory.
[group('python')]
py-format:
    ruff format .

# Sort Python imports in all Python scripts in the project.
[group('python')]
py-sort-imports:
    ruff check . -n --select=I --fix

# Run all Python recipes in sequence.
[group('python')]
python: py-lints py-format py-sort-imports

# === Python Testing Commands ===

# Run all Python tests with pytest
[group('python-test')]
py-test:
    pytest bin/ -v

# Run Python tests with coverage report
[group('python-test')]
py-test-cov:
    pytest bin/ --cov=bin --cov-report=html --cov-report=term-missing

# Run Python tests for a specific module
[group('python-test')]
py-test-module MODULE:
    pytest bin/test_{{ MODULE }}.py -v

# Run Python tests across all supported Python versions
[group('python-test')]
py-test-tox:
    tox

# Run Python tests in parallel
[group('python-test')]
py-test-parallel:
    pytest bin/ -n auto -v

# Clean Python test artifacts
[group('python-test')]
py-test-clean:
    rm -rf .coverage
    rm -rf htmlcov/
    rm -rf .pytest_cache/
    rm -rf .tox/
    find bin -name "__pycache__" -type d -exec rm -rf {} +

# D.A.N.C.E.
[group('fun')]
ice:
    @echo "🕺 1, 2, 3, 4, fight! 🕺"
    @open "https://www.youtube.com/watch?v=sy1dYFGkPUE" || xdg-open "https://www.youtube.com/watch?v=sy1dYFGkPUE" || echo "Please open: https://www.youtube.com/watch?v=sy1dYFGkPUE"

# === Testing Commands ===

# Run all nf-test tests
[group('nf-test')]
test:
    nf-test test --profile test,docker

# Run tests with verbose output
[group('nf-test')]
test-verbose:
    nf-test test --verbose --profile test,docker

# Run only module tests
[group('nf-test')]
test-modules:
    nf-test test --tag modules --profile test,docker

# Run only workflow tests
[group('nf-test')]
test-workflows:
    nf-test test --tag workflows --profile test,docker

# Run only pipeline tests
[group('nf-test')]
test-pipeline:
    nf-test test --tag pipeline --profile test,docker

# Run a specific test file
[group('nf-test')]
test-file FILE:
    nf-test test {{ FILE }} --profile test,docker

# Update test snapshots after intentional changes
[group('nf-test')]
test-update:
    nf-test test --update-snapshot --profile test,docker

# Clean test outputs and working directories
[group('nf-test')]
test-clean:
    rm -rf .nf-test/
    rm -rf tests/output/
    rm -rf tests/work/
    rm -rf tests/.nextflow/

# === Globus Integration Commands ===

# Set up Globus environment (copies template to .env if not exists)
[group('globus')]
globus-setup:
    @if [ ! -f globus/config/.env ]; then \
        echo "Creating Globus configuration from template..."; \
        cp globus/config/.env.template globus/config/.env; \
        echo "Please edit globus/config/.env with your values"; \
    else \
        echo "Globus config already exists at globus/config/.env"; \
    fi

# Deploy the Globus action provider
[group('globus')]
globus-deploy: globus-setup
    cd globus/scripts && ./deploy.sh

# Deploy action provider as systemd service (requires sudo)
[group('globus')]
globus-deploy-systemd: globus-setup
    cd globus/scripts && sudo ./deploy.sh --systemd

# Deploy action provider with Docker
[group('globus')]
globus-deploy-docker: globus-setup
    cd globus/scripts && ./deploy.sh --docker

# Register the OneRoof flow with Globus
[group('globus')]
globus-register: globus-setup
    cd globus/scripts && ./register_flow.py

# Test the Globus flow with sample Nanopore data
[group('globus')]
globus-test-nanopore flow_id="" input="/test/pod5" primers="/test/primers.bed" ref="/test/ref.fasta" gbk="/test/ref.gbk" output="/test/output":
    #!/usr/bin/env bash
    cd globus/scripts
    if [ -n "{{ flow_id }}" ]; then
        ./test_flow.py --flow-id {{ flow_id }} --platform nanopore --input-data {{ input }} --primer-bed {{ primers }} --refseq {{ ref }} --ref-gbk {{ gbk }} --output-path {{ output }}
    else
        ./test_flow.py --platform nanopore --input-data {{ input }} --primer-bed {{ primers }} --refseq {{ ref }} --ref-gbk {{ gbk }} --output-path {{ output }}
    fi

# Test the Globus flow with sample Illumina data
[group('globus')]
globus-test-illumina flow_id="" input="/test/fastq" primers="/test/primers.bed" ref="/test/ref.fasta" gbk="/test/ref.gbk" output="/test/output":
    #!/usr/bin/env bash
    cd globus/scripts
    if [ -n "{{ flow_id }}" ]; then
        ./test_flow.py --flow-id {{ flow_id }} --platform illumina --input-data {{ input }} --primer-bed {{ primers }} --refseq {{ ref }} --ref-gbk {{ gbk }} --output-path {{ output }}
    else
        ./test_flow.py --platform illumina --input-data {{ input }} --primer-bed {{ primers }} --refseq {{ ref }} --ref-gbk {{ gbk }} --output-path {{ output }}
    fi

# Check status of Globus action provider
[group('globus')]
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
[group('globus')]
globus-logs:
    @if systemctl is-active --quiet oneroof-action-provider; then \
        journalctl -u oneroof-action-provider -f; \
    elif docker ps | grep -q oneroof-action-provider; then \
        docker logs -f oneroof-action-provider; \
    else \
        echo "Action provider is not running"; \
    fi

# Stop Globus action provider
[group('globus')]
globus-stop:
    @if systemctl is-active --quiet oneroof-action-provider; then \
        sudo systemctl stop oneroof-action-provider; \
    elif docker ps | grep -q oneroof-action-provider; then \
        docker stop oneroof-action-provider; \
    else \
        echo "Action provider is not running"; \
    fi

# Clean up Globus temporary files and logs
[group('globus')]
globus-clean:
    rm -f globus/scripts/temp_*.json
    rm -f globus/config/flow_id.txt

# Full Globus setup: configure, deploy, and register
[group('globus')]
globus-init: globus-setup globus-deploy globus-register
    @echo "Globus integration setup complete!"
    @echo "Test with: just globus-test-nanopore or just globus-test-illumina"

# ALIASES! MORE ALIASES!

alias readme := make-readme
[private]
alias render := render-readme
[private]
alias zip := compress_html
alias doc := docs
alias env := setup-env
alias install := setup-env
alias local := setup-env
alias dev := setup-env
alias all_docker := docker
alias container_prep := docker
alias py := python

# Short aliases for common commands

alias t := test
alias tv := test-verbose
alias tm := test-modules
alias tw := test-workflows
alias tp := test-pipeline
alias d := docs
alias b := docker-build
alias p := docker-push
alias fmt := py-format
alias lint := py-lints
alias cov := py-test-cov

# More intuitive names

alias build := docker-build
alias push := docker-push
alias format := py-format
alias coverage := py-test-cov
alias help := default

# Common mistypings

alias tests := test
alias pytest := py-test
alias tox := py-test-tox
alias build-docker := docker-build
alias docker-built := docker-build
alias foramt := py-format
alias formta := py-format
alias pytset := py-test
alias pytes := py-test

# Workflow shortcuts

alias g := globus-status
alias gs := globus-setup
alias gi := globus-init
alias clean := test-clean
alias cleanpy := py-test-clean
alias tc := test-clean
alias pc := py-test-clean

# Ultra-short aliases for power users

[private]
alias r := render-readme
alias e := setup-env
alias f := py-format
alias l := py-lints
alias s := py-sort-imports
alias c := py-test-cov
alias pt := py-test
alias ptp := py-test-parallel
alias tu := test-update
alias tf := test-file

# Intuitive alternatives

alias setup := setup-env
alias environment := setup-env
alias venv := setup-env
alias virtualenv := setup-env
alias pyenv := setup-env
alias check := py-lints
alias fix := py-lints
alias imports := py-sort-imports
alias sort := py-sort-imports
alias parallel := py-test-parallel
alias update := test-update
alias snapshot := test-update
alias snapshots := test-update

# More common typos

alias dcos := docs
alias dosc := docs
alias dcso := docs
alias tetst := test
alias tets := test
alias tes := test
alias tst := test
alias nftest := test
alias nf-test := test
alias buil := docker-build
alias biuld := docker-build
alias puhs := docker-push
alias psuh := docker-push
alias formatt := py-format
alias fomat := py-format
alias fromat := py-format
alias lnit := py-lints
alias litn := py-lints
alias covreage := py-test-cov
alias coverge := py-test-cov

# Pipeline-specific shortcuts

alias nanopore := globus-test-nanopore
alias illumina := globus-test-illumina
alias ont := globus-test-nanopore
alias nano := globus-test-nanopore
alias ill := globus-test-illumina

# Documentation shortcuts

[private]
alias rd := render-dev
[private]
alias ra := render-architecture
[private]
alias rw := render-whats-that-file
[private]
alias rdm := render-data-mgmt
[private]
alias compress := compress_html
[private]
alias gzip := compress_html
alias site := render-site
alias preview := preview-site
alias publish := publish-docs
alias docs-clean := clean-site
alias rs := render-site
alias ps := preview-site
alias pd := publish-docs
alias cs := clean-site

# Testing shortcuts

alias test-all := test
alias testall := test
alias ta := test
alias test-v := test-verbose
alias test-m := test-modules
alias test-w := test-workflows
alias test-p := test-pipeline
alias test-u := test-update

# Python testing shortcuts

alias py-t := py-test
alias py-c := py-test-cov
alias py-p := py-test-parallel
alias py-m := py-test-module

# Container shortcuts

alias container := docker
alias containers := docker
alias image := docker-build
alias images := docker-build

# Globus shortcuts

alias gd := globus-deploy
alias gr := globus-register
alias gl := globus-logs
alias gst := globus-stop
alias gc := globus-clean

# Quick status checks

alias status := globus-status
alias logs := globus-logs
alias stop := globus-stop
