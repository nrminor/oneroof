@default:
    just --list

# Render only docs/index.qmd to generate README
[group('docs')]
[private]
render-readme:
    @# Temporarily move _quarto.yml to avoid website mode
    @mv _quarto.yml _quarto.yml.bak
    quarto render docs/index.qmd --to gfm --output-dir . --no-execute-daemon
    @mv _quarto.yml.bak _quarto.yml

# Convert the resulting markdown file from `docs/index.qmd` to `README.md` in the project root and fix badge placement.
[group('docs')]
make-readme: render-readme
    #!/usr/bin/env python3
    import re
    import shutil

    # First move the file from root (where it was rendered) to README.md
    shutil.move('docs/index.md', 'README.md')

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

# Render all documentation (website, PDFs, and markdown)
[group('docs')]
docs: render-site make-readme

# Render the Quarto documentation website (HTML and markdown)
[group('docs')]
render-site:
    quarto render
    @echo "✓ Website rendered to _site/"

# Preview the Quarto documentation website locally
[group('docs')]
preview-site:
    quarto preview

# Publish documentation to GitHub Pages (requires gh-pages branch)
[group('docs')]
publish-docs:
    PRE_COMMIT_ALLOW_NO_CONFIG=1 quarto publish gh-pages --no-prompt

# Clean Quarto build artifacts
[group('docs')]
clean-site:
    rm -rf _site .quarto

# Copy docs/index.html to site root (internal recipe)
[group('docs')]
[private]
copy-index:
    @if [ -f "_site/docs/index.html" ]; then \
        cp _site/docs/index.html _site/index.html && echo "✓ Set docs/index.html as homepage"; \
    else \
        echo "⚠️  No _site directory found - using default homepage"; \
    fi

# Fix paths in the copied index.html and copy back PDFs/markdown (internal recipe)
[group('docs')]
[private]
fix-index-paths:
    @python docs/fix-index-paths.py

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

# === Rust Development Commands ===

# Check Rust code compiles without building
[group('rust')]
rs-check:
    cargo check

# Build Rust scripts in debug mode
[group('rust')]
rs-build:
    cargo build

# Build Rust scripts in release mode (with LTO)
[group('rust')]
rs-build-release:
    cargo build --release

# Run all Rust tests
[group('rust')]
rs-test:
    cargo test

# Run Rust tests with output shown
[group('rust')]
rs-test-verbose:
    cargo test -- --nocapture

# Run Clippy lints (strict, matching CI)
[group('rust')]
rs-lint:
    cargo clippy --all-targets

# Format Rust code
[group('rust')]
rs-format:
    cargo fmt

# Check Rust formatting without modifying
[group('rust')]
rs-format-check:
    cargo fmt -- --check

# Run all Rust checks (format, lint, test)
[group('rust')]
rust: rs-format rs-lint rs-test
    @echo "✅ All Rust checks passed!"

# Clean Rust build artifacts
[group('rust')]
rs-clean:
    cargo clean

# D.A.N.C.E.
[group('fun')]
ice:
    @echo "🕺 1, 2, 3, 4, fight! 🕺"
    @open "https://www.youtube.com/watch?v=sy1dYFGkPUE" || xdg-open "https://www.youtube.com/watch?v=sy1dYFGkPUE" || echo "Please open: https://www.youtube.com/watch?v=sy1dYFGkPUE"

# === Testing Commands ===

LOGDIR := 'test_logs'

[group('testing')]
test-all: test-illumina test-nanopore
    @echo "🎉 All tests (Illumina + Nanopore) completed successfully!"

alias test := test-all
alias t := test-all

# Run all illumina tests
[group('testing')]
test-illumina: clean-logs test-illumina-with-primers test-illumina-without-primers test-illumina-with-phylo test-illumina-missing-fastq test-illumina-dedup
    @echo "✅ All Illumina tests completed successfully!"

# Run all nanopore tests
[group('testing')]
test-nanopore: clean-logs test-nanopore-with-primers test-nanopore-without-primers test-nanopore-with-haplo test-nanopore-with-phylo test-nanopore-missing-input
    @echo "✅ All Nanopore tests completed successfully!"

# Individual Illumina test targets
[group('testing')]
test-illumina-with-primers:
    @echo "🧪 Running Illumina test with primers..."
    @mkdir -p {{ LOGDIR }}
    @if nextflow run . -profile illumina_test_with_primers > {{ LOGDIR }}/illumina_with_primers.log 2>&1; then \
        echo "  ✅ illumina_test_with_primers PASSED"; \
    else \
        echo "  ❌ illumina_test_with_primers FAILED"; \
        echo "     Check {{ LOGDIR }}/illumina_with_primers.log for details"; \
        exit 1; \
    fi

[group('testing')]
test-illumina-without-primers:
    @echo "🧪 Running Illumina test without primers..."
    @mkdir -p {{ LOGDIR }}
    @if nextflow run . -profile illumina_test_without_primers > {{ LOGDIR }}/illumina_without_primers.log 2>&1; then \
        echo "  ✅ illumina_test_without_primers PASSED"; \
    else \
        echo "  ❌ illumina_test_without_primers FAILED"; \
        echo "     Check {{ LOGDIR }}/illumina_without_primers.log for details"; \
        exit 1; \
    fi

# [group('testing')]
# test-illumina-with-metagenomics:
#     @echo "🧪 Running Illumina test with metagenomics..."
#     @mkdir -p {{ LOGDIR }}
#     @if nextflow run . -profile illumina_test_with_metagenomics > {{ LOGDIR }}/illumina_with_metagenomics.log 2>&1; then \
#         echo "  ✅ illumina_test_with_metagenomics PASSED"; \
#     else \
#         echo "  ❌ illumina_test_with_metagenomics FAILED"; \
#         echo "     Check {{ LOGDIR }}/illumina_with_metagenomics.log for details"; \
#         exit 1; \
#     fi

[group('testing')]
test-illumina-with-phylo:
    @echo "🧪 Running Illumina test with phylo..."
    @mkdir -p {{ LOGDIR }}
    @if nextflow run . -profile illumina_test_with_phylo > {{ LOGDIR }}/illumina_with_phylo.log 2>&1; then \
        echo "  ✅ illumina_test_with_phylo PASSED"; \
    else \
        echo "  ❌ illumina_test_with_phylo FAILED"; \
        echo "     Check {{ LOGDIR }}/illumina_with_phylo.log for details"; \
        exit 1; \
    fi

[group('testing')]
test-illumina-missing-fastq:
    @echo "🧪 Running pipeline expecting failure due to bad input..."
    @mkdir -p {{ LOGDIR }}
    @if ! nextflow run . -profile illumina_test_missing_fastq > {{ LOGDIR }}/illumina_too_small_input.log 2>&1; then \
        echo "  ✅ illumina_test_missing_fastq FAILED as expected"; \
    else \
        echo "  ❌ illumina_test_missing_fastq unexpectedly SUCCEEDED"; \
        echo "     Check {{ LOGDIR }}/illumina_test_missing_fastq.log for details"; \
        exit 1; \
    fi

[group('testing')]
test-illumina-bad-primers:
    @echo "🧪 Running pipeline expecting failure due to bad input..."
    @mkdir -p {{ LOGDIR }}
    @if ! nextflow run . -profile illumina_test_bad_primers > {{ LOGDIR }}/illumina_test_bad_primers.log 2>&1; then \
        echo "  ✅ illumina_test_bad_primers FAILED as expected"; \
    else \
        echo "  ❌ illumina_test_bad_primers unexpectedly SUCCEEDED"; \
        echo "     Check {{ LOGDIR }}/illumina_test_bad_primers.log for details"; \
        exit 1; \
    fi

[group('testing')]
test-illumina-dedup:
    @echo "🧪 Running Illumina test with dedup..."
    @mkdir -p {{ LOGDIR }}
    @if nextflow run . -profile illumina_test_dedup > {{ LOGDIR }}/illumina_test_dedup.log 2>&1; then \
        echo "  ✅ illumina_test_dedup PASSED"; \
    else \
        echo "  ❌ illumina_test_dedup FAILED"; \
        echo "     Check {{ LOGDIR }}/illumina_test_dedup.log for details"; \
        exit 1; \
    fi

# Individual Nanopore test targets
[group('testing')]
test-nanopore-with-primers:
    @echo "🧪 Running Nanopore test with primers..."
    @mkdir -p {{ LOGDIR }}
    @if nextflow run . -profile nanopore_test_with_primers > {{ LOGDIR }}/nanopore_with_primers.log 2>&1; then \
        echo "  ✅ nanopore_test_with_primers PASSED"; \
    else \
        echo "  ❌ nanopore_test_with_primers FAILED"; \
        echo "     Check {{ LOGDIR }}/nanopore_with_primers.log for details"; \
        exit 1; \
    fi

[group('testing')]
test-nanopore-without-primers:
    @echo "🧪 Running Nanopore test without primers..."
    @mkdir -p {{ LOGDIR }}
    @if nextflow run . -profile nanopore_test_without_primers > {{ LOGDIR }}/nanopore_without_primers.log 2>&1; then \
        echo "  ✅ nanopore_test_without_primers PASSED"; \
    else \
        echo "  ❌ nanopore_test_without_primers FAILED"; \
        echo "     Check {{ LOGDIR }}/nanopore_without_primers.log for details"; \
        exit 1; \
    fi

[group('testing')]
test-nanopore-with-haplo:
    @echo "🧪 Running Nanopore test with haplotyping..."
    @mkdir -p {{ LOGDIR }}
    @if nextflow run . -profile nanopore_test_with_haplo > {{ LOGDIR }}/nanopore_with_haplo.log 2>&1; then \
        echo "  ✅ nanopore_test_with_haplo PASSED"; \
    else \
        echo "  ❌ nanopore_test_with_haplo FAILED"; \
        echo "     Check {{ LOGDIR }}/nanopore_with_haplo.log for details"; \
        exit 1; \
    fi

# [group('testing')]
# test-nanopore-with-metagenomics:
#     @echo "🧪 Running Nanopore test with metagenomics..."
#     @mkdir -p {{ LOGDIR }}
#     @if nextflow run . -profile nanopore_test_with_metagenomics > {{ LOGDIR }}/nanopore_with_metagenomics.log 2>&1; then \
#         echo "  ✅ nanopore_test_with_metagenomics PASSED"; \
#     else \
#         echo "  ❌ nanopore_test_with_metagenomics FAILED"; \
#         echo "     Check {{ LOGDIR }}/nanopore_with_metagenomics.log for details"; \
#         exit 1; \
#     fi

[group('testing')]
test-nanopore-with-phylo:
    @echo "🧪 Running Nanopore test with phylo..."
    @mkdir -p {{ LOGDIR }}
    @if nextflow run . -profile nanopore_test_with_phylo > {{ LOGDIR }}/nanopore_with_phylo.log 2>&1; then \
        echo "  ✅ nanopore_test_with_phylo PASSED"; \
    else \
        echo "  ❌ nanopore_test_with_phylo FAILED"; \
        echo "     Check {{ LOGDIR }}/nanopore_with_phylo.log for details"; \
        exit 1; \
    fi

[group('testing')]
test-nanopore-missing-input:
    @echo "🧪 Running pipeline expecting failure due to bad input..."
    @mkdir -p {{ LOGDIR }}
    @if ! nextflow run . -profile nanopore_test_missing_input > {{ LOGDIR }}/illumina_too_small_input.log 2>&1; then \
        echo "  ✅ nanopore_test_missing_input FAILED as expected"; \
    else \
        echo "  ❌ nanopore_test_missing_input unexpectedly SUCCEEDED"; \
        echo "     Check {{ LOGDIR }}/nanopore_test_missing_input.log for details"; \
        exit 1; \
    fi

[group('testing')]
test-nanopore-bad-primers:
    @echo "🧪 Running pipeline expecting failure due to bad input..."
    @mkdir -p {{ LOGDIR }}
    @if ! nextflow run . -profile nanopore_test_bad_primers > {{ LOGDIR }}/nanopore_test_bad_primers.log 2>&1; then \
        echo "  ✅ nanopore_test_bad_primers FAILED as expected"; \
    else \
        echo "  ❌ nanopore_test_bad_primers unexpectedly SUCCEEDED"; \
        echo "     Check {{ LOGDIR }}/nanopore_test_bad_primers.log for details"; \
        exit 1; \
    fi

# Clean up test logs
[group('testing')]
clean-logs:
    @echo "🧹 Cleaning previous test logs..."
    @if [ -d "{{ LOGDIR }}" ]; then rm -rf {{ LOGDIR }}; fi

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
alias render := render-readme
alias doc := docs
alias env := setup-env
alias install := setup-env
alias local := setup-env
alias dev := setup-env
alias all_docker := docker
alias container_prep := docker
alias py := python

# Short aliases for common commands
# alias t := test
# alias tv := test-verbose
# alias tm := test-modules
# alias tw := test-workflows
# alias tp := test-pipeline

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
# alias tests := test

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
alias clean := clean-logs
alias cleanpy := py-test-clean

# alias tc := test-clean

alias pc := py-test-clean

# Ultra-short aliases for power users

alias r := render-readme
alias e := setup-env
alias f := py-format
alias l := py-lints
alias s := py-sort-imports
alias c := py-test-cov
alias pt := py-test
alias ptp := py-test-parallel

# alias tu := test-update
# alias tf := test-file
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

# alias update := test-update
# alias snapshot := test-update
# alias snapshots := test-update
# More common typos

alias dcos := docs
alias dosc := docs
alias dcso := docs

# alias tetst := test
# alias tets := test
# alias tes := test
# alias tst := test
# alias nftest := test
# alias nf-test := test

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

alias site := render-site
alias preview := preview-site
alias publish := publish-docs
alias docs-clean := clean-site
alias rs := render-site
alias ps := preview-site
alias pd := publish-docs
alias cs := clean-site

# Testing shortcuts
# alias test-all := test
# alias testall := test
# alias ta := test
# alias test-v := test-verbose
# alias test-m := test-modules
# alias test-w := test-workflows
# alias test-p := test-pipeline
# alias test-u := test-update
# Python testing shortcuts
# alias py-t := py-test
# alias py-c := py-test-cov
# alias py-p := py-test-parallel
# alias py-m := py-test-module
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

# Rust shortcuts

alias cargo := rs-build
alias clippy := rs-lint
alias rustfmt := rs-format
alias rsc := rs-check
alias rsb := rs-build
alias rsbr := rs-build-release
alias rst := rs-test
alias rsl := rs-lint
alias rsf := rs-format
alias cargotest := rs-test
alias cargobuild := rs-build
