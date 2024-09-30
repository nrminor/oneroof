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

# Compress HTML files rendered from both the main and the developer documentation.
compress_html:
    @gzip -f docs/index.html
    @gzip -f docs/developer.html

# Render the main docs and the developer docs.
qmd: render render-dev

# Run all quarto recipes in sequence.
docs: render make-readme render-dev compress_html

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

# Freeze the project's PyPI dependencies into a standard pip `requirements.txt` file.
py-freeze:
    uv pip freeze > requirements.txt

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
python: py-lints py-format py-sort-imports # py-freeze

# Run all recipes in sequence with one another.
all: docs setup-env python docker
