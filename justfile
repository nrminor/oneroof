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

render:
    quarto render docs/index.qmd

make-readme:
    @mv docs/index.md ./README.md

render-dev:
    quarto render docs/developer.qmd

compress_html:
    @gzip -f docs/index.html
    @gzip -f docs/developer.html

qmd: render render-dev

docs: render make-readme render-dev compress_html

setup-env:
    pixi install

docker-build:
    docker build -t nrminor/dorado-and-friends:v0.1.0 .

docker-push:
    docker push nrminor/dorado-and-friends:v0.1.0

docker: docker-build docker-push

py-freeze:
    uv pip freeze > requirements.txt

py-lints:
    ruff check . --exit-zero --fix --unsafe-fixes

py-format:
    ruff format .

py-sort-imports:
    ruff check . -n --select=I --fix

python: py-lints py-format py-sort-imports # py-freeze

all: docs setup-env python docker
