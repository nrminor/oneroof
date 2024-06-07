@default:
    just --list

alias readme := make_readme
alias zip := compress_html
alias install := setup_local_env
alias local := setup_local_env
alias dev := setup_local_env
alias all_docker := docker
alias container_prep := docker

render:
    quarto render docs/base_and_variant_call.qmd

make_readme:
    @mv docs/base_and_variant_call.md ./README.md

render_dev:
    quarto render docs/developer.qmd

compress_html:
    @gzip -f docs/base_and_variant_call.html
    @gzip -f docs/developer.html

qmd: render render_dev

docs: render make_readme render_dev compress_html

setup_local_env:
    pixi install

docker_build:
    docker build -t nrminor/dorado-and-friends:v0.1.0 .

docker_push:
    docker push nrminor/dorado-and-friends:v0.1.0

docker: docker_build docker_push

all: docs docker
