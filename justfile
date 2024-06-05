@default:
    just --list

alias readme := make_readme
alias zip := compress_html
alias install := setup_local_env
alias local := setup_local_env
alias dev := setup_local_env
alias all_docker := docker
alias container_prep := docker

make_readme:
    @mv docs/base_and_variant_call_on_chtc.md ./README.md

compress_html:
    @gzip docs/base_and_variant_call_on_chtc.html

docs: make_readme compress_html

setup_local_env:
    pixi install

docker_build:
    docker build -t nrminor/dorado-and-friends:v0.1.0 .

docker_push:
    docker push nrminor/dorado-and-friends:v0.1.0

docker: docker_build docker_push
