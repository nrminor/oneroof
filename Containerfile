FROM ubuntu:20.04

# Set working directory
WORKDIR /scratch

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York
ENV HOME=/home

# Set default command to be the bash shell
ENTRYPOINT ["bash"]

# Install a few ubuntu dependencies
RUN apt-get update && \
    apt install --fix-broken && \
    apt-get install -y \
    build-essential \
    curl \
    wget \
    make \
    gcc \
    cmake \
    libxml2-dev \
    libxslt-dev \
    libffi-dev \
    git && \
    apt install --fix-broken && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    mkdir /dependencies && \
    dpkg -l > /dependencies/apt-get.lock

# pull the latest dorado executable
RUN cd $HOME && \
    wget --quiet https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.7.1-linux-x64.tar.gz && \
    tar -xvf dorado-0.7.1-linux-x64.tar.gz && \
    rm -rf dorado-0.7.1-linux-x64.tar.gz

# add the dorado files to $PATH
ENV PATH=$PATH:$HOME/dorado-0.7.1-linux-x64/bin:$HOME/dorado-0.7.1-linux-x64/lib:$HOME/dorado-0.7.1-linux-x64

# predownload dorado models so that dorado can basecall offline
RUN cd $HOME && dorado download

# Install everything else with Pixi:
# --------------------------------
# 1) copy the required dependency and configuration files into the image
#    Note: bin/oneroof_cli/ is needed for the editable install declared in pyproject.toml
COPY pyproject.toml $HOME/pyproject.toml
COPY pixi.lock $HOME/pixi.lock
COPY bin/oneroof_cli $HOME/bin/oneroof_cli

# 2) install pixi
RUN cd $HOME && PIXI_ARCH=x86_64 curl -fsSL https://pixi.sh/install.sh | bash

# 3) make sure pixi and pixi installs are on the $PATH
ENV PATH=$PATH:$HOME/.pixi/bin

# 4) install everything else with pixi (rust-script brings in the full Rust toolchain)
RUN cd $HOME && pixi install --verbose --color=always --frozen && pixi clean cache --assume-yes

# 5) modify the shell config so that each container launches within the pixi env
RUN echo "export PATH=$PATH:$HOME/.pixi/envs/default/bin" >> $HOME/.bashrc

# 6) modify some nextflow environment variables
RUN echo "export NXF_CACHE_DIR=/scratch" >> $HOME/.bashrc
RUN echo "export NXF_HOME=/scratch" >> $HOME/.bashrc

# Pre-compile Rust scripts for container use
# --------------------------------
# This avoids rust-script JIT compilation overhead on every process invocation.
# Local (non-container) runs still use rust-script directly via the .rs files.
# Note: rust-script depends on the full Rust toolchain, so cargo is available.
# We must add the pixi bin to PATH so cargo can find rustc.
COPY Cargo.toml $HOME/Cargo.toml
COPY Cargo.lock $HOME/Cargo.lock
COPY bin/find_and_trim_amplicons.rs $HOME/bin/find_and_trim_amplicons.rs
RUN cd $HOME && \
    export PATH="$HOME/.pixi/envs/default/bin:$PATH" && \
    RUSTFLAGS="-C target-cpu=native" cargo build --release && \
    cp $HOME/target/release/find_and_trim_amplicons $HOME/.pixi/envs/default/bin/
