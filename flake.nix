{
  description = "Reproducible dev shell for the `oneroof` bioinformatic processing pipeline";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    flake-utils = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };

        dorado = pkgs.stdenv.mkDerivation {
          name = "dorado";
          src = pkgs.fetchurl {
            url = "https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.7.1-linux-x64.tar.gz";
            sha256 = "<fill-me-in>"; # run nix with fake sha256 to get it
          };

          unpackPhase = "tar -xvf $src";
          installPhase = ''
            set -euo pipefail
            mkdir -p $out
            cp -r dorado-0.7.1-linux-x64/* $out/
          '';
        };

      in
      {
        devShells.default = pkgs.mkShell {
          name = "oneroof-env";

          buildInputs = [
            pkgs.stdenv
            pkgs.gcc
            pkgs.curl
            pkgs.wget
            pkgs.make
            pkgs.git
            pkgs.cmake
            pkgs.libxml2
            pkgs.libxslt
            pkgs.libffi
            pkgs.pixi
            dorado
          ];

          shellHook = ''
            export PATH=$PATH:${dorado}/bin:${dorado}/lib
            export NXF_CACHE_DIR=/scratch
            export NXF_HOME=/scratch
            echo "Activating pixi..."
            cd ~/  # or $HOME if you define it
            pixi install --verbose --color=always --frozen
            pixi shell hook >> ~/.bashrc

          '';
        };
      }
    );
}
