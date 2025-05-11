{
  description = "Reproducible dev shell for the `oneroof` bioinformatic processing pipeline";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
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

        dorado =
          if system == "x86_64-linux" then
            pkgs.stdenv.mkDerivation {
              name = "dorado";
              src = pkgs.fetchurl {
                url = "https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.7.1-linux-x64.tar.gz";
                sha256 = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855";
              };

              unpackPhase = "tar -xvf $src";
              installPhase = ''
                set -euo pipefail
                mkdir -p $out
                cp -r dorado-0.7.1-linux-x64/* $out/
              '';
            }
          else
            null;

      in
      {
        devShells.default = pkgs.mkShell {
          name = "oneroof";

          buildInputs = [
            pkgs.stdenv
            pkgs.gcc
            pkgs.curl
            pkgs.wget
            pkgs.openjdk
            pkgs.git
            pkgs.cmake
            pkgs.libxml2
            pkgs.libxslt
            pkgs.libffi
            pkgs.pixi
            dorado
          ] ++ pkgs.lib.optional (dorado != null) dorado;

          shellHook = ''
            ${pkgs.lib.optionalString (dorado != null) ''
              export PATH=$PATH:${dorado}/bin:${dorado}/lib
            ''}
            echo "🔧 Entering oneroof dev shell"
            export PS1="(oneroof) $PS1"
            if [ ! -d .pixi/envs/default ]; then
              echo "Pixi env not found. Running install..."
              pixi install --frozen
            fi

            export PATH="$PWD/.pixi/envs/default/bin:$PATH"
          '';
        };
      }
    );
}
