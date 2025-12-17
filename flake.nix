{
  description = "Reproducible dev shell for the `oneroof` bioinformatic processing pipeline";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    rust-overlay.url = "github:oxalica/rust-overlay";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      rust-overlay,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
          overlays = [ rust-overlay.overlays.default ];
        };

        # Rust toolchain with components needed for development
        rustToolchain = pkgs.rust-bin.stable.latest.default.override {
          extensions = [
            "rust-src" # Required for rust-analyzer go-to-definition on std
            "rust-analyzer" # LSP server
            "clippy" # Linter
            "rustfmt" # Formatter
          ];
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

        # Build refman using cargo-binstall
        refman = pkgs.stdenv.mkDerivation {
          name = "refman";
          version = "latest";

          nativeBuildInputs = with pkgs; [
            cargo-binstall
            pkg-config
            openssl
            cacert
          ];

          buildInputs = with pkgs; [
            openssl
          ];

          dontUnpack = true;
          dontBuild = true;

          installPhase = ''
            export HOME=$TMPDIR
            export CARGO_HOME=$TMPDIR/.cargo
            export RUSTUP_HOME=$TMPDIR/.rustup
            export SSL_CERT_FILE=${pkgs.cacert}/etc/ssl/certs/ca-bundle.crt

            # Create output directory
            mkdir -p $out/bin

            # Install refman using cargo-binstall
            ${pkgs.cargo-binstall}/bin/cargo-binstall \
              --no-confirm \
              --no-symlinks \
              --root $TMPDIR \
              refman

            # Move the binary to the output
            mv $TMPDIR/bin/refman $out/bin/refman
            chmod +x $out/bin/refman
          '';
        };

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
            pkgs.just
            pkgs.just-lsp
            pkgs.pre-commit
            rustToolchain
            pkgs.cargo-binstall
            refman
          ]
          ++ pkgs.lib.optional (dorado != null) dorado;

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

            # Install pre-commit hooks if not already installed
            if [ -f .pre-commit-config.yaml ]; then
              if [ ! -f .git/hooks/pre-commit ] || [ ! -f .pre-commit-installed.stamp ]; then
                echo "📝 Installing pre-commit hooks..."
                pre-commit install
                touch .pre-commit-installed.stamp
                echo "✅ Pre-commit hooks installed successfully"
              fi
            fi
          '';
        };
      }
    );
}
