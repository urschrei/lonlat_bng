#!/bin/bash
set -ex

. /io/ci/utils.sh

export TARGET=x86_64-unknown-linux-gnu
export PROJECT_NAME=convert_bng
export PATH="$PATH:$HOME/.cargo/bin"
export TRAVIS_RUST_VERSION=stable

install_rustup() {
    # uninstall the rust toolchain installed by travis, we are going to use rustup
    # sh ~/rust/lib/rustlib/uninstall.sh
    curl https://sh.rustup.rs -sSf | sh -s -- -y
    rustc -V
    cargo -V
}

# Generate artifacts for release
mk_artifacts() {
    cargo build --manifest-path=/io/Cargo.toml --target $TARGET --release
}

mk_tarball() {
    # create a "staging" directory
    local td=$(echo $(mktemp -d 2>/dev/null || mktemp -d -t tmp))
    local out_dir=/io$(pwd)

    # TODO update this part to copy the artifacts that make sense for your project
    # NOTE All Cargo build artifacts will be under the 'target/$TARGET/{debug,release}'
    for lib in /io/target/$TARGET/release/liblonlat_bng.*; do
        strip -s $lib
    done
    cp /io/target/$TARGET/release/liblonlat_bng.* $td

    pushd $td
    # release tarball will look like 'rust-everywhere-v1.2.3-x86_64-unknown-linux-gnu.tar.gz'
    tar czf /io/${PROJECT_NAME}-${TRAVIS_TAG}-${TARGET}.tar.gz *

    popd
    rm -r $td
}

main() {
    install_rustup
    mk_artifacts
    mk_tarball
}
main
