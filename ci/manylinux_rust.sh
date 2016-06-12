set -ex

. /io/ci/utils.sh

export TARGET=x86_64-unknown-linux-gnu
export PROJECT_NAME=convert_bng

install_rustup() {
    # uninstall the rust toolchain installed by travis, we are going to use rustup
    # sh ~/rust/lib/rustlib/uninstall.sh
    curl https://sh.rustup.rs -sSf | sh -s -- -y
    rustc -V
    cargo -V
}

# Generate artifacts for release
mk_artifacts() {
    cargo build --target $TARGET --release
}

mk_tarball() {
    # create a "staging" directory
    local td=$(echo $(mktemp -d 2>/dev/null || mktemp -d -t tmp))
    local out_dir=/io/$(pwd)

    # TODO update this part to copy the artifacts that make sense for your project
    # NOTE All Cargo build artifacts will be under the 'target/$TARGET/{debug,release}'
    for lib in target/$TARGET/release/liblonlat_bng.*; do
        strip -s $lib
    done
    cp target/$TARGET/release/liblonlat_bng.* $td

    pushd $td
    # release tarball will look like 'rust-everywhere-v1.2.3-x86_64-unknown-linux-gnu.tar.gz'
    tar czf $out_dir/${PROJECT_NAME}-${TRAVIS_TAG}-${TARGET}.tar.gz *

    popd
    rm -r $td
}

main() {
    install_rustup
    mk_artifacts
    mk_tarball
}
main
