# This script takes care of building your crate and packaging it for release

set -ex

main() {
    local src=$(pwd) \
          stage=

    case $TRAVIS_OS_NAME in
        linux)
            stage=$(mktemp -d)
            ;;
        osx)
            stage=$(mktemp -d -t tmp)
            ;;
    esac

    test -f Cargo.lock || cargo generate-lockfile
    RUSTFLAGS='-C rpath' cross rustc --target $TARGET --features=headers --release
    if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
        for lib in target/$TARGET/release/*.so; do
            strip -s $lib
        done
        cp target/$TARGET/release/*.so $stage
    fi
    if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
        for lib in target/$TARGET/release/*.dylib; do
            strip -ur $lib
        done
        cp target/$TARGET/release/*.dylib $stage
    fi
    cp -r target/$TARGET/release/*.dSYM $stage 2>/dev/null || :
    cp include/header.h $stage
    cd $stage
    tar czf $src/$CRATE_NAME-$TRAVIS_TAG-$TARGET.tar.gz *
    cd $src
    rm -rf $stage
}

main
