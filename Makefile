ifeq ($(shell uname),Darwin)
    EXT := dylib
else
    EXT := so
endif

all: target/debug/lonlat_bng.$(EXT)
	python src/convert.py

target/debug/lonlat_bng.$(EXT): src/lib.rs Cargo.toml
	cargo build

clean:
	rm -rf target
