# CMake wrapper Makefile - preserves `make` UX for existing users

BUILD_DIR ?= build
CMAKE_ARGS ?=

all:
	cmake -B $(BUILD_DIR) $(CMAKE_ARGS)
	cmake --build $(BUILD_DIR) --parallel

clean:
	rm -rf $(BUILD_DIR)

distclean: clean
	rm -f bin/raytrace bin/trim bin/e2adc bin/instrument bin/atmosphere bin/phosimcatgen bin/version

.PHONY: all clean distclean
