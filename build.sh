#!/bin/bash

# Set build type
BUILD_TYPE=${1:-release}  # Default to 'release' if no argument is given


# Configure flags based on build type
if [ "$BUILD_TYPE" == "debug" ]; then
    export FC=mpiifx
    echo ">>> Building in DEBUG mode with compiler $FC"
    export FPM_FFLAGS="-qopenmp -g -O0 -traceback -xHost -cpp -g -O0 -traceback -check all -fpe0 -fp-stack-check -init=snan,arrays"
    export FPM_LDFLAGS="-qopenmp"
elif [ "$BUILD_TYPE" == "traceback" ]; then
    export FC=mpiifx
    echo ">>> Building in RELEASE mode with compiler $FC"
    export FPM_FFLAGS="-qopenmp -g -traceback -O3 -xHost -cpp"
    export FPM_LDFLAGS="-qopenmp "
elif [ "$BUILD_TYPE" == "release" ]; then
    export FC=mpiifx
    echo ">>> Building in RELEASE mode with compiler $FC"
    export FPM_FFLAGS="-qopenmp  -O3 -xHost -cpp"
    export FPM_LDFLAGS="-qopenmp "
else
    echo "!!! Unknown build type: $BUILD_TYPE"
    echo "!!! Usage: ./build.sh [release|debug|profile]"
    exit 1
fi

# Build using FPM with the specified compiler
fpm build  --verbose --compiler "$FC"
fpm install --prefix=./  --verbose --compiler "$FC"

# Copy the executable
cp bin/*_out .
chmod +x *_out
