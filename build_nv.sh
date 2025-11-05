#!/bin/bash

# Set build type
BUILD_TYPE=${1:-release}  # Default to 'release' if no argument is given

# Configure flags based on build type
if [ "$BUILD_TYPE" == "debug" ]; then
    export FC=nvfortran
    echo ">>> Building in DEBUG mode with compiler $FC for NVIDIA GPU with OpenMP+OpenACC"
    export FPM_FFLAGS="-acc=gpu -mp -Minfo=accel,mp -g -O0 -traceback -Mbounds -Mchkptr -Mchkstk -Ktrap=fp -DP32"
    export FPM_LDFLAGS="-acc=gpu -mp"
elif [ "$BUILD_TYPE" == "traceback" ]; then
    export FC=nvfortran
    echo ">>> Building in RELEASE mode with traceback with compiler $FC for NVIDIA GPU with OpenMP+OpenACC"
    export FPM_FFLAGS="-acc=gpu -mp -Minfo=accel,mp -g -traceback -O3 -fast -DP32"
    export FPM_LDFLAGS="-acc=gpu -mp"
elif [ "$BUILD_TYPE" == "release" ]; then
    export FC=nvfortran
    echo ">>> Building in RELEASE mode with compiler $FC for NVIDIA GPU with OpenMP+OpenACC"
    export FPM_FFLAGS="-acc=gpu -mp -Minfo=accel,mp -g -O3 -fast -DP32"
    export FPM_LDFLAGS="-acc=gpu -mp"
elif [ "$BUILD_TYPE" == "profile" ]; then
    export FC=nvfortran
    echo ">>> Building in PROFILE mode with compiler $FC for NVIDIA GPU with OpenMP+OpenACC"
    export FPM_FFLAGS="-acc=gpu -mp -Minfo=accel,mp -O3 -fast -pg -DP32"
    export FPM_LDFLAGS="-acc=gpu -mp -pg"
else
    echo "!!! Unknown build type: $BUILD_TYPE"
    echo "!!! Usage: ./build.sh [release|debug|traceback|profile]"
    exit 1
fi

# Build using FPM with the specified compiler
fpm build  --verbose --compiler "$FC"
fpm install --prefix=./  --verbose --compiler "$FC"

# Copy the executable
cp bin/*_out .
chmod +x *_out
