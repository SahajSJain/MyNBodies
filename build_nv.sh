#!/bin/bash

# Set build type
BUILD_TYPE=${1:-release}  # Default to 'release' if no argument is given

# Common settings
export FC=nvfortran

# Common flags for all build types
GPU_ARCH="-gpu=cc70,cc75,cc80,cc86,cc89,cc90"  # or specific arch like -gpu=cc80
COMMON_FFLAGS="$GPU_ARCH -mp -Minfo=accel,mp -acc -cuda -cudalib=cufft,cublas,cusolver -DP32"
COMMON_LDFLAGS="$GPU_ARCH -mp -acc -cuda -cudalib=cufft,cublas,cusolver"


# Configure flags based on build type
if [ "$BUILD_TYPE" == "debug" ]; then
    echo ">>> Building in DEBUG mode with compiler $FC for NVIDIA GPU with OpenMP+OpenACC"
    export FPM_FFLAGS="$COMMON_FFLAGS -g -O0 -traceback -Mbounds -Mchkptr -Mchkstk -Ktrap=fp -gpu=debug,lineinfo -Minfo=accel"
    export FPM_LDFLAGS="$COMMON_LDFLAGS"
elif [ "$BUILD_TYPE" == "traceback" ]; then
    echo ">>> Building in RELEASE mode with traceback with compiler $FC for NVIDIA GPU with OpenMP+OpenACC"
    export FPM_FFLAGS="$COMMON_FFLAGS -g -traceback -O3 -fast"
    export FPM_LDFLAGS="$COMMON_LDFLAGS"
elif [ "$BUILD_TYPE" == "release" ]; then
    echo ">>> Building in RELEASE mode with compiler $FC for NVIDIA GPU with OpenMP+OpenACC"
    export FPM_FFLAGS="$COMMON_FFLAGS -g -O3"
    export FPM_LDFLAGS="$COMMON_LDFLAGS"
elif [ "$BUILD_TYPE" == "profile" ]; then
    echo ">>> Building in PROFILE mode with compiler $FC for NVIDIA GPU with OpenMP+OpenACC"
    export FPM_FFLAGS="$COMMON_FFLAGS -O3 -fast -pg"
    export FPM_LDFLAGS="$COMMON_LDFLAGS -pg"
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
