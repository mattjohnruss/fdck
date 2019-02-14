#!/bin/bash
pushd build/drivers
make -j4 "$1"
strip "$1"
popd
pushd build-mingw/drivers
make -j4 "$1"
strip "$1".exe
popd
zip -j "$1".zip build/drivers/"$1" build-mingw/drivers/"$1".exe
