#!/bin/bash

mkdir release
cd release
cmake -DCMAKE_BUILD_TYPE=Debug -S .. -B .
make
