#!/bin/bash

mkdir debug
cd debug
cmake -DCMAKE_BUILD_TYPE=Debug -S .. -B .
make
