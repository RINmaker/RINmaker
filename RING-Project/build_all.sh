#!/usr/bin/bash

rm -rf out
mkdir out
make -j4 all
echo "done."
