#!/bin/bash

name=${1}

# Compile
gcc -I/users/sferrar2/data/sferrar2/CWork/include -O3 jedi${name}.c /users/sferrar2/data/sferrar2/CWork/lib/libcfitsio.a  /oscar/runtime/opt/zlib/1.2.11/lib/libz.a -o jedi${name}_s -lm -lpthread

# Move dir
cd ../jedisim/

# Copy over new file
cp ../jedisim_source/jedi${name}_s ./jedi${name}_s

# Run
#./jedi${name}_s bg138000_z09_tmp/bg138000_z09_tmp_catalog.txt bg138000_z09_tmp/bg138000_z09_tmp_distlist.txt 100000
