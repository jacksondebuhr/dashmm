#!/bin/bash
make clean
make
./makepoints --kernel=laplace --data=cube
./makepoints --kernel=yukawa --data=cube
./makepoints --kernel=helmholtz --data=cube
./makepoints --kernel=laplace --data=sphere
./makepoints --kernel=yukawa --data=sphere
./makepoints --kernel=helmholtz --data=sphere

./makepoints --kernel=laplace --data=cube --singlesign=yes
./makepoints --kernel=laplace --data=sphere --singlesign=yes

cd ../combinepoints
make clean
make
./combine ../makepoints/prepared.laplace.cube.dat ../makepoints/prepared.yukawa.cube.dat ../makepoints/prepared.helmholtz.cube.dat prepared.all.cube.dat
./combine ../makepoints/prepared.laplace.sphere.dat ../makepoints/prepared.yukawa.sphere.dat ../makepoints/prepared.helmholtz.sphere.dat prepared.all.sphere.dat
