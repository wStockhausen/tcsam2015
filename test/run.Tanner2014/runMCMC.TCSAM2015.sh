#!/bin/sh
echo on
DIR="$( cd "$( dirname "$0" )" && pwd )"
cd ${DIR}
cp ../../dist/Debug-MacOSX/CLang-MacOSX/tcsam2015 ./tcsam2015
./tcsam2015 -rs -nox  -mcmc 1000000 -mcscale 2000 -mcsave 1000 -configFile ../input.Tanner2014/TCSAM2015_ModelConfig.dat
./tcsam2015 -mceval                                            -configFile ../input.Tanner2014/TCSAM2015_ModelConfig.dat
echo off
