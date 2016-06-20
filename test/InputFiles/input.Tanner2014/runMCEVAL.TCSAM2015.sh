#!/bin/sh
echo on
DIR="$( cd "$( dirname "$0" )" && pwd )"
cd ${DIR}
cp ../../dist/Debug-MacOSX/CLang-MacOSX/tcsam2015 ./tcsam2015
./tcsam2015 -mceval -configFile ../input.Tanner2014/Model.Config.dat
echo off
