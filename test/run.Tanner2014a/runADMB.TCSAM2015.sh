#!/bin/sh
echo on
DIR="$( cd "$( dirname "$0" )" && pwd )"
cd ${DIR}
cp ../../dist/Debug-MacOSX/CLang-MacOSX/tcsam2015 ./tcsam2015
./tcsam2015 -rs -nox  -configFile ../input.Tanner2014a/Model.Config.dat -nohess
echo off
