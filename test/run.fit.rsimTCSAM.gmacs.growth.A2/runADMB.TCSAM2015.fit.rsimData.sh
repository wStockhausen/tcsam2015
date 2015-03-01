#!/bin/sh
DIR="$( cd "$( dirname "$0" )" && pwd )"
cd ${DIR}
cp ../../dist/Debug-MacOSX/CLang-MacOSX/tcsam2015 ./tcsam2015
echo ./tcsam2015 -rs -nox  -configFile Model.Config.dat -nohess -crit 0.000001
./tcsam2015 -rs -nox  -configFile Model.Config.dat -crit 0.000001 -nohess
