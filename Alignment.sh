#!/bin/zsh

### THIS SCRIPT IS FOR ANALYSIS OF LENGTH MEASURED BY OPTICAL COMB LASER INTERFEROMETER ###

### INPUT FILE NAME HERE ###

FILENAME="/Users/YASUDA/Data/Alignment/171010/print_001"
FILE_ENV="/Users/YASUDA/Data/Environment/171010_174041/171010_174041"
ENV_DATE="${FILE_ENV}_date.dat"
echo ${ENV_DATE}

### INPUT PARAMETER ###

fCut=100
fCut_scan=100
frep=59.45
N=20

# ANALYZE THE FRINGE DATA TAKEN BY OSCILLOSCOPE #

echo "===FRINGE ANALYSIS==="
root -q '/Users/YASUDA/macro/Alignment/Alignment.C("'${FILENAME}'", '${fCut}', '${fCut_scan}')'

echo "===REPETITION FREQUENCY==="
echo "frep is ${frep}, N(integer) is ${N}"

echo "===ENVIRONMENT ANALYSIS==="
# echo "PLEASE INPUT THE FILENAME"
# read FILE_ENV

read env_date < ${ENV_DATE}
echo ${env_date}
root -q '/Users/YASUDA/macro/Environment/Environment.C("'${FILE_ENV}'",'${env_date}')'
