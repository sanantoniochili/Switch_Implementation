#!/bin/bash

usage() { echo "Usage: $0 [-s <method>] [-t <time_out>] [-e <stepsize_and_criterion>]" 1>&2; exit 1; }

# Get flags and arguments
while getopts ":s:t:o:" opt; do
  case $opt in
    s)
      echo "-s was triggered, Parameter: $OPTARG" >&2
      SFLAG="--switch"
      method=${OPTARG}
      METHOD_NM="switch"
      ;;
    t)
      echo "-t was triggered, Parameter: $OPTARG" >&2
      TFLAG="-t"
      TIMEOUT=${OPTARG}
      ;;
    o)
      echo "-o was triggered, Parameter: $OPTARG" >&2
      OFLAG="-o"
      CWD=$(pwd)
      OPTIONS="${CWD}/${OPTARG}"
      METHOD_NM="${METHOD_NM}_o"
      if [[ ! -f  $OPTIONS ]]; then
        echo "Options file not found."  # check if file exists
        exit 127
      fi
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

if [ -z "$method" ]
then
    # Method to execute
    read -p "Choose method (conj or bfgs): " method
    method="${method}"
    METHOD_NM="${method}"
fi

##################################################
################## USER INPUT ####################

# Catch user input with files list
read -p "Enter file with list of inputs : " ifiles

# Structure name is based on structure sequential number
read -p "Enter structures name : " SNAME
read -p "Enter first structure number : " counter

##################################################
################### VARIABLES ####################

# Passed arguments (options for GULP)
ARGS="$SFLAG $TFLAG $TIMEOUT $OFLAG $OPTIONS"
echo $ARGS

# IO DIRS
INPUT_DIR="input"
OUTPUT_DIR="output"

# Data DIR
DATA_DIR="/users/phd/tonyts/Desktop/Data"

# Map file
MAP="map_files.txt"

# Log file
LOG="${METHOD_NM}_stoplog.txt"

# Find file with inputs
ifiles="${DATA_DIR}/${ifiles}"
if [[ ! -f  $ifiles ]]; then
	echo "File not found."	# check if file exists
	exit 127
fi
readarray -t filesList < $ifiles

# Copy GULP IO to directories
GIN="${INPUT_DIR}/${SNAME}${counter}.gin"
GOT="${OUTPUT_DIR}/${SNAME}${counter}.got"

##################################################
################# DIRECTORIES ####################

# Appearance
BLUE='\033[0;34m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

COLUMNS=$(tput cols) 
title="DIRECTORIES" 
printf "${BLUE}%*s\n${NC}" $(((${#title}+$COLUMNS)/2)) "$title"

# Work inside new test directory
read -p "Choose test directory: " testdir
if [ ! -d "$testdir" ]; then
    echo "Creating test directory.."
    mkdir tests/$testdir
fi
TEST_DIR="tests/${testdir}/${METHOD_NM}"

CWD=$(pwd)
# Print Current directory
printf "\nInside ${CWD}. Moving to ${GREEN}${TEST_DIR}${NC}..\n"

# Check existence of method DIR
if [ ! -d "${TEST_DIR}" ]; then
    echo "Creating method directory.."
    mkdir $TEST_DIR
fi

# Copy script to method directory 
# to produce .gin, .got inside it
cp method.py $TEST_DIR/method.py
cp read_gulp.py $TEST_DIR/read_gulp.py
cd $TEST_DIR

# Check existence of IO method dirs
if [ ! -d "${INPUT_DIR}" ]; then
    echo "Creating input directory.."
    mkdir $INPUT_DIR
fi
if [ ! -d "${OUTPUT_DIR}" ]; then
    echo "Creating output directory.."
    mkdir $OUTPUT_DIR
fi

##################################################
################## EXECUTION #####################

# Empty files
> $MAP
> $LOG

title="EXECUTION" 
printf "${BLUE}%*s\n${NC}" $(((${#title}+$COLUMNS)/2)) "$title"

# Try random init
# Read every input file in files list
# Run python script and get .gin
# Run GULP and save output
for file in "${filesList[@]}"; do

    # Check if file exists
    if [ ! -f $file ]; then
      echo "File not found"
      continue
    fi

    # Log file
    printf "\n`date`\n File : %s\n" "$file" >> $LOG

    # Run GULP relaxation
    python method.py $method $file $ARGS || {
        printf "\n`date` Python script failed with file \"%s\".\n" "$file" >> $LOG
    }
    cp "gulp.gin" "${GIN}"
    cp "gulp.got" "${GOT}"

    # Map structure to initial file
    printf "${file} : ${SNAME}${counter}\n" >> $MAP

    # Add headers to csv file
    HFLAG=""
    if [[ counter -eq 1 ]]; then
      HFLAG="-c"
    fi

    # Add results to csv
    printf "Reading GULP output..."
    python read_gulp.py $GOT results.csv $METHOD_NM $HFLAG
    printf "..${GREEN}DONE${NC}\n\n"

    # Count total
    ((counter++))
done

rm gulp.gin
rm gulp.got

rm method.py
rm read_gulp.py