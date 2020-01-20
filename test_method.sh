#!/bin/bash

usage() { echo "Usage: $0 [-s \"<from_method> <to_method>\"] [-t <time_out>] [-o <additional_options_file>]" 1>&2; exit 1; }

# Get flags and arguments
while getopts ":s:t:o:" opt; do
  case $opt in
    s)
      echo "-s was triggered, Parameter: $OPTARG" >&2
      
      set -f # disable glob
      IFS=' ' # split on space characters
      switch2=($OPTARG) # keep 2 methods

      method="${switch2[0]}"
      method2="${switch2[1]}"

      # For later operations
      SFLAG="-switch"
      METHOD_NM="switch_${method2}"
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
      usage
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
    read -p "Choose method (e.g. conj): " method
    method="${method}"
    METHOD_NM="${method}"
fi

##################################################
################## USER INPUT ####################

# Catch user input with files list
read -p "Enter file with list of inputs : " ifiles

# Structure name is based on structure sequential number
read -p "Enter samples' folder name : " SAMPLES

# Choose testing directory
read -p "Choose test folder: " testdir

##################################################
################### VARIABLES ####################

# Appearance
BLUE='\033[0;34m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color
COLUMNS=$(tput cols) 

# Passed arguments (options for GULP)
ARGS="$SFLAG $method2 $TFLAG $TIMEOUT $OFLAG $OPTIONS"

# IO DIRS
INPUT_DIR="input"
OUTPUT_DIR="output"

# Data DIR
DATA_DIR="/users/phd/tonyts/Desktop/Data"

# Map file
MAP="map_files.txt"

# Log file
LOG="${METHOD_NM}_stoplog.txt"

# Find file with list of inputs
ifiles="${DATA_DIR}/${ifiles}"
if [[ ! -f  $ifiles ]]; then
	echo "File with input list not found."	# check if file exists
	exit 127
fi
readarray -t filesList < $ifiles

##################################################
################# DIRECTORIES ####################

# Work inside new test directory
if [ ! -d "$testdir" ]; then
    echo "Creating test directory.."
    mkdir tests/$testdir
fi
TEST_DIR="tests/${testdir}/${METHOD_NM}"

# Check existence of method DIR
if [ ! -d "${TEST_DIR}" ]; then
    echo "Creating method directory.."
    mkdir $TEST_DIR
fi

CWD=$(pwd)
# Print Current directory
printf "\nInside ${CWD}. Moving to ${GREEN}${TEST_DIR}${NC}..\n"

# Copy script to method directory 
# to produce .gin, .got inside it
cp method.py $TEST_DIR/method.py
cp read_gulp.py $TEST_DIR/read_gulp.py
cd $TEST_DIR

############ INSIDE TEST DIRECTORY ###############

title="${METHOD_NM}" 
printf "${BLUE}%*s\n\n${NC}" $(((${#title}+$COLUMNS)/2)) "$title"

# Check existence of input directories
if [ ! -d "${INPUT_DIR}" ]; then
    echo "Creating input directory.."
    mkdir $INPUT_DIR
fi
if [ ! -d "${INPUT_DIR}/${SAMPLES}" ]; then
    echo "Creating samples directory.."
    mkdir $INPUT_DIR/$SAMPLES
fi

# Check existence of output directories
if [ ! -d "${OUTPUT_DIR}" ]; then
    echo "Creating output directory.."
    mkdir $OUTPUT_DIR
fi
if [ ! -d "${OUTPUT_DIR}/${SAMPLES}" ]; then
    echo "Creating samples directory.."
    mkdir $OUTPUT_DIR/$SAMPLES
fi

IN_DIR="${INPUT_DIR}/${SAMPLES}"
OUT_DIR="${OUTPUT_DIR}/${SAMPLES}"

##################################################
################## EXECUTION #####################

# Try random init
# Read every input file in files list
# Run python script and get .gin
# Run GULP and save output
counter=1
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

    # Copy GULP IO to directories
    GIN="${IN_DIR}/structure${counter}.gin"
    GOT="${OUT_DIR}/structure${counter}.got"

    cp "gulp.gin" "${GIN}"
    cp "gulp.got" "${GOT}"

    # Map structure to initial file
    printf "${file} : ${SAMPLES}/structure${counter}\n" >> $MAP

    # Add results to csv
    printf "Reading GULP output"
    python read_gulp.py $GOT results.csv $METHOD_NM
    printf "${GREEN}DONE${NC}\n\n"

    # Count total
    ((counter++))
done

title="${METHOD_NM}" 
printf "${BLUE}%*s\n\n${NC}" $(((${#title}+$COLUMNS)/2)) "$title"

if [[ -f  temp.txt ]]; then
  rm temp.txt
fi

rm gulp.gin
rm gulp.got

rm method.py
rm read_gulp.py