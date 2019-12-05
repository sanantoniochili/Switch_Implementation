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

##################################################
################## USER INPUT ####################

# Catch user input with files list (random initialisation)
read -p "Enter file with list of inputs [random init] : " ifiles
ifiles="${DATA_DIR}/${ifiles}"
if [[ ! -f  $ifiles ]]; then
	echo "File not found."	# check if file exists
	exit 127
fi
readarray -t random_filesList < $ifiles

# Catch user input with files list (rattled initialisation)
read -p "Enter file with list of inputs [rattled init]: " r_ifiles
r_ifiles="${DATA_DIR}/${r_ifiles}"
if [[ ! -f  $r_ifiles ]]; then
    echo "File not found."  # check if file exists
    exit 127
fi
readarray -t rattled_filesList < $r_ifiles

##################################################
################# DIRECTORIES ####################

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
testdir="tests/${testdir}/${METHOD_NM}"

CWD=$(pwd)
# Print Current directory
printf "\nInside ${CWD}. Moving to ${GREEN}${testdir}${NC}..\n"

# Check existence of method DIR
if [ ! -d "${testdir}" ]; then
    echo "Creating method directory.."
    mkdir $testdir
fi

# Copy script to method directory 
# to produce .gin, .got inside it
cp method.py $testdir/method.py
cp read_gulp.py $testdir/read_gulp.py
cd $testdir

# Check existence of IO method dirs
if [ ! -d "${INPUT_DIR}" ]; then
    echo "Creating input directory.."
    mkdir $INPUT_DIR
fi
if [ ! -d "${OUTPUT_DIR}" ]; then
    echo "Creating output directory.."
    mkdir $OUTPUT_DIR
fi

counter=1

##################################################
################## EXECUTION #####################

title="RANDOM" 
printf "${BLUE}%*s\n${NC}" $(((${#title}+$COLUMNS)/2)) "$title"

# Try random init
# Read every input file in files list
# Run python script and get .gin
# Run GULP and save output
for file in "${random_filesList[@]}"; do

    # Check if file exists
    if [ ! -f $file ]; then
      echo "File not found"
      continue
    fi

    # Log file
    LOG="${METHOD_NM}_stoplog.txt"
    printf "\n`date`\n File : %s\n" "$file" >> $LOG

    # Make GULP input file
    python method.py $method $file $ARGS || {
        printf "\n`date` Python script failed with file \"%s\".\n" "$file" >> $LOG
    }

    # GULP input filename
    GIN="${INPUT_DIR}/structure${counter}.gin"
    GOT="${OUTPUT_DIR}/structure${counter}.got"

    # Map structure to initial file
    printf "${file} : structure${counter}\n" >> $MAP

    # GULP relaxation
    title="RELAXATION" 
    printf "${GREEN}%*s\n${NC}" $(((${#title}+$COLUMNS)/2)) "$title"
    cp "gulp.gin" "${GIN}"
    gulp < "${GIN}" > "${GOT}" || {
    	echo "Failed to execute GULP properly"
        exit 1
    }

    # Add headers to csv file
    HFLAG=""
    if [[ counter -eq 1 ]]; then
      HFLAG="-c"
    fi

    # Add results to csv
    printf "\nReading GULP output..."
    python read_gulp.py $GOT results.csv $METHOD_NM $HFLAG
    printf "..${GREEN}DONE${NC}"

    # Count total
    ((counter++))
done

title="RATTLED" 
printf "${BLUE}%*s\n${NC}" $(((${#title}+$COLUMNS)/2)) "$title"

# Initialise counter 
# to match no. of samples
counter=201

# Try rattled init
# Read every input file in files list
# Run python script and get .gin
# Run GULP and save output
for file in "${rattled_filesList[@]}"; do

    # Check if file exists
    if [ ! -f $file ]; then
      echo "File not found"
      continue
    fi

    # Log file
    LOG="${METHOD_NM}_stoplog.txt"
    printf "\n`date`\n File : %s\n" "$file" >> $LOG

    # Make GULP input file
    python method.py $method $file $ARGS || {
        printf "\n`date` Python script failed with file \"%s\".\n" "$file" >> $LOG
    }

    # GULP input filename
    GIN="${INPUT_DIR}/rat_structure${counter}.gin"
    GOT="${OUTPUT_DIR}/rat_structure${counter}.got"

    # Map structure to initial file
    printf "${file} : rat_structure${counter}\n" >> $MAP

    # GULP relaxation
    title="RELAXATION" 
    printf "${GREEN}%*s\n${NC}" $(((${#title}+$COLUMNS)/2)) "$title"
    echo "Running GULP relaxation with ${GIN}.."
    cp "gulp.gin" "${GIN}"
    gulp < "${GIN}" > "${GOT}" || {
        echo "Failed to execute GULP properly"
        exit 1
    }

    # Add headers to csv file
    HFLAG=""
    if [[ counter -eq 1 ]]; then
      HFLAG="-c"
    fi

    # Add results to csv
    printf "\nReading GULP output..."
    python read_gulp.py $GOT results.csv $METHOD_NM $HFLAG
    printf "..${GREEN}DONE${NC}"

    # Count total
    ((counter++))
done

rm gulp.gin
rm gulp.got

rm method.py
rm read_gulp.py