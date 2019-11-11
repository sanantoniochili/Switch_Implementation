#!/bin/bash

while getopts ":s:" opt; do
  case $opt in
    s)
      echo "-s was triggered, Parameter: $OPTARG" >&2
      flag="--switch"
      method=$2
      METHOD_NM="switch"
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

# IO DIRS
INPUT_DIR="input"
OUTPUT_DIR="output"

# Data DIR
DATA_DIR="/users/phd/tonyts/Desktop/Data"

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

# Check existence of method DIR
if [ ! -d "$METHOD_NM" ]; then
    echo "Creating method directory.."
    mkdir $METHOD_NM
fi

# Copy script to method directory 
# to produce .gin, .got inside it
cp method.py $METHOD_NM/method.py
cd $METHOD_NM

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
    echo "Running Python script for gulp.gin.."
    python method.py $method $file $flag || {
        printf "\n`date` Python script failed with file \"%s\".\n" "$file" >> $LOG
    }

    # GULP input filename
    GIN="${INPUT_DIR}/structure${counter}.gin"
    GOT="${OUTPUT_DIR}/structure${counter}.got"

    # GULP relaxation
    echo "Running GULP relaxation with ${GIN}.."
    cp "gulp.gin" "${GIN}"
    gulp < "${GIN}" > "${GOT}" || {
    	echo "Failed to execute GULP properly"
        exit 1
    }
    ((counter++))
done

# Try rattled init
# Read every input file in files list
# Run python script and get .gin
# Run GULP and save output
# for file in "${rattled_filesList[@]}"; do

    # Check if file exists
    if [ ! -f $file ]; then
      echo "File not found"
      continue
    fi

    # Log file
    LOG="${METHOD_NM}_stoplog.txt"
    printf "\n`date`\n File : %s\n" "$file" >> $LOG

    # Make GULP input file
    echo "Running Python script for gulp.gin.."
    python method.py $method $file $flag || {
        printf "\n`date` Python script failed with file \"%s\".\n" "$file" >> $LOG
    }

    # GULP input filename
    GIN="${INPUT_DIR}/rat_structure${counter}.gin"
    GOT="${OUTPUT_DIR}/rat_structure${counter}.got"

    # GULP relaxation
    echo "Running GULP relaxation with ${GIN}.."
    cp "gulp.gin" "${GIN}"
    gulp < "${GIN}" > "${GOT}" || {
        echo "Failed to execute GULP properly"
        exit 1
    }
    ((counter++))
done

rm gulp.gin
rm gulp.got

rm method.py