#!/bin/bash

# IO DIRS
INPUT_DIR="input"
OUTPUT_DIR="output"

# Check existence
if [ ! -d "$INPUT_DIR" ]; then
    echo "Creating input directory.."
    mkdir input
fi
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Creating output directory.."
    mkdir output
fi

# Catch user input with files list
read -p "Enter file with list of inputs: " ifiles
ifiles="${INPUT_DIR}/${ifiles}"
if [[ ! -f  $ifiles ]]; then
	echo "File not found."	# check if file exists
	exit 127
fi
readarray -t filesList < $ifiles

# Script to execute
read -p "Choose method (conjugate or switch or BFGS): " method
pyfile="${method}.py"
if [[ ! -f $pyfile ]]; then
	echo "File not found."	# check if file exists
	exit 127
fi

# Check existence of IO method dirs
if [ ! -d "${INPUT_DIR}/${method}" ]; then
    echo "Creating directory inside input.."
    mkdir input/$method
fi
if [ ! -d "${OUTPUT_DIR}/${method}" ]; then
    echo "Creating directory inside output.."
    mkdir output/$method
fi

counter=1

# Try random init
# Read every input file in files list
# Run python script and get .gin
# Run GULP and save output
for file in "${filesList[@]}"; do

    # Log file
    LOG="${method}_stoplog.txt"
    printf "\n`date`\n File : %s\n" "$file" >> $OUTPUT_DIR/$LOG

    # Make GULP input file
    echo "Running Python script for gulp.gin.."
    python $pyfile $file || {
        printf "\n`date` Python script failed with file \"%s\".\n" "$file" >> $OUTPUT_DIR/$LOG
        # ((counter++))
        # continue
    }

    # GULP input filename
    GIN="${INPUT_DIR}/${method}/structure${counter}.gin"
    GOT="${OUTPUT_DIR}/${method}/structure${counter}.got"

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