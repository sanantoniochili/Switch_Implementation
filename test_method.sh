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

cp method.py $method/method.py

# # Catch user input with files list (random initialisation)
# read -p "Enter file with list of inputs [random init] : " ifiles
# ifiles="${INPUT_DIR}/${ifiles}"
# if [[ ! -f  $ifiles ]]; then
# 	echo "File not found."	# check if file exists
# 	exit 127
# fi
# readarray -t random_filesList < $ifiles

# # Catch user input with files list (rattled initialisation)
# read -p "Enter file with list of inputs [rattled init]: " r_ifiles
# r_ifiles="${INPUT_DIR}/${r_ifiles}"
# if [[ ! -f  $r_ifiles ]]; then
#     echo "File not found."  # check if file exists
#     exit 127
# fi
# readarray -t rattled_filesList < $r_ifiles

# IO DIRS
INPUT_DIR="input"
OUTPUT_DIR="output"

# Check existence of method DIR
if [ ! -d "$METHOD_NM" ]; then
    echo "Creating method directory.."
    mkdir $METHOD_NM
fi

# Check existence of IO method dirs
if [ ! -d "${METHOD_NM}/${INPUT_DIR}" ]; then
    echo "Creating input directory.."
    mkdir $METHOD_NM/input
fi
if [ ! -d "${METHOD_NM}/${OUTPUT_DIR}" ]; then
    echo "Creating output directory.."
    mkdir $METHOD_NM/output
fi

counter=1

# # Try random init
# # Read every input file in files list
# # Run python script and get .gin
# # Run GULP and save output
# for file in "${random_filesList[@]}"; do

#     # Log file
#     LOG="${METHOD_NM}_stoplog.txt"
#     printf "\n`date`\n File : %s\n" "$file" >> $OUTPUT_DIR/$LOG

#     # Make GULP input file
#     echo "Running Python script for gulp.gin.."
#     python method.py $method $file $flag || {
#         printf "\n`date` Python script failed with file \"%s\".\n" "$file" >> $OUTPUT_DIR/$LOG
#         # ((counter++))
#         # continue
#     }

#     # GULP input filename
#     GIN="${INPUT_DIR}/${METHOD_NM}/structure${counter}.gin"
#     GOT="${OUTPUT_DIR}/${METHOD_NM}/structure${counter}.got"

#     # GULP relaxation
#     echo "Running GULP relaxation with ${GIN}.."
#     cp "gulp.gin" "${GIN}"
#     gulp < "${GIN}" > "${GOT}" || {
#     	echo "Failed to execute GULP properly"
#         exit 1
#     }
#     ((counter++))
# done

# # Try rattled init
# # Read every input file in files list
# # Run python script and get .gin
# # Run GULP and save output
# # for file in "${rattled_filesList[@]}"; do

#     # Log file
#     LOG="${METHOD_NM}_stoplog.txt"
#     printf "\n`date`\n File : %s\n" "$file" >> $OUTPUT_DIR/$LOG

#     # Make GULP input file
#     echo "Running Python script for gulp.gin.."
#     python method.py $method $file $flag || {
#         printf "\n`date` Python script failed with file \"%s\".\n" "$file" >> $OUTPUT_DIR/$LOG
#         # ((counter++))
#         # continue
#     }

#     # GULP input filename
#     GIN="${INPUT_DIR}/${METHOD_NM}/rat_structure${counter}.gin"
#     GOT="${OUTPUT_DIR}/${METHOD_NM}/rat_structure${counter}.got"

#     # GULP relaxation
#     echo "Running GULP relaxation with ${GIN}.."
#     cp "gulp.gin" "${GIN}"
#     gulp < "${GIN}" > "${GOT}" || {
#         echo "Failed to execute GULP properly"
#         exit 1
#     }
#     ((counter++))
# done

# rm gulp.gin
# rm gulp.got

rm $method/method.py