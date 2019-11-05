# !/bin/bash

RED='\e[1;31m'
GREEN='\e[1;32m'
NC='\e[0m' # No Color

# Output DIR
OUTPUT_DIR="output"

if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Directory does not exist"
    mkdir output
fi

# Method to assess
read -p "Choose method (conjugate or switch or BFGS): " method
DIR="${OUTPUT_DIR}/${method}"

# FILE with results
FILE="${OUTPUT_DIR}/${method}_results.txt"
> $FILE

for file in ${DIR}/*.got; do
	if grep -q "Optimisation achieved" $file; then
		out="${file}: \t Optimisation ${GREEN}achieved${NC}"
	elif grep -q "Too many failed attempts to optimise" $file; then
		out="${file}: \t Optimisation ${RED}failed${NC}"
	else
		out="${file}: \t Undefined"
	fi
	printf "$out\n" >> $FILE
done
