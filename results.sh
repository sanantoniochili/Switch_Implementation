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

total=0			# random init
rat_total=0		# rattled init

counter=0		# random init
rat_counter=0	# rattled init

for file in ${DIR}/*.got; do

	# Check for optimisation
	flag=false
	if grep -q "Optimisation achieved" $file; then
		out="${file}:\tOptimisation ${GREEN}Achieved${NC}"
		flag=true
	elif grep -q "Too many failed attempts to optimise" $file; then
		out="${file}:\t${RED}Too many attempts${NC}"
	elif grep -q "Conditions for a minimum have not been satisfied" $file; then
		out="${file}:\t${RED}No lower point found${NC}"
	else
		out="${file}:\tUndefined"
	fi

	# Count success
	if [[ $file == *rat_* ]]; then	# rattled
		((rat_total++))
		if $flag; then				# success
			((rat_counter++))
		fi
	else							# random
		((total++))
		if $flag; then				# success
			((counter++))
		fi
	fi

	# Print file result
	printf "$out\n" >> $FILE
done

# Print percentages
printf "\n***** Optimised *****\n" >> $FILE
printf "Random: ${counter}/${total}\n" >> $FILE
printf "Rattled: ${rat_counter}/${rat_total}\n" >> $FILE