#!/bin/bash

# Catch user input with files list
read -p "Enter file with list of inputs: " ifiles
if [[ ! -f $ifiles ]]; then
	echo "File not found"	# check if file exists
	exit 127
fi
readarray -t filesList < $ifiles

# Script to execute
read -p "Enter script file: " pyfile
if [[ ! -f $pyfile ]]; then
	echo "File not found"	# check if file exists
	exit 127
fi

counter=1
# Pass files as input to program
# Redirect error messages
for file in "${filesList[@]}"; do
    printf "\n File : %s\n" "$file" >> output/stoplog.txt
    python $pyfile $file || {
        printf "\n Python script failed with file \"%s\".\n" "$file" >> output/stoplog.txt
    }
    cp "gulp.gin" "input/CG${counter}.gin"
    gulp < "input/CG${counter}.gin" > "output/CG${counter}.got" || {
    	echo "Failed to execute GULP properly"
    }
    ((counter++))
done