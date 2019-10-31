#!/bin/bash

# Catch user input with files list
read -p "Enter file with list of inputs: " ifiles
readarray -t filesList < $ifiles

read -p "Enter script file: " pyfile

counter=1
# Pass files as input to program
# Redirect error messages
for file in "${filesList[@]}"; do
    printf "\n File : %s\n" "$file" >> output/stoplog.txt
    python $pyfile $file || {
        printf "\n Python script failed with file \"%s\".\n" "$file" >> output/stoplog.txt
    }
    cp "gulp.got" "output/gulp${counter}.got"
    ((counter++))
done