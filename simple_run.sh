#!/bin/bash

# Choose testing directory
read -p "Choose test environment: " testdir

# Find input folder in every method folder
for method in `ls -d -- $testdir/*/`; do
	if [ -d "${method}/input" ]; then
		# Structures' folders e.g. "random"
		for d in `ls ${method}/input`; do
			for file in `ls ${method}/input/$d`; do

				# Create missing output folders in same hierarchy
				if [ ! -d "${method}/output" ]; then
					echo "Creating ouput directory.."
				    mkdir ${method}/output
				fi
				if [ ! -d "${method}/output/$d" ]; then
					echo "Creating ${d} directory.."
				    mkdir ${method}/output
				fi

				# Get structure name
				struct="$(cut -d'.' -f1 <<< $file )"

				# Run GULP
				gulp < ${method}/input/$d/file > ${method}/output/$d/${struct}.got
				
			done
		done
	fi
done