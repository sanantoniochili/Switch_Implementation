#!/bin/bash

val=0.0015
for x in {2..40}; do
	python rattled.py $val $x
	# val=$(awk "BEGIN {print $val+0.0001; exit}")
done

# for x in {1..5}; do
# 		python rattled.py $val
# 	val=$(awk "BEGIN {print $val+0.0002; exit}")
# done

# for x in {1..14}; do
# 		python rattled.py $val
# 	val=$(awk "BEGIN {print $val+0.0005; exit}")
# done

# val=0.01
# for x in {1..10}; do
# 		python rattled.py $val
# 	val=$(awk "BEGIN {print $val+0.001; exit}")
# done

# for x in {1..15}; do
# 		python rattled.py $val
# 	val=$(awk "BEGIN {print $val+0.002; exit}")
# done

# for x in {1..10}; do
# 		python rattled.py $val
# 	val=$(awk "BEGIN {print $val+0.005; exit}")
# done

# val=0.1
# for x in {1..10}; do
# 		python rattled.py $val
# 	val=$(awk "BEGIN {print $val+0.01; exit}")
# done

# for x in {1..15}; do
# 		python rattled.py $val
# 	val=$(awk "BEGIN {print $val+0.02; exit}")
# done

# for x in {1..10}; do
# 		python rattled.py $val
# 	val=$(awk "BEGIN {print $val+0.05; exit}")
# done

# python rattled.py $val
# val=$(awk "BEGIN {print $val+1; exit}")
# python rattled.py $val
