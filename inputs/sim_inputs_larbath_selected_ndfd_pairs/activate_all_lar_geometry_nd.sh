#!/bin/bash

input_file=$2
output_file=$4

# First copy the original geometry file
cp $input_file $output_file

# Next make all the materials LAr
sed -i 's@<materialref ref=\".*\"/>@<materialref ref=\"LAr\"/>@g' $output_file

# Finally activate all regions of the detector
IFS_tmp=$IFS
IFS=" "
lines_added=0
input_str='<auxiliary auxtype="SensDet" auxvalue="SimEnergyDeposit"/>'
while IFS= read -r line_no # l=Look at all the lines that have a solidref
do
	# Get next line
	next_line_no=$((line_no + 1 + lines_added))
	next_line="$(sed -n "$next_line_no"p $output_file)"
	# The value of the next line determines where the input string should go
	# If the next line already has a SensDet, don't add anything
	if [[ $next_line =~ "<auxiliary auxtype=\"SensDet\" auxvalue=\"".*"\"/>" ]]
	then
		continue
	# If the next line is <physvol>, then the input string should be added 
	# right above the next instance of "</volume>" (and indented the same as 	 # it)
	elif [[ $next_line =~ "<physvol>" ]]
	then
		# Get line number of the next instance of "</volume>"
		vol_line_no=$(awk -v awk_line_no=$next_line_no '/<\/volume>/ { if (NR > awk_line_no) { print NR; exit } }' $output_file)

		# Add the input string right above it, at the same indentation
		vol_line="$(sed -n "$vol_line_no"p $output_file)"
		indentation_len=$(echo "$vol_line" | awk '{print match($0,/[^ ]/)-1}')
		indentation=$(printf "%*s%s" $((indentation_len - 1)) '')
		sed -i "${vol_line_no} i \ $indentation$input_str" $output_file
		
		# Increment the lines_added by 1
		lines_added=$((lines_added + 1))
	# If the next line is </volume>, add the input string right above it 
	# with two additional spaces
	elif [[ $next_line =~ "</volume>" ]]
	then
		vol_line_no=$next_line_no
		vol_line="$(sed -n "$vol_line_no"p $output_file)"
		indentation_len=$(echo "$vol_line" | awk '{print match($0,/[^ ]/)-1}')
		indentation=$(printf "%*s%s" $((indentation_len + 1)) '')
		sed -i "${vol_line_no} i \ $indentation$input_str" $output_file
		
		lines_added=$((lines_added + 1))
	# If the next line is some other auxiliary, add the string at that line 
	# with no additional spaces
	elif [[ $next_line =~ "<auxiliary auxtype=\"".*"\" auxvalue=\"".*"\"/>" ]]
	then
		vol_line_no=$next_line_no
		vol_line="$(sed -n "$vol_line_no"p $output_file)"
		indentation_len=$(echo "$vol_line" | awk '{print match($0,/[^ ]/)-1}')
		indentation=$(printf "%*s%s" $((indentation_len - 1)) '')
		sed -i "${vol_line_no} i \ $indentation$input_str" $output_file
		
		lines_added=$((lines_added + 1))
	else
		continue
	fi
done <<< $(grep -n '<solidref ref=\".*\"/>' $output_file | awk -F ":" '{print $1}')
IFS=$IFS_tmp
