#!/bin/bash

### Define functions
# Function to pad text with characters to make sentences stand out more
function string_padder () {
	string=$1
	count=${#string}
	v=$(printf "%-${count}s" "=")
	padding="${v// /=}"
	echo "

	$padding
	$string
	$padding

	"
}

finished="Finished project analysis on $(date)"
echo $finished
string_padder "$finished"
