#!/bin/bash
#
# cs6900-01
# Project 01
# Ryan Miller
# w051rem
# Due 22 Jan 2021
# System = bender 
# Compiler syntax = NA
# Job Control File = NA
# Additional File   = NA
# Results file = results01.txt
#

# Color formatting declarations

errorRed="$(tput bold)$(tput setaf 1)"
boldWhite="$(tput setaf 7)$(tput bold)"
clearCyan="$(tput sgr0)$(tput setaf 6)"

# Validating that the program has been called correctly, if no inputs are provided, error

totalItems=$#
if [[ "$totalItems" == '0' ]]; then
	echo "${errorRed}Error: No Inputs Given!" >&2; exit 
fi

# Regex definitions for input identification and verification

regexInt='^[-]?[0-9]+$'
regexHelp='^\?|-h|-help|help$'
regexVersion='^-v|-version$'

# Check to see if the program was given the help flag.  If so print help page

if [[ $1 =~ $regexHelp ]]; then
  	printf "${boldWhite} ---------\n Help Page\n ---------\n ${clearCyan}Accepted Inputs: One or more space seperated integers\n Example Input: proj01.sh 1 2 3\n"
	exit
fi

# Check to see if the program was given the version flag.  If so print version page

if [[ $1 =~ $regexVersion ]]; then
  	printf "${boldWhite} ------------\n Version Info\n ------------\n${clearCyan} proj01: v1.0\n Integer Sum/Average Calculator\n Made By Ryan Miller\n"
	exit
fi

# Validate given inputs.  If invalid inputs are given, keep track of these, print them and throw error

errorArr=()
for i in "${@}"; do 
	if ! [[ $i =~ $regexInt ]]; then
  		errorArr+="($i) "
	fi
done

if [[ "${#errorArr[@]}" != '0' ]]; then
	echo "${errorRed}Error: Non-Integer Inputs: ${errorArr[@]}" >&2; exit 
fi

# Loop through the given inputs and add them together, storing the results in the sum variable
# After the sum is found, find the average of the numbers and store this in the average variable

sum=0
average=0

for i in "${@}"; do 
	sum=$(echo "$i + $sum" | bc)
done

average=$(echo "$sum / $#" | bc -l)

# Print the results to the screen, rounding to three decimal places

printf " ${boldWhite}-Sum:${clearCyan} $sum\n ${boldWhite}-Average:${clearCyan} %0.3f\n" $average

