import os,sys # Importing necessary libraries

file = open("boku_qc_KI_Compartment.sh", "w") # Initialize new file boku_qc_KI_Compartment.sh. "w" means write a new file, Note: will clobber files with the same name!
file.write("#!/bin/bash\nset -e\nset -u\nset -o pipefail\n\n") # Write bash header to file

for file in os.listdir("."): # For each file listed in the current directory
	s="iu-filter-quality-bokulich ./" # Initialize string s with Illumina-utils command
	if file.endswith("_config"): # If the file ends with _config ...
		f=file # f equals the file name that you are currently looking at 
		file = open("boku_qc_KI_Compartment.sh", "a") # open boku_qc_KI_Compartment.sh, "a" means it will be appending to the file, not overwriting what's already there
		file.write(s+f+"\n") # write to the file s (Illumina-utils command) plus f (the file that you are currently looking at) plus \n an end of line character.

