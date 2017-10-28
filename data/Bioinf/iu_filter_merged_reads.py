import os,sys # Importing necessary libraries

file = open("iu-filter-merged-reads_KI_Compartment.sh", "w") # Initialize new file iu-filter-merged-reads_KI_Compartment.sh. "w" means write a new file, Note: will clobber files with the same name!
file.write("#!/bin/bash\nset -e\nset -u\nset -o pipefail\n\n") # Write bash header to file

for root, dirs, files in os.walk("."): # Use os.walk to walk through all files and directories within the enclosing directory
	for name in files: # Iterate through all files within the enclosing directory
		s="iu-filter-merged-reads -m 3 " # Initialize string s with Illumina-utils command
		path=os.path.abspath('.') # Create string with absolute path to current folder
		if name.endswith("MERGED"): # If the file ends with MERGED ...
			file = open("iu-filter-merged-reads_KI_Compartment.sh", "a") # open iu-filter-merged-reads_KI_Compartment.sh, "a" means it will be appending to the file, not overwriting what's already there
			file.write(s+path+"/"+root+"/"+name+"\n") # write to the file s (Illumina-utils command) plus f (the file that you are currently looking at) plus \n an end of line character.


