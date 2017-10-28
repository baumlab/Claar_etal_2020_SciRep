import os,sys # Importing necessary libraries

for root, dirs, files in os.walk("."): # Use os.walk to walk through all files and directories within the enclosing directory
	for name in files: # Iterate through all files within the enclosing directory
		s=os.path.join(root,name) # Concatenate root directory and filename into string
		path=os.path.abspath('.') # Create string with absolute path to current folder
		t=s.split('.') # Split s on '.' to index individual chunks of the string
		n=len(t) # Make n = the length of t, useful for indexing within the array
		if t[n-1] == "fastq": # If the last chunk within array s is "fastq"; This chooses only files which are .fastq within the file structure
			v=t[1].split('_') # Split t on '_' to index individual chunks of the array
			m=len(v) # Make m = the length of v, useful for indexing within the array
			if v[m-1] == "R1" and v[m-2] == "PASSED": # If string at v[m-2] equals R1, and the string at v[m-3] equals PASSED (i.e. if this file is read 1 after boku_qc)
				name2 = "." # Initialize name2, with a period to stay consistent with remainder of output
				for x in range(0,m): # For each string in array v, 
					if x == m-1: # If string x equals m-2 ...
						name2 = name2+"_R2" # ... name2 equals name2 plus _R2 ...
					elif x == 0: # ... else if string x equals 0 ...
						name2 = name2+v[x] # ... name2 equals name2 plus v[0] ...
					else: # ... otherwise ...
						name2 = name2+"_"+v[x] # ... name2 equals name2 plus v[x]
				name2 = name2+".fastq" # Once you're done iterating through each string in array v, then add .fastq to the end of name2
				#print(name2)
				filename="" # Initialize filename
				f=t[1].split('/') # Split the first string in array t on '/' to index individual chunks of string t[1]
				# print(f)
				g=f[len(f)-1].split('_') # Split the last string in array f on '_' to index individual chunks of the last string in array f
				# print(g)
				h=g[6].split('-')
				filename=g[2]+"_"+g[3]+"_"+g[4]+"_"+g[5]+"_"+h[0] # Reconstruct filename using strings in array g
				file = open(filename+"_merge_config", "w") # Initialize new file based on filename created above concatenated with _config. "w" means write a new file, Note: will clobber files with the same name!
				file.write("[general]\nproject_name = KI_Compartment_"+filename+"\nresearcher_email = dclaar@uvic.ca\ninput_directory = "+path+"/"+root[2:]+"\noutput_directory = "+path+"/"+root[2:]+"\n\n[files]\npair_1= "+path+"/"+root[2:]+"/"+name+"\npair_2= "+path+"/"+name2[2:]) # Write the Illumina-utils config file ready for boku_qc
				file.close() # Close the file that was opened above.
