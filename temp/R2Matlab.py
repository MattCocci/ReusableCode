#!/usr/bin/python

################################################################
## R2Matlab ####################################################  
################################################################

# Converts an R script into Matlab code

################################################################

import sys, os

# Arguments are the list of files to convert
file_list = sys.argv
file_num= len(file_list) - 1

# Loop through convertable files
for i in range(file_num):

    split = file_list[i+1].split('.')
    

    # Check for correct file type before converting
    if split[-1] != 'R':

	print file_list[i+1] + ' not an .R file'


    # Check that the R file in question exists
    elif os.path.isfile(file_list[i+1]) == False:

	print file_list[i+1] + ' is not in the current directory'

    else:
	
	# Create matlab file and open up the relevant files
	Mfile_name = "".join(split[0:len(split)-1]) + '.M'
	Rfile = open(file_list[i+1], 'r')
	Mfile = open(Mfile_name, 'w')

	grouped_code = []
	group = []
	
	# March through R code, grouping by brackets
	for line in Rfile:  
	    
	    # Check if line includes } for end of  R code block
	    #	and that it's not a comment
	    if line.find('}') != -1 and line.find('}') < line.find('#'):
		split_line = line.split('}')
		group.append(split_line[0])
		grouped_code.append(group) 
		group = [split_line[1]]

	    else:

		group.append(line)

	grouped_code.append(group)
	
	print grouped_code


	


    




	Rfile.close()
	Mfile.close()
	

