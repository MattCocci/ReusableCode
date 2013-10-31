#!/usr/bin/python

################################################################
## R2Matlab ####################################################  
################################################################

# Converts an R script into Matlab code

################################################################

import sys, os, re

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
	
	## FILE AND CONTENT HANDLING #######################

	Mfile_name = "".join(split[0:len(split)-1]) + '.m'
	Rfile = open(file_list[i+1], 'r')
	Mfile = open(Mfile_name, 'w')
	code = Rfile.read()
	Rfile.close()


	## SIMPLE FIXES AND SUBSTITUTIONS #################

	R2m = {	r'#'	: r'%'	 ,  # comment char
		r'<-'	: r'='	 ,  # assignment 
		r'\[\s,': r'(:,' ,  # x[,1] to x(:,1]
		r',\s\]': r',:)' ,  # x[1,] to x[,1)
		r'\['	: r'('	 ,  # x[1] to x(1)	    
		r'\]'	: r')'	 ,
		r'NA'	: r'NaN' ,
		r'rnorm': r'randn'}

	for key in R2m:
	    code = re.sub(key, R2m[key], code)


	## LINE BY LINE FIXES #############################

	# split on lines
	code = code.split('\n')

	## FIX BROKEN LINES THAT CONTINUOUE ONTO NEXT LINE
	#   which requires ... in matlab
	code_fixed = []
	checking = ''
	for line in code:
	    checking = checking + line	

	    # If ending paren after open paren in line
	    if checking.find(')') >= checking.find('('):
		code_fixed.append(checking)
		checking = ''
	code = code_fixed


	## ALL OTHER LINE-BY-LINE FIXES
	for ind, line in enumerate(code):

	    # rep(0, n) to zeros(n,1)
	    if line.find('rep(0') != -1:
		split = line.split('rep(0,')
		start = split[0].strip()
		end = split[1].strip()
		re.sub('\)', ', 1)', end)
		new_line = start + 'zeros(' + end
		line = new_line

	    # rep(1, n) to ones(n,1)
	    if line.find('rep(1') != -1:
		split = line.split('rep(1,')
		start = split[0].strip()
		end = split[1].strip()
		re.sub('\)', ', 1)', end)
		new_line = start + 'ones(' + end
		line = new_line

#	    # randn(n) to randn(n,1)
#	    if line.find('randn(') != -1:
#		split = line.split('rep(1,')
#		start = split[0].strip()
#		end = split[1].strip()
#		re.sub('\)', ', 1)', end)
#		new_line = start + 'ones(' + end
#		line = new_line

	    # Replace variables like a.1 with a_1
	    pattern = r'[a-zA-Z_]\w*(\.\w)+'   # matches a.b, a.b.c, etc.
	    matchObj = re.search(pattern, line)
	    if matchObj != None:
		varname = matchObj.group(0)
		new_varname = re.sub(r'\.', r'_', varname)
		line = re.sub(pattern, new_varname, line)

	    # Fix for loops 
	    pattern = r'for\s*\(\s*\w+\s+in\s+\w+:\w+\s*\)\s*\{'
	    matchObj  = re.search(pattern, line)
	    if matchObj != None:
		line = re.sub(r'\(', ' ', line) 
		line = re.sub(r'\)', ' ', line) 
		line = re.sub(r'\{', '\n', line) 
		line = re.sub(r'in', '=', line) 
	    # checks that we're not replacing } in comments
	    a = line.find(r'}')
	    b = line.find(r'%')
	    if a > -1 and (b == -1 or a < b):
		line = re.sub(r'\}', '\nend', line)

	    # If statements

	    # = 1:n to seq

	    # = 1:x:y to seq

	    code[ind] = line

	code = "\n".join(code)
	Mfile.write(code)
	Mfile.close()
	

