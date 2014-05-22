##########################################

# Simple function to do base conversion 
# by "repeated division"

# That process says to convert from base
# 10 to base "b" where b < 10

# It expresses the numeral N in base b,
# the number for which will be denoted
# N_b

# An ex to convert base 10 to 2:
# 147 = 73 * 2 + 1 = c_0
# 73  = 36 * 2 + 1 = c_1 
# 36  = 18 * 2 + 0 = c_2 
# 18  = 9  * 2 + 0 = c_3 
# 9   = 4  * 2 + 1 = c_4
# 4   = 2  * 2 + 0 = c_5
# 2   = 1  * 2 + 0 = c_6
# 1   = 0  * 2 + 1 = c_7

# So that N_2 = 10010011

# This code generalizes that

##########################################



def BaseConvert(toConvert, newBase):
    if type(toConvert) != int or type(newBase) != int:
        print "Integers only please\n"
    elif newBase > 10:
        print 'New base must be less than 10\n'
    else:
        toReturn = []
        current = toConvert 
        while current != 0:
            toReturn.append(current%newBase)
            current = current / newBase
        toReturn.reverse()
        for i, x in enumerate(toReturn):
            toReturn[i] = str(toReturn[i])
        return ''.join(toReturn)


            
    


