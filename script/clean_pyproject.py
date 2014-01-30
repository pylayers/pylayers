import os
import pdb

#######################################################
#           Clean Pyrpoject directory
#
#   remove .lch from pyproject/launch
#   remove .tra from pyproject/trace
#   remove .field from pyproject/tud
#   remove .tud from pyproject/tud
#   remove .tauk from pyproject/tud
#   remove .tang from pyproject/tud
#   remove .rang from pyproject/tud
#
#   and remove [output] entries into the given 'inifile.ini' file
##########################################################


inifile='default.ini'
try:
    path=os.getenv('BASENAME')
except:
    print('Error : there is no project  directory in $BASENAME')



dirlist=['output']
extension=['.lch','.field','.tra','.tud','.tang','.rang','.tauk']
rindic=False


# remove file

for d in dirlist:
    for ex in extension:
        files = os.listdir(path +'/' +d )
        for f in files:
            if not os.path.isdir(path+d+f) and ex in f:
                rindic=True
                print f
                os.remove(path+d+f)
#                print path+d+f

        if rindic:
            print 'removed *' + ex +' from ' +d +'\n'
            rindic=False


# remove output into the where2.ini file

flagout=False
file=open(path +'/simul/' +inifile,'r')
lines=file.readlines()
file.close()
file=open(path +'/simul/' +inifile,'w')
for line in lines:
    if line =='[output]\n':
        file.write(line)
        flagout=True
    elif flagout:
        if line =='[tud]\n':
            file.write('\n')
            file.write(line)
            flagout=False
    else :
        file.write(line)
file.close()
print 'removed [output] entries into ' +inifile +'\n'
