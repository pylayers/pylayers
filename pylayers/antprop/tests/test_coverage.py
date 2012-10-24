# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Initialization

# <codecell>

from pylayers.antprop.coverage import *
import time

# <markdowncell>

# Instantiate a coverage object.
# By defaut, the TA-Office.str layout strucure is loaded.
# 
# This information can be modified into the **coverage.ini** file intothe project directory.

# <codecell>

C=Coverage()
C.L.filename

# <markdowncell>

# # Coverage

# <markdowncell>

# Coverage required Layout graphs.
# 
# - If the graphs for the given layout have already been build, they are load from the project directory. 
# - Otherwise, they are build now **(it may take few minutes )** and will be stored into the project directory for the next execution.

# <codecell>

try :
    C.L.dumpr()
except:
    C.L.build()


# <markdowncell>

# Transmitter position can be modified.
# After, the coverage operation can be started and siplayed:

C.tx = np.array((39,1))
start = time.time()
C.cover()
finish = time.time()
C.showPr()
print 'Tx position: ',C.tx 
print 'Coverage in %1.2f seconds' % (finish-start)

