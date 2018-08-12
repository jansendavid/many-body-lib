#!/usr/bin/env python
import sys
import math
import os
import shutil
tests=['basistest', 'operatortest', 'diagtest', 'timeevtest']
for ex in tests:    
    cmd='make %s && ./%s' % (ex, ex)
    failure= os.system(cmd)
    if failure:
        print( 'expirienced error while running %s ' % ex)
        sys.exit(1)
    cmd='rm %s ' % (ex)
    failure= os.system(cmd)
    if failure:
        print( 'expirienced error while deleteing %s ' % ex)
        sys.exit(1)
 print( 'All modules of the Many Body library were successfully tested ')
