from __future__ import division
import functools
#import sys
import os

mydir2 = os.path.expanduser("~/GitHub/MicrobialScaling2/data/Louca/")

IN = mydir2 + 'GPC_OTU_Table_97sim.biom'

SSADs = []

OUT = open(mydir2 + 'SSADdata.txt','w+')
OUT.close()

OUT = open(mydir2 + 'SSADdata.txt','a')

c1 = '{'

row = str()
switch = 'off'

clist = ['1','2','3','4','5','6','7','8','9','0','.']

with open(IN) as f:

    f_read_ch = functools.partial(f.read, 1)
    for c in iter(f_read_ch, ''):

        #if c == c1: test block
        #    print '\nc = c1'
        #    break

        if switch == 'on':

            if c == ']' and c1 == ']':
                print>> OUT, row
                break

            if c in clist: row+=c

            elif c == ',': row+=' '

            elif c == '[': row = str()

            elif c == ']':
                print>> OUT, row
                row = str()

        elif c == '0' and c1 == '[' and switch == 'off':
            row+=c
            switch = 'on'

        c1 = c

OUT.close()
