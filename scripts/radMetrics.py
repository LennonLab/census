from __future__ import division
import sys
import os
import random
import numpy as np
from scipy import stats


mydir = os.path.expanduser("~/GitHub")

sys.path.append(mydir + "/tools/metrics")
import metrics
sys.path.append(mydir + "/tools/getRADs")
import getRADs


OUT = open(mydir + '/MicrobialScaling2/data/Louca/Louca-SADMetricData.txt','w+')
RADs = []

name = 'Louca'
RADs = getRADs.EMP_SADs(mydir +'/MicrobialScaling2/data/Louca', name)

print len(RADs), 'RADs'

ct = 0
numRADs = len(RADs)
for RAD in RADs:

    RAD = list([x for x in RAD if x > 0])

    N = sum(RAD)
    S = len(RAD)

    if S < 2: continue
    if max(RAD) == min(RAD): continue

    # Evenness
    Var = np.var(RAD, ddof = 1)
    Evar = metrics.e_var(RAD)
    ESimp = metrics.e_simpson(RAD)
    EQ = metrics.EQ(RAD)
    O = metrics.OE(RAD)
    #Camargo = 0.0 # metrics.camargo(RAD)   # Takes too long
    ENee = metrics.NHC(RAD)
    EPielou = metrics.e_pielou(RAD)
    EHeip = metrics.e_heip(RAD)


    # Dominance
    BP = metrics.Berger_Parker(RAD)
    SimpDom = metrics.simpsons_dom(RAD)
    Nmax = max(RAD)
    McN = metrics.McNaughton(RAD)

    # Rarity
    skew = stats.skew(RAD)
    logskew = metrics.Rlogskew(RAD)
    #p_ones = metrics.r_singletons(RAD)
    #p_zpt1 = metrics.p_ZPtOne(RAD)

    # Preston's alpha and richness, from Curtis and Sloan (2002).
    # Estimating prokaryotic diversity and its limits. PNAS.
    preston_a, preston_S = metrics.Preston(RAD)

    # Richness estimators
    chao1, ace, jknife1, jknife2 = metrics.EstimateS1(RAD)
    margalef = metrics.Margalef(RAD)
    menhinick = metrics.Menhinick(RAD)

    ct+=1

    kind = 'micro'
    print>>OUT, name, kind, N, S, Var, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S
    #print name, kind, N, S, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S

    print name, numRADs - ct
OUT.close()
