from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import sys
from os.path import expanduser
import random

mydir = expanduser("~/GitHub/asm")


def rounded(n):
    rem = n % 10
    if rem < 5:
        n = int(n / 10) * 10
    else:
        n = int((n + 10) / 10) * 10
    return n


def randcolor():
    r = lambda: random.randint(0,255)
    return '#%02X%02X%02X' % (r(),r(),r())


def EstimateS2(SiteList):

    """ Chao2 and ICE estimators of S for two or more samples. These metrics
    account for the occurrence (presence/absence) of taxa in a sample, but not
    the observed abundance """

    s = np.sum(np.array(SiteList), axis=0).tolist()
    S = len(s) - s.count(0)

    m = len(SiteList)

    m_inf = 0
    SpDict = {}

    for site in SiteList:

        if min(site) <= 10: m_inf += 1

        for sp in site:
            if sp in SpDict:
                SpDict[sp] += 1
            else: SpDict[sp] = 1

    IncVals = SpDict.values()

    qs = [0]*10

    for i, q in enumerate(qs):

        qs[i] = IncVals.count(i+1)

    # Chao2
    q1 = qs[0]
    q2 = qs[1]
    chao2 = S + (((m-1)/m) * ((q1*(q1-1)) / (2*(q2+1))))

    var = 'und'
    if q1 > 0 and q2 > 0:
        var = q2 * (0.5*(q1/q2)**2 + (q1/q2)**3 + 0.25*(q1/q2)**4)

    # ICE
    num = 0
    n_inf = 0

    for i, qk in enumerate(qs):
        num += (i+1)*i*qk
        n_inf += (i+1)*qk

    if m_inf-1 == 0 or n_inf**2 == 0 or 1 - (q1/n_inf) == 0:
        ice = 0
    else:
        ci = 1 - (q1/n_inf)

        gamma = (sum(qs)/ci) * (m_inf/(m_inf-1)) * (num/(n_inf**2)) - 1
        cv = max(0, gamma)

        ice = (S-sum(qs)) + (sum(qs)/ci) + ((q1/ci) * cv)

    return [chao2, var, ice, S]



################################# SET UP SYSTEM ################################

fs = 10

N = 10**5
Si = 10**3

alphas = [0, 0.999] # rows
std_maxs = [0, 10] # cols
numSites = 100

ct = 1
fig = plt.figure()
for std_max in std_maxs:
    for alpha in alphas:
        xs = []
        ys = []
        ss = []
        clrs = []

        sps = range(Si)
        for i in range(Si):
            clr = randcolor()
            clrs.append(clr)

        iclrs = []
        xmeans = np.random.uniform(0, 100, Si)
        ymeans = np.random.uniform(0, 100, Si)
        stds = np.random.uniform(0, std_max, Si)

        for i in range(N):
            s = 10**10
            while s >= Si:
                if alpha > 0:
                    s = np.random.logseries(alpha)
                else: s = np.random.randint(0, Si)

            ss.append(s)
            if std_max == 0:
                x = np.random.uniform(0, 100)
                xs.append(x)
                y = np.random.uniform(0, 100)
                ys.append(y)
            else:
                x = np.random.normal(xmeans[s], stds[s])
                xs.append(x)
                y = np.random.normal(ymeans[s], stds[s])
                ys.append(y)

            iclrs.append(clrs[s])


        indices = random.sample(range(len(xs)), int(N*0.1))
        xs2 = []
        ys2 = []
        ss2 = []
        iclrs2 = []
        for i in indices:
            xs2.append(xs[i])
            ys2.append(ys[i])
            ss2.append(ss[i])
            iclrs2.append(iclrs[i])



        ################################# START FIGURE #################################
        '''
        figure = False
        if figure == True:
            fig2 = plt.figure()
            fig2.add_subplot(2,2,1)
            plt.xlim(0,100)
            plt.ylim(0,100)
            sz = 0.1
            plt.scatter(xs, ys, s=sz, color=iclrs, alpha = 0.5, linewidth=0)
            plt.gca().set_title("Total microbiome", fontsize=fs)
            plt.tick_params(
                axis='both',          # changes apply to the both axes
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                left=False,
                labelbottom=False, # labels along the bottom edge are off
                labelleft=False)

            fig2.add_subplot(2,2,2)
            sites = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
            for x in sites:
                plt.axvline(x, ymin=0, ymax=100, color='w', linestyle='--', linewidth=1)
                plt.axhline(x, xmin=0, xmax=100, color='w', linestyle='--', linewidth=1)

            plt.xlim(0,100)
            plt.ylim(0,100)

            plt.scatter(xs2, ys2, s=sz, color=iclrs2, alpha = 0.5, linewidth=0)
            plt.gca().set_title("Sampled microbiome", fontsize=fs)
            plt.tick_params(
                axis='both',          # changes apply to the both axes
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                left=False,
                labelbottom=False, # labels along the bottom edge are off
                labelleft=False)


            plt.savefig(mydir+'/Fig1-'+str(std_max)+'-'+str(alpha)+'.png', dpi=200, bbox_inches = "tight")
            plt.close()
            sys.exit()
        '''

        ################################## END FIGURE ##################################


        ################### START SITE BY SPECIES CONSTRUCTION #########################
        s_by_s = [ [0]*Si for i in range(numSites)]

        for i, val in enumerate(xs):
            rx = rounded(val)
            ry = rounded(ys[i])

            if rx < 0 or rx > 99: continue
            if ry < 0 or ry > 99: continue
            s = ss[i]
            s_by_s[int(round(rx/10 + ry,0))][s] += 1
        s = np.sum(np.array(s_by_s), axis=0).tolist()
        St = len(s) - s.count(0)
        print ct
        print 'True S:', St

        xs = []
        ys = []
        ss = []
        iclrs = []

        s_by_s = [ [0]*Si for i in range(numSites)]

        for i, val in enumerate(xs2):
            rx = rounded(val)
            ry = rounded(ys2[i])

            if rx < 0 or rx > 99: continue
            if ry < 0 or ry > 99: continue
            s = ss2[i]
            s_by_s[int(round(rx/10 + ry,0))][s] += 1


        s = np.sum(np.array(s_by_s), axis=0).tolist()
        S = len(s) - s.count(0)
        print 'Sampled S:', S


        ICEs = []
        Chao2s = []

        samp_size = [2, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95]
        for samp in samp_size:
            I = []
            C = []
            ct2 = 0
            while ct2 < 20:
                s_by_s2 = []
                sites = random.sample(range(len(s_by_s)), samp)
                for j in sites:
                    s_by_s2.append(s_by_s[j])

                chao2, var, ice, S = EstimateS2(s_by_s2)
                if ice > 0 and np.isnan(ice) == False:
                    I.append(ice)

                    if chao2 > 0 and np.isnan(chao2) == False:
                        C.append(chao2)
                        print ct,' : ',samp,' : ',ice
                        ct2 += 1


            I = np.mean(I)
            C = np.mean(C)

            ICEs.append(int(round(I,0)))
            Chao2s.append(int(round(C,0)))

        print 'ICE:', ICEs[-1]
        print 'Chao2:', Chao2s[-1],'\n'

        fig.add_subplot(2,2,ct)

        plt.scatter(samp_size, Chao2s, s=20, color='c', label='Chao2')
        plt.scatter(samp_size, ICEs, s=20, facecolor='None', edgecolor='m', label='ICE')

        plt.axhline(y = St, color='k', linestyle='--', linewidth=1)
        plt.text(2, 1.04*St, "True richness", fontsize=fs-2)
        plt.xlim(-1,101)
        plt.ylim(0, 1.15*St)
        leg = plt.legend(loc=4,prop={'size':fs-2})
        leg.draw_frame(False)

        if ct == 1:
            plt.ylabel("Spatially uniform\n\nRichness", fontsize=fs)
            plt.gca().set_title("Even distribution of abundance\n(high detectability)\n", fontsize=fs)
            plt.xlabel("Randomly selected subplots")
        if ct == 3:
            plt.ylabel("Spatially aggregated\n\nRichness", fontsize=fs)
            plt.xlabel("Randomly selected subplots")
        if ct == 2:
            plt.gca().set_title("Uneven distribution of abundance\n(low detectability)\n", fontsize=fs)
            plt.ylabel("Richness", fontsize = fs)
            plt.xlabel("Randomly selected subplots")
        if ct == 4:
            plt.ylabel("Richness", fontsize = fs)
            plt.xlabel("Randomly selected subplots")
        plt.tick_params(axis='both', labelsize=fs-2)
        ct += 1


plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir+'/Fig2_'+str(Si)+'.png', dpi=200, bbox_inches = "tight")
plt.close()
