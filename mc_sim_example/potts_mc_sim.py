#!/usr/bin/env python3
import math
from collections import Counter
import numpy as np
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
import scipy.io as sio
import random
from copy import deepcopy
from multiprocessing import Pool
import time

path='/home-3/msternk1@jhu.edu/scratch/mcmc_potts/hd/simulations/fit_200203/200916_lJ0p01_hj_gauge_transformed/'
msaFile='/home-3/msternk1@jhu.edu/scratch/mcmc_potts/hd/msa_files/hd_full_gapStrip.txt'
hWeight = 1
jWeight = 1
hFile='h_sw80_lJ0p01_gauge_transformed.npy'
jFile='Jtrans_sw80_lJ0p01.mat'
jMatName = 'Jtrans'

random.seed(4832)
nRounds = 100000
nSeqMSAs = 1000
betaSet = 0.1
betaInc = 0.2
betaStepInc = nRounds / 20
jWeight = 1
hWeight = 1
numCores = 24

def readHFile(file):
    mat = sio.loadmat(file)
    hs = mat['h']
    print(hs.shape)
    return(hs)

def readJFile(file, jName):
    mat = sio.loadmat(file)
    j_trans = mat[jName]
    j_in = np.transpose(j_trans)
    lSeq = int(-0.5*(1 - math.sqrt(8*len(j_in) + 1)))
    js = []
    jCount = 0
    for i in range(lSeq+1):
        js.append([])
        for j in range(lSeq+1):
            if i >= j:
                js[i].append(np.full((21,21),0))
            else:
                js[i].append(np.transpose(j_in[jCount][:]))
                jCount += 1
    for i in range(len(js)):
        for j in range(i+1, len(js)):
            js[j][i] = np.transpose(js[i][j])
    print(np.array(js).shape)
    return(np.array(js))

def alignmentSequences(inSequences):
    seqs = []
    testSeq= ''
    for line in inSequences:
        if line[0] == '>':
            if len(testSeq) > 0 and 'X' not in testSeq:
                seqs.append(testSeq)
            testSeq = ''
        else:
            testSeq += line.rstrip()
    seqs.append(testSeq)
    return(seqs)

def gapStrip(f,r):
    fGSnorm = []
    rGSnorm = []
    for i in range(len(f)):
        if '-' in r[i]:
            x = r[i].index('-')
            fx = f[i][x]
            del f[i][x]
            del r[i][x]
            for j in range(len(f[i])):
                f[i][j] = f[i][j] / (1 - fx)
        rGSnorm.append(r[i])
        fRow = []
        x = 0
        for j in f[i]:
            x += j
            fRow.append(x)
        fGSnorm.append(fRow)
    return([fGSnorm,rGSnorm])

def randHDSeq(res):
    seq = []
    for i in range(len(res)):
        x = random.choice(res[i])
        seq.append(x)
    return(seq)

def makeSeqSet(nSeq, res):
    seqs = []
    for i in range(nSeq):
        seqs.append(randHDSeq(res))
    return(seqs)

def mutate(seq, resList):
    while True:
        mutPosition = random.randint(0, len(seq) - 1)
        if len(resList[mutPosition]) > 1:
            break
    while True:
        mutRes = random.choice(resList[mutPosition])
        if seq[mutPosition] != mutRes:
            break
    mutResIndex = residues.index(mutRes)
    oldResIndex = residues.index(seq[mutPosition])
    return(mutPosition, mutRes, mutResIndex, oldResIndex)

def makeFasta(seqs):
    fasta = []
    for i in range(len(seqs)):
        fasta.append('>sequence' + str(i) + '\n' + ''.join(seqs[i]) + '\n')
    return(fasta)

def calcPotts(seq, hs, js, res, alphaH, alphaJ):
    hScore = 0
    jScore = 0
    for i in range(len(seq)):
        x = res.index(seq[i])
        hScore -= hs[i][x]
        for j in range(i+1,len(seq)):
            y = res.index(seq[j])
            jScore -= js[i][j][x][y]
    eScore_raw = hScore + jScore
    eScore_weighted = alphaH * hScore + alphaJ * jScore
    return((hScore, jScore, eScore_raw, eScore_weighted))

def mutationEnergy(seq, position, resIndex, j_params):
    energy = 0
    for i in range(len(seq)):
        if i != position:
            a = residues.index(seq[i])
            energy -= j_params[i][position][a][resIndex]
    return(energy)

def mcSim(num, seq):
    eDiffs = []
    
    t0 = time.time()
    
    beta = betaSet
    tempCount = 0
    seqScore = []
    maxSeq = []
    maxSeqScore = 0

    for j in range(nRounds):
        if j == 0:
            hSeq, jSeq, eSeq_raw, eSeq_weighted = calcPotts(seq,h_params,j_params,residues,hWeight, jWeight)
            seqScore.append([hSeq, jSeq, eSeq_raw, eSeq_weighted])
            maxSeq = seq
            maxSeqScore = eSeq_weighted

        else:
            mutationPos, mutationRes, mutationResIndex, oldResIndex = mutate(seq, rGS)

            newLocal = -h_params[mutationPos][mutationResIndex]
            newCoupling = mutationEnergy(seq, mutationPos, mutationResIndex, j_params)
            oldLocal = -h_params[mutationPos][oldResIndex]
            oldCoupling = mutationEnergy(seq, mutationPos, oldResIndex, j_params)
            
            dH = newLocal - oldLocal
            dJ = newCoupling - oldCoupling
            dE_raw = dH + dJ
            dE_weighted = hWeight * dH + jWeight * dJ

            eDiffs.append(dE_weighted)
            
            if dE_weighted < 0:
                seq[mutationPos] = mutationRes
                hSeq += dH
                jSeq += dJ
                eSeq_raw += dE_raw
                eSeq_weighted += dE_weighted
                
                if eSeq_weighted < maxSeqScore:
                    maxSeq = seq
                    maxSeqScore = eSeq_weighted
            else:
                pAccept = math.exp(-beta * dE_weighted)
                y = random.random()
                if y < pAccept:
                    seq[mutationPos] = mutationRes
                    hSeq += dH
                    jSeq += dJ
                    eSeq_raw += dE_raw
                    eSeq_weighted += dE_weighted
                    
            seqScore.append([hSeq, jSeq, eSeq_raw, eSeq_weighted])

        tempCount += 1
        if tempCount > betaStepInc:
            beta += betaInc
            tempCount = 0
    
    t1 = time.time()
    tTotal = t1 - t0
    
    print(f"Completed sequence {num} in {tTotal} seconds.")
    return(seqScore,seq,maxSeq,maxSeqScore,eDiffs)

if __name__ == '__main__':

    print('Running...')

    #FIGURE SETTINGS
    mpl.rcParams['axes.titlesize'] = 20
    mpl.rcParams['axes.labelsize'] = 18
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12
    mpl.rcParams['axes.facecolor'] = 'FFFFFF'
    mpl.rcParams['axes.edgecolor'] = '000000'
    mpl.rcParams['axes.linewidth'] = 1.0 
    mpl.rcParams['axes.labelweight'] = 'regular'
    mpl.rcParams['xtick.major.pad'] = 3
    mpl.rcParams['ytick.major.pad'] = 3

    h_params = np.load(hFile)
    j_params = readJFile(jFile, jMatName)
    
    with open(msaFile, 'r') as n:
        sequences = alignmentSequences(n)

    print(len(sequences))
    
    residues = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    fAll = []
    resAll = []

    for i in range(len(sequences[0])):
        fRow = []
        resRow = []
        Pi = Counter(sequence[i] for sequence in sequences)
        for z in Pi:
            resRow.append(z)
            fRow.append(Pi[z]/sum(Pi.values()))
        resAll.append(resRow)
        fAll.append(fRow)

    fGS = deepcopy(fAll)
    rGS = deepcopy(resAll)
    gs = gapStrip(fGS,rGS)
    fGS = gs[0]
    rGS = gs[1]

    seqNums = [i for i in range(1, nSeqMSAs + 1)]
    inSeqs = makeSeqSet(nSeqMSAs, rGS)
    inSeqsTup = zip(seqNums, inSeqs)
    
    p = Pool(processes=numCores)
    y1 = p.starmap(mcSim, inSeqsTup)
    p.close()

    y2 = []
    finalSeqs = []
    maxSeqs = []
    maxSeqsScores = []
    eDiffsAll = []
    for i in y1:
        y2.append(i[0])
        finalSeqs.append(i[1])
        maxSeqs.append(i[2])
        maxSeqsScores.append(i[3])
        eDiffsAll.append(i[4])

    y = np.array(y2)
    x = np.swapaxes(y,0,1)

    eDiffsArray = np.array(eDiffsAll)
    eDiffsTrans = np.swapaxes(eDiffsArray,0,1)
    eDiffMedian = np.mean(eDiffsTrans)
    print(f"Median mutation cost = {eDiffMedian}")
    
    maxInd = np.argmin(maxSeqsScores)
    maxSeq = maxSeqs[maxInd]
    maxSeqScore = maxSeqsScores[maxInd]
    
    b1 = np.arange(betaSet,betaInc*20,betaInc)
    betaTrack = np.repeat(b1,betaStepInc)

    maxSeq = ''.join(maxSeq)

    print(f">Lowest energy sequences\n{maxSeq}")
    print(maxSeqScore)

    fig, axes = plt.subplots(2,1,sharex=True,gridspec_kw={'height_ratios':[1,4]})
    axes[0].plot(betaTrack,'r-')
    axes[1].plot(x[:,:,0])
    axes[0].set_ylabel('1 / kT')
    axes[1].set_ylabel('H seq')
    axes[1].set_xlabel('MC step')
    fig.savefig(path + 'MC_htot.png', bbox_inches='tight', dpi = 400)
    plt.close()

    fig, axes = plt.subplots(2,1,sharex=True,gridspec_kw={'height_ratios':[1,4]})
    axes[0].plot(betaTrack,'r-')
    axes[1].plot(x[:,:,1])
    axes[0].set_ylabel('1 / kT')
    axes[1].set_ylabel('J seq')
    axes[1].set_xlabel('MC step')
    fig.savefig(path + 'MC_Jtot.png', bbox_inches='tight', dpi = 400)
    plt.close()

    fig, axes = plt.subplots(2,1,sharex=True,gridspec_kw={'height_ratios':[1,4]})
    axes[0].plot(betaTrack,'r-')
    axes[1].plot(x[:,:,2])
    axes[0].set_ylabel('1 / kT')
    axes[1].set_ylabel('Raw E seq')
    axes[1].set_xlabel('MC step')
    fig.savefig(path + 'MC_Etot_raw.png', bbox_inches='tight', dpi = 400)
    plt.close()
    
    fig, axes = plt.subplots(2,1,sharex=True,gridspec_kw={'height_ratios':[1,4]})
    axes[0].plot(betaTrack,'r-')
    axes[1].plot(x[:,:,3])
    axes[0].set_ylabel('1 / kT')
    axes[1].set_ylabel('E seq')
    axes[1].set_xlabel('MC step')
    fig.savefig(path + 'MC_Etot.png', bbox_inches='tight', dpi = 400)
    plt.close()

    fig, ax = plt.subplots()
    ax.hist(np.array(eDiffsTrans).flatten(), bins = 'doane', edgecolor = 'black', facecolor = 'blue', alpha = 0.5)
    ax.set_xlabel('Emut - E')
    ax.set_ylabel('Count')
    fig.savefig(path + 'e_diffs_hist.png', bbox_inches='tight',dpi = 300)
    plt.close()

    inSeqsFasta = makeFasta(inSeqs)
    outSeqsFasta = makeFasta(finalSeqs)

    with open(path + 'input_seqs.txt','w') as f:
        for i in range(len(inSeqsFasta)):
            f.write(inSeqsFasta[i])
    with open(path + 'final_seqs.txt','w') as f:
        for i in range(len(outSeqsFasta)):
            f.write(outSeqsFasta[i])

