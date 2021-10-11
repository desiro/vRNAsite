#!/usr/bin/env python3
# script: dinucShuffling.py
# author: Daniel Desiro'
# dependencies: numpy, RNAfold, VARNAv3-93.jar, inkscape
script_usage="""
usage
    dinucShuffling.py -vst <vRNAsite table> -pfx <out_prefix> [options]

version
    dinucShuffling.py 0.0.1 (alpha)

dependencies
    x

description
    Validates significance of hotspot profiles in a sequence.

--prefix,-pfx
    output directory and prefix for result files

--vRNAsite,-vst
    input vRNAsite table

--candidatePeak,-cdp
    define minimum peak MFE value (default: -10.0)

--candidateMFE,-cdm
    define minimum MFE value (default: -10.0)

--candidateSPLASH,-cds
    define the minimum SPLASH count value (default: 0)

--candidateLength,-cdl
    minimum length for an interaction (default: 5)

--intra,-tra
    do intra molecular interactions from vRNAsite predictions (default: False)

--plotSize,-pls
    defines the plot size modifier (default: 1.0)

-tickScale,-tcs
    set x and y tick scale (default: 3)

--overwrite,-ovr
    overwrite data with folder named in prefix (default: False)

--criticalPercentile,-ctp
    define critical percentile for the dinucleotide shuffling Z-score (default: 0.05)

--criticalSides,-cts
    define number of sides for the dinucleotide shuffling Z-score (default: 1) (choices: 1,2)

--subsetSequences,-sub
    only use a subset of sequences for analysis; these should match the aSeq and 
    bSeq columns and should be divided by commas (default: )

--shuffles,-shf
    number of dinucleotide shuffles or random snippets (default: 1000)

--reducePlot,-rdp
    reduces total number of plot points to reduce complexity and size of the plot (default: False)

--loadShuffleData,-lsd
    load data from result folder to skip shuffling (default: False)

--seed,-sed
    set the seed for random (default: 1337)

--dangles,-dng
    use dangling ends for foldings (default: 2) (choices: 0,1,2,3)

--noLP,-nlp
    disable lonely pairs for RNAcofold (default: False)

--threads,-thr
    number of threads to use for RNAcofold (default: 1)

reference
    Reference.
"""

import argparse as ap
import sys
import os
import re
import time
import pickle
from operator import attrgetter, itemgetter
from itertools import product
from math import ceil, floor, log, isnan, log10
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from numpy import arange, mean, zeros, quantile, concatenate, median, std, cov, corrcoef, array
from scipy.stats.stats import pearsonr, spearmanr, kendalltau
import scipy.stats as st
from scipy.signal import savgol_filter
from random import seed, randint, random

# multiprocessing
try:
    import multiprocessingPickle4
    import multiprocessing as mp
    ctx = mp.get_context()
    ctx.reducer = multiprocessingPickle4.MultiprocessingPickle4()
except:
    import multiprocessing as mp
from multiprocessing import Pool

try:
    from RNA import fold_compound, cvar, CONSTRAINT_DB, CONSTRAINT_DB_DEFAULT, bp_distance
except:
    print("Error: The mandatory library ViennaRNA 2.4 is missing, please read the README.md for more information!")
    exit()


################################################################################
## main
################################################################################

def main(opt):
    ############################################################################
    ## load fasta file
    opt["var_rvp"] = False
    time_s = getTime()
    ## create output folder
    opt = makeDir(**opt)
    ############################################################################
    print(f"Status: Read tables ...")
    #tlist, segments, genomes, header, opt = readTable(**opt)
    tlist = readTable(**opt)
    time_s = getTime(time_s, f"Read table")
    ############################################################################
    ## inter and intra
    if not opt["var_tra"]:
        als = "inter"
        als_tlist = removeIntra(tlist)
    else:
        als = "intra"
        als_tlist = removeInter(tlist)
        if not als_tlist:
            print("Error: There are no intra interactions in the dataset.")
            sys.exit()
    ############################################################################
    ## get do shuffles
    print(f"Status: Do {als} shuffles ...")
    #genome_dict = doShuffles(als_tlist, segments, genomes, als, genome_dict, **opt)
    doShuffles(als_tlist, als, **opt)
    time_s = getTime(time_s, f"Do {als} shuffles")
    ############################################################################
    return opt["var_pfx"]




################################################################################
## functions
################################################################################

def readTable(**opt):
    ## read table file
    if not opt["var_tra"]: name = "inter"
    else:                  name = "intra"
    if opt["var_lsd"]:
        out_path = os.path.splitext(os.path.basename(os.path.abspath(opt["var_pfx"])))[0]
        infile = os.path.join(opt["var_pfx"], f"{out_path}_dinucleotide_{name}.tsv")
    else:
        infile = opt['var_vst']
    tlist, segments = list(), dict()
    with open(infile, "r") as tabin:
        header = next(tabin)
        header = header.strip().split()
        for i,line in enumerate(tabin):
            lnk = dict()
            line = line.strip().split()
            if len(header) != len(line):
                print(f"Error: Not the same number of header and list elements!")
                sys.exit()
            for name,item in zip(header,line):
                lnk[name] = item
            if len(lnk["aS"].split("-")) > 2: continue
            if len(lnk["bS"].split("-")) > 2: continue
            an,at = lnk["aS"].split("-")
            bn,bt = lnk["bS"].split("-")
            if at != "WT" or bt != "WT": continue
            lk = links(**lnk)
            tlist.append(lk)
    tlist = sorted(tlist, key=attrgetter("peak"), reverse=True)
    return tlist

def readTableNew(**opt):
    ## read table file
    #--orientation,-ort
    #only allow specific orientations of the vRNPs; A = arbitrary, U = upright, R = reverse (default: A) (choices: A,U,R,UR,AU,AR,AUR)
    tlist, segments = list(), dict()
    with open(opt["var_vst"], "r") as tabin:
        header = next(tabin)
        header = header.strip().split()
        for i,line in enumerate(tabin):
            lnk = dict()
            line = line.strip().split()
            if len(header) != len(line):
                print(f"Error: Not the same number of header and list elements!")
                sys.exit()
            for name,item in zip(header,line):
                lnk[name] = item
            ## check for orientation
            if "aType" not in lnk.keys():
                asplit = lnk["aSeq"].split("-")
                if len(asplit) == 2:
                    lnk["aSeq"] = asplit[0]
                    lnk["aType"] = asplit[1]
                else:
                    lnk["aSeq"] = asplit[0]
                    lnk["aType"] = "-".join(asplit[1:2])
            if "bType" not in lnk.keys():
                bsplit = lnk["bSeq"].split("-")
                if len(bsplit) == 2:
                    lnk["bSeq"] = bsplit[0]
                    lnk["bType"] = bsplit[1]
                else:
                    lnk["bSeq"] = bsplit[0]
                    lnk["bType"] = "-".join(bsplit[1:2])
            if len(lnk["aType"].split("-")) == 1: aSequence = f"{lnk['aType']}-A"
            if len(lnk["bType"].split("-")) == 1: bSequence = f"{lnk['bType']}-A"
            aOrt = aSequence.split("-")[1]
            bOrt = bSequence.split("-")[1]
            #if opt["var_ort"]:
            #    if aOrt not in list(opt["var_ort"]) or bOrt not in list(opt["var_ort"]): continue
            if aOrt != "A" or bOrt != "A": continue
            ## create class and remove items
            lk = links(**lnk)
            if lk.aj < lk.ai:
                opt["var_rvp"] = True
                lk.ai, lk.aj, lk.bi, lk.bj = lk.alen-lk.ai, lk.alen-lk.aj+1, lk.blen-lk.bi, lk.blen-lk.bj+1 # removed -1
            else:
                lk.ai, lk.aj, lk.bi, lk.bj = lk.ai-1, lk.aj, lk.bi-1, lk.bj
            if opt["var_sub"]:
                if lk.aSeq not in opt["var_sub"].split(",") or lk.bSeq not in opt["var_sub"].split(","): continue
            sta = segments.get(lk.aSeq,(0,set()))[1]
            stb = segments.get(lk.bSeq,(0,set()))[1]
            sta.add(lk.aType)
            stb.add(lk.bType)
            segments[lk.aSeq] = (lk.alen,sta)
            segments[lk.bSeq] = (lk.blen,stb)
            if type(lk.peak) is int:
                lk.peak = float(lk.peak)
            if type(lk.mfe) is int:
                lk.mfe = float(lk.mfe)
            if type(lk.sp_mean) is int:
                lk.sp_mean = float(lk.sp_mean)
            if lk.aj - lk.ai < opt["var_cdl"]: continue
            if lk.bj - lk.bi < opt["var_cdl"]: continue
            if type(lk.peak) is float:
                if lk.peak > opt["var_cdp"]: continue
            if type(lk.mfe) is float:
                if lk.mfe > opt["var_cdm"]: continue
            if "sp_mean" in lk.__dict__.keys():
                if type(lk.sp_mean) is float:
                    if lk.sp_mean < opt["var_cds"]: continue
            else:
                lk.sp_mean = float("nan")
            tlist.append(lk)
            ##########################
            #TODO: remove optional interactions, only take one with highest peak
            ##########################
    #genomes = product(*[{(s,x) for x in t} for s,(l,t) in segments])
    #segments = {s:l for s,(l,t) in segments}
    tlist = sorted(tlist, key=attrgetter("peak"), reverse=True)
    if not opt["var_tra"] and len(segments) == 0:
        print(f"Error: Not enough sequences for inter molecular interactions!")
        sys.exit()
    segments, genomes = getGenomes(segments, **opt)
    return tlist, segments, genomes, header, opt

def getGenomes(segments, **opt):
    ## get all complete genomes
    genomes = list()
    #segments["SC35M_NA"] = (1461, {"WT","mutx"})
    segments = sorted(segments.items(), key=itemgetter(1,0), reverse=True)
    itr = [{(s,x) for x in t} for s,(l,t) in segments]
    #itr = [t for s,(l,t) in segments]
    segments = {s:l for s,(l,t) in segments}
    genomes = list(product(*itr))
    return segments, genomes

class links(object):
    def __init__(self, **data):
        self.__dict__.update((k,trans(v)) for k,v in data.items())
    def trans(s):
        if isinstance(s, str):
            try: return int(s)
            except ValueError:
                try: return float(s)
                except ValueError:
                    if s in ["True", "False"]: return s == "True"
                    else: return s
        else: return s
    def plot(self, sep):
        ldat = sep.join([f"{var}" for key,var in vars(self).items()])
        return ldat

def makeDir(**opt):
    ## create directory
    dir_name, dir_base = opt["var_pfx"], opt["var_pfx"]
    if not opt["var_ovr"]:
        i = 1
        while os.path.isdir(dir_name):
            dir_name = f"{dir_base}_{i}"
            i += 1
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
    opt["var_pfx"] = dir_name
    return opt

def removeIntra(tlist):
    ## remove intra interactions
    return [tx for tx in tlist if tx.aS != tx.bS]

def removeInter(tlist):
    ## remove inter interactions
    return [tx for tx in tlist if tx.aS == tx.bS and tx.aj < tx.bi]

def getTime(time_s=0, name=""):
    if time_s and name:
        time_e = time.time()-time_s
        time_e = time.strftime("%H:%M:%S", time.gmtime(time_e))
        time_c = time.strftime('%x %X')
        print(f"Status: {name} finished at {time_c} in {time_e}")
    return time.time()

def doShuffles(tlist, als, **opt):
    ## get hotspots
    # genomes = [(("S1,"WT"),("S2,"WT"),("S3,"WT")),(("S1,"WT"),("S2,"mut"),("S3,"WT")),...]
    #for gen in genomes:
    #names = "".join([t for s,t in list(gen) if t != "WT"])
    #print(f"Status: Do {als} {names} links ...")
    #if names: names = "_" + names
    seed(opt["var_sed"])
    if opt["var_lsd"]:
        link_list = sorted(tlist, key=attrgetter("peak"), reverse=True)
    else:
        print(f"Status: Do dinucleotide shuffling ... ")
        t = len(tlist)
        pool_list = [(lk, i, t, opt) for i,lk in enumerate(tlist)]
        with Pool(processes=opt["var_thr"]) as p:
            link_list = p.starmap(getDinucQuantile, pool_list)
        link_list = sorted(link_list, key=attrgetter("peak"), reverse=True)
    plotScatter(link_list, f"{als}", **opt)
    if not opt["var_lsd"]:
        link_list = sorted(link_list, key=attrgetter("peak"), reverse=False)
        writeData(link_list, f"{als}", **opt)
    analyzeData(link_list, f"{als}_stats", **opt)

def getDinucQuantile(lk, i, t, opt):
    ## do dinucleotide shuffle and get quantile
    #aS	ai	aj	bS	bi	bj	peak	mfe	RNA	structure	a_mfe	a_structure	b_mfe	b_structure	free_mfe	free_structure
    s1,s2 = lk.RNA.split("&")
    diNuc_list = [f"{dinuclShuffle(s1)}&{dinuclShuffle(s2)}" for x in range(opt["var_shf"])]
    constraint = f"{'<'*len(s1)}{'>'*len(s2)}"
    energy_list = [doCofold(dinuc, constraint, opt) for dinuc in diNuc_list]
    #linePlot(sorted(energy_list), f"dinuc_{lk.aS}-{lk.ai}:{lk.aj}_{lk.bS}-{lk.bi}:{lk.bj}", **opt)
    #lk.qt_mfe = round(quantile(energy_list, opt["var_cfd"]),1)
    # Interval  Z
    # 80%   1.282
    # 85%   1.440
    # 90%   1.645
    # 95%   1.960
    # 99%   2.576
    # 99.5% 2.807
    # 99.9% 3.291
    m, s = mean(energy_list), std(energy_list)
    mstd = m-getZ(opt["var_ctp"], opt["var_cts"])*s
    #c1,c2 = st.norm.interval(alpha=opt["var_cfd"], loc=mean(energy_list), scale=st.sem(energy_list))
    lk.di_mfe, lk.mean_mfe, lk.std_mfe = mstd, m, s
    print(f"Status: Dinucleotide shuffling {(i/t)*100:6.2f} % ... done.                     ", end="\r")
    return lk

def getZ(a, s):
    return st.norm.ppf(1-a/s)

def getAlpha(z, s):
    return st.norm.sf(z) * s

def doCofold(RNA, constraint, opt):
    ## do Cofold
    cvar.dangles = opt["var_dng"]
    cvar.noLP = int(opt["var_nlp"])
    fc = fold_compound(RNA)
    fc.constraints_add(constraint, CONSTRAINT_DB | CONSTRAINT_DB_DEFAULT)
    pattern, mfe = fc.mfe_dimer()
    return mfe

def plotScatter(link_list, name, **opt):
    ## pyplot matrix plot
    out_path = os.path.splitext(os.path.basename(os.path.abspath(opt["var_pfx"])))[0]
    outname = os.path.join(opt["var_pfx"], f"{out_path}_dinucleotide_{name}")
    mfe_list  = [lk.mfe for lk in link_list]
    qte_list  = [lk.di_mfe for lk in link_list]
    peak_list = [lk.peak for lk in link_list]
    if opt["var_rdp"]:
        test_dict = dict()
        for m,q,p in zip(mfe_list,qte_list,peak_list):
            mx, qx = int(round(m*3,0)),int(round(q*3,0))
            px = test_dict.get((mx,qx),0.0)
            if p < px: test_dict[(mx,qx)] = p
        test_dict = {k:v for k,v in sorted(test_dict.items(), key=lambda item: item[1], reverse=True)}
        mfe_list  = [m/3  for (m,q),p in test_dict.items()]
        qte_list  = [q/3  for (m,q),p in test_dict.items()]
        peak_list = [p    for (m,q),p in test_dict.items()]
    norm = Normalize(vmin=min(peak_list), vmax=opt["var_cdp"])
    mapper = ScalarMappable(norm=norm, cmap=get_cmap("viridis"))
    color_list = [mapper.to_rgba(peak) for peak in peak_list]
    ## create plots
    s = opt["var_pls"]
    pmin, pmax = floor(min(mfe_list+qte_list)), opt["var_cdp"]+opt["var_tcs"]
    z = ceil((abs(pmin)-10) / 5)
    mx = (4.8 / z) * z + 1.6
    pz = plt.figure(figsize=(mx*s,4.8*s))
    plt.xlim([pmin,pmax])
    plt.ylim([pmin,pmax])
    plt.xticks(arange(pmin,pmax,opt["var_tcs"]))
    plt.yticks(arange(pmin,pmax,opt["var_tcs"]))
    #p1 = plt.scatter(mfe_list, qte_list, s=6, c=color_list, alpha=1.0)
    p1 = plt.scatter(mfe_list, qte_list, s=6, c=peak_list, vmin=min(peak_list), vmax=opt["var_cdp"], cmap=get_cmap('viridis'), alpha=1.0)
    plt.plot([pmin, pmax], [pmin, pmax], '-', color = 'r')
    plt.xlabel("fold mfe")
    perc = int(round(opt['var_ctp']*100,0))
    plt.ylabel(f"{perc}$^{{th}}$ percentile mfe")
    plt.colorbar(p1, label="peak mfe")
    #plt.colorbar()
    plt.subplots_adjust(bottom=0.2, left=0.2)
    #outProfile = f"{outname}_{l}_{dname}_{seg}-{typ}"
    pz.savefig(f"{outname}.svg", bbox_inches = 'tight', pad_inches = 0.1*s)
    pz.savefig(f"{outname}.pdf", bbox_inches = 'tight', pad_inches = 0.1*s)
    plt.close(pz)

def writeData(link_list, name, **opt):
    ## write scatter data
    out_path = os.path.splitext(os.path.basename(os.path.abspath(opt["var_pfx"])))[0]
    outname = os.path.join(opt["var_pfx"], f"{out_path}_dinucleotide_{name}")
    with open(f"{outname}.tsv", "w") as outdinuc:
        outdinuc.write(f"aS\tai\taj\tbS\tbi\tbj\tpeak\tmfe\tRNA\tstructure\ta_mfe\ta_structure\tb_mfe\tb_structure\tfree_mfe\tfree_structure\tdi_mfe\tmean_mfe\tstd_mfe\n")
        for lk in link_list:
            outdinuc.write(lk.plot("\t")+"\n")

def analyzeData(link_list, name, **opt):
    ## analyze data
    out_path = os.path.splitext(os.path.basename(os.path.abspath(opt["var_pfx"])))[0]
    outname = os.path.join(opt["var_pfx"], f"{out_path}_dinucleotide_{name}")
    lkn = min(link_list, key=attrgetter("peak"))
    lkx = max(link_list, key=attrgetter("peak"))
    minlk = floor(lkn.peak)
    maxlk = ceil(lkx.peak)
    p_range = {p:[0,0] for p in arange(minlk,maxlk+1,1)}
    p_range["total"] = [0,0]
    for lk in link_list:
        p = int(round(lk.peak-0.001,0))
        if lk.mfe < lk.di_mfe:
            p_range[p][0] += 1
            p_range["total"][0] += 1
        else:
            p_range[p][1] += 1
            p_range["total"][1] += 1
    with open(f"{outname}.tsv", "w") as outdinuc:
        outdinuc.write(f"peak\tbelow_qt\tabove_qt\tfraction\n")
        for p,(a,b) in p_range.items():
            frac = a/(a+b) if (a+b) != 0 else 0
            outdinuc.write(f"{p}\t{a}\t{b}\t{frac:.4f}\n")
    with open(f"{outname}_tex.tsv", "w") as outdinuc:
        outdinuc.write(f"peak\ttex\tbelow_qt\tabove_qt\tfraction\n")
        for p,(a,b) in p_range.items():
            frac = a/(a+b) if (a+b) != 0 else 0
            outdinuc.write(f"{p}\t${{}}^{{{a}}}{{\mskip -5mu/\mskip -3mu}}_{{{a+b}}}$\t{a}\t{b}\t{frac:.2f}\n")




################################################################################
## dinucleotide shuffle
################################################################################

def doEuler(s):
    dint_cnt = getCounts(s)[0]
    nts = [nt for nt in ["A","C","G","U"] if nt in s]
    edge_list = [(nt, getEdge(nt, dint_cnt)) for nt in nts if nt != s[-1]]
    isConnected = checkLastConnection(edge_list, nts, s[-1])
    return isConnected, edge_list, nts

def getCounts(s):
    dint_cnt = dict()
    nts = ["A","C","G","U"]
    nt_list  = {nt:list() for nt in nts}
    dint_cnt = {(a,b):0 for a in nts for b in nts}
    dint_test = 1
    for i in range(len(s)-1):
        a, b = s[i], s[i+1]
        nt_list[a].append(b)
        dint_cnt[(a,b)] += 1; dint_test += 1
    return dint_cnt, nt_list

def getEdge(a, dint_cnt):
    rnd = random()
    nts = ["A","C","G","U"]
    nt_sum = sum([dint_cnt[(a,b)] for b in nts])
    nt_num, nt = 0, "U"
    for b in nts:
        nt_num += dint_cnt[(a,b)]
        if rnd < float(nt_num)/float(nt_sum):
            nt = b
            break
    dint_cnt[(a,nt)] -= 1
    return nt

def checkLastConnection(edge_list, nts, last_nt):
    last = {n:0 for n in nts}
    for a,b in edge_list:
        if b == last_nt: last[a] = 1
    for i in range(2):
        for a,b in edge_list:
            if last[b] == 1: last[a] = 1
    for n in nts:
        if n != last_nt and last[n] == 0:
            return False
    return True

def shuffleEdges(shuffle_list):
    end, n = len(shuffle_list), len(shuffle_list)
    for i in range(n-1):
        rnd = int(random() * end)
        tmp  = shuffle_list[rnd]
        shuffle_list[rnd] = shuffle_list[end-1]
        shuffle_list[end-1] = tmp
        end -= 1
    return shuffle_list

def dinuclShuffle(s):
    s = s.upper().replace("T","U")
    isConnected = False
    while not isConnected:
        isConnected, edge_list, nts = doEuler(s)
    nt_list = getCounts(s)[1]
    for a,b in edge_list: nt_list[a].remove(b)
    for a in nts: shuffleEdges(nt_list[a])
    for a,b in edge_list: nt_list[a].append(b)
    last_nt, shuffle_list = s[0], [s[0]]
    for i in range(len(s)-2):
        nt = nt_list[last_nt][0]
        shuffle_list.append(nt)
        del nt_list[last_nt][0]
        last_nt = nt
    shuffle_list.append(s[-1])
    dinuc = "".join(shuffle_list)
    return dinuc




################################################################################
## parser
################################################################################

if __name__ == "__main__":

    ############################################################################
    ## get time and save call
    sscript = sys.argv[0]
    start_time = time.time()
    current_time = time.strftime('%x %X')
    scall = " ".join(sys.argv[1:])
    with open(f"{sscript}.log", "a") as calllog:
        calllog.write(f"Start : {current_time}\n")
        calllog.write(f"Script: {sscript}\n")
        calllog.write(f"Call  : {scall}\n")
    print(f"Call: {scall}")
    print(f"Status: Started at {current_time}")
    ############################################################################
    ## transform string into int, float, bool if possible
    def trans(s):
        if isinstance(s, str):
            try: return int(s)
            except ValueError:
                try: return float(s)
                except ValueError:
                    if s in ["True", "False"]: return s == "True"
                    else: return s
        else: return s
    ############################################################################
    ## save documentation
    rx_text = re.compile(r"\n^(.+?)\n((?:.+\n)+)",re.MULTILINE)
    rx_oneline = re.compile(r"\n+")
    rx_options = re.compile(r"\((.+?)\:(.+?)\)")
    help_dict, type_dict, text_dict, mand_list = {}, {}, {}, []
    for match in rx_text.finditer(script_usage):
        argument = match.groups()[0].strip()
        text = " ".join(rx_oneline.sub("",match.groups()[1].strip()).split())
        argopts = {"action":"store", "help":None, "default":None, "choices":None}
        for option in rx_options.finditer(text):
            key = option.group(1).strip()
            var = option.group(2).strip()
            if var == "False": argopts["action"] = "store_true"
            if var == "True": argopts["action"] = "store_false"
            if key == "choices": var = [vs.strip() for vs in var.split(",")]
            if key == "default": var = trans(var)
            argopts[key] = var
        if argopts["default"]: add_default = f" (default: {str(argopts['default'])})"
        else: add_default = ""
        argopts["help"] = rx_options.sub("",text).strip()+add_default
        argnames = argument.split(",")
        if len(argnames) > 1:
            if argopts["default"] == None:
                mand_list.append(f"var_{argnames[1][1:]}")
            type_dict[f"var_{argnames[1][1:]}"] = argopts["default"]
            argopts["argshort"] = argnames[1]
            help_dict[argnames[0]] = argopts
        else:
            text_dict[argnames[0]] = argopts["help"]
    ############################################################################
    ## get arguments
    if text_dict["dependencies"]:
        desc = f"{text_dict['description']} (dependencies: {text_dict['dependencies']})"
    else:
        desc = text_dict['description']
    p = ap.ArgumentParser(prog=sscript, prefix_chars="-", usage=text_dict["usage"],
                          description=desc, epilog=text_dict["reference"])
    p.add_argument("-v", "--version", action="version", version=text_dict["version"])
    for argname,argopts in help_dict.items():
        argshort = argopts["argshort"]
        if argopts["choices"]:
            p.add_argument(argshort, argname,            dest=f"var_{argshort[1:]}",\
                           action=argopts["action"],     help=argopts["help"],\
                           default=argopts["default"],   choices=argopts["choices"])
        else:
            p.add_argument(argopts["argshort"], argname, dest=f"var_{argshort[1:]}",\
                           action=argopts["action"],     help=argopts["help"],\
                           default=argopts["default"])
    p._optionals.title = "arguments"
    opt = vars(p.parse_args())
    ############################################################################
    ## validate arguments
    if None in [opt[mand] for mand in mand_list]:
        print("Error: Mandatory arguments missing!")
        print(f"Usage: {text_dict['usage']} use -h or --help for more information.")
        sys.exit()
    for key,var in opt.items():
        if key not in mand_list:
            arg_req, arg_in = type_dict[key], trans(var)
            if type(arg_req) == type(arg_in):
                opt[key] = arg_in
            else:
                print(f"Error: Argument {key} is not of type {type(arg_req)}!")
                sys.exit()
    ############################################################################
    ## call main function
    try:
        #saved = main(**opt)
        saved = main(opt)
    except KeyboardInterrupt:
        print("Error: Interrupted by user!")
        sys.exit()
    except SystemExit:
        print("Error: System exit!")
        sys.exit()
    except Exception:
        print("Error: Script exception!")
        traceback.print_exc(file=sys.stderr)
        sys.exit()
    ############################################################################
    ## finish
    started_time = current_time
    elapsed_time = time.time()-start_time
    elapsed_time = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    current_time = time.strftime('%x %X')
    if saved:
        with open(f"{sscript}.log", "a") as calllog,\
             open(os.path.join(saved,f"call.log"), "a") as dirlog:
            calllog.write(f"Save  : {os.path.abspath(saved)}\n")
            calllog.write(f"Finish: {current_time} in {elapsed_time}\n")
            ## dirlog
            dirlog.write(f"Start : {started_time}\n")
            dirlog.write(f"Script: {sscript}\n")
            dirlog.write(f"Call  : {scall}\n")
            dirlog.write(f"Save  : {os.path.abspath(saved)}\n")
            dirlog.write(f"Finish: {current_time} in {elapsed_time}\n")
    print(f"Status: Saved at {saved}")
    print(f"Status: Finished at {current_time} in {elapsed_time}")
    sys.exit(0)
