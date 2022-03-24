#!/usr/bin/env python3
# script: vRNAsite.py
# author: Daniel Desiro'
script_usage="""
usage
    vRNAsite.py -fsa <in_fasta> -pfx <out_prefix> [options]

version
    vRNAsite.py 0.0.1 (alpha)

dependencies
    Python v3.9.7, NumPy v1.22.2, Pandas v1.4.0, Bokeh v2.4.2, 
    ViennaRNA v2.5.0, Matplotlib v3.5.1, VARNA v3.93, Circos v0.69.8

description
    Predicts potential intermolecular long-range RNA-RNA interactions between 
    two or more RNA sequences. Can also be applied to predict intramolecular 
    long-range RNA-RNA interactions. The tool predicts short consecutive 
    stable interactions. Example call: python vRNAsite.py -pfx example -fsa 
    example.fa -thr 4 -ovr -rev -cmp -rvp -nex peak -clp -13.0

################################################################

--prefix,-pfx
    output directory and prefix for result files

--fasta,-fsa
    fasta file with all mutants; one entry for each RNA from WT and all
    mutants; will create all possible combinations between mutants and WT;
    naming: \">name WT\" or \">name mut_name:start-end\"; start and end
    defines which part in the WT to be replaced with the mutant string

--reverse,-rev
    creates reverse of each strain if set (default: False)

--complement,-cmp
    creates complements of each strain if set (default: False)

--matrix,-mat
    load pickled matrix object from a previous run (default: False)

--threads,-thr
    number of threads to use for RNAcofold (default: 1)

--intra,-tra
    also do intra RNA combinations (default: False)

--slice,-slc
    window size for standard window search or minimum length for a seed to be
    accepted if the seed option is specified (default: 20)

--step,-stp
    step size for window search, mainly for debugging, should be always set 
    to 1 (default: 1)

--overwrite,-ovr
    overwrite data with folder named in prefix (default: False)

--dangles,-dng
    use dangling ends for foldings (default: 2) (choices: 0,1,2,3)

--temperature,-tmp
    changes the RNA folding temperature (default: 37.0)

--noLP,-nlp
    disable lonely pairs for RNAcofold (default: False)

--bigGenome,-big
    use this option for an alternative multiprocessing algorithm for big genomes (default: False)

--onlyPlotting,-opl
    only to plotting, requires the candidates.pcl and the matrices.pcl to be
    in the prefix directory (default: False)

--reversePositions,-rvp
    reverse all positions, useful for negative stranded RNA (default: False)

--namingExtension,-nex
    use the this parameter as an extension for naming (default: none) (choices: none,peak,mfe,)

################################################################

--splashData,-spd
    optionally read SPLASH, SHAPE or similar data from matrices in directory, names have to match the fasta names (default: )

--splashReads,-spr
    set mean minimum threshold for SPLASH, SHAPE or similar data (default: 0.0)

################################################################

--weightMatrix,-wgm
    use weight matrix to consider the relative positioning (default: False)

--weightDistance,-dst
    distance threshold modifier for links to be in the same area from the center
    of their respective RNAs (default: 0.1)

--weightDescent,-dsc
    defines how the weight should descent (default: 0.0)

--saveOrientation,-svo
    also save orientated matrices, this will use up more space (default: False)

################################################################

--clusterPeak,-clp
    maximum peak energy for clusters (default: -10.0)

--clusterBarrier,-clb
    maximum barrier energy for cluster (default: -10.0)

--clusterValue,-clv
    sets the cluster extraction value for plotting (default: 1.0)

################################################################

--candidateDangling,-cdd
    trimming will retain a minimum amount of free bases (default: 2)

--candidateProcesses,-cdp
    turn off multi processing for candidate extraction (default: True)

--candidateLoad,-cdl
    load candidates from file to reduce extraction calculations (default: False)

--candidateCutoff,-cdc
    maximum free energy for a structure to be accepted (default: -10.0)

--candidateMinlen,-cdm
    define minimum length of RNA interaction (default: 8)

################################################################

--varnaFold,-vrf
    make VARNA plots (default: False)

--varnaPython,-vrh
    make VARNA python scripts (default: False)

--varnaSingle,-vrs
    also do single sequence VARNA plots (default: False)

--varnaAlgorithm,-vra
    defines the VARNA drawing algorithm (default: radiate) (choices: line,circular,radiate,naview)

--varnaPeakEnergy,-vre
    maximum peak average free energy for a structure to be plotted with VARNA 
    (default: -10.0)

--varnaMinimumEnergy,-vrm
    maximum minimum free energy for a structure to be plotted with VARNA 
    (default: -10.0)

--varnaPath,-vrp
    use this VARNA path; example: VARNAv3-93.jar 
    (default: VARNAv3-93.jar)

--varnaPdf,-vrd
    use this Inkscape path to create pdf files; example: inkscape
    (default: )

################################################################

--plotHeats,-plh
    plot interaction heat maps (default: False)

--plotProcesses,-plp
    turn on multi processing for matrix plotting (default: False)

--plotSize,-pls
    defines the plot size modifier (default: 1.0)

--plotMax,-plm
    set maximum heat map plot energy for candidate plots (default: -10.0)

--plotStepX,-plx
    step size for x axis ticks (default: 100)

--plotStepY,-ply
    step size for y axis ticks (default: 100)

--plotSvg,-plv
    also plot svg heat maps (default: False)

--plotColor,-plc
    reverse plot color (default: False)

--plotNan,-pln
    plot nan values, else they will be set to zero (default: False)

################################################################

--bokehHeats,-bkh
    enables bokeh printing (default: False)

--bokehProcesses,-bkp
    turn off multi processing for bokeh plotting (default: True)

--bokehSvg,-bkv
    use svg option for bokeh plots, interactive plots will run less smooth
    (default: False)

--bokehSize,-bks
    defines the bokeh plot size modifier (default: 1.0)

--bokehMax,-bkm
    set maximum bokeh plot energy (default: -10.0)

################################################################

--distributionPlots,-dsh
    enables peak energy distribution printing (default: False)

--distributionSize,-dss
    defines the distribution plot size modifier (default: 1.0)

--distributionSvg,-dsv
    also plot svg distribution plots (default: False)

################################################################

--circosPlots,-cih
    enables circos printing (default: False)

--circosPeakEnergy,-cie
    maximum peak average free energy for a link to be plotted with circos 
    (default: -10.0)

--circosMinimumEnergy,-cim
    maximum minimum free energy for a link to be plotted with circos 
    (default: -10.0)

--circosPath,-cip
    use alternative circos binary (default: circos)

--circosSingle,-cis
    create plots for every segment (default: False)

-circosIntraInter,-cit
    specify the interaction type which should be shown with circos (default: all) 
    (choices: all,intra,inter)

-circosIntraDist,-cid
    define the minimum distance of intra circos intreactions (default: 10)

--circosRange,-cir
    set circos position plot range for one sequence; sequence name, start and 
    end position has to be divided by a minus symbol; e.g. <name>-<start>-<end> 
    (default: )

--circosMutants,-ciu
    create mutant circos plots; the mutant segments can be defined by this 
    parameter; sequences name and mutant name have to be separated by a comma, 
    multiple mutants by semicolon; e.g. <name1>,<mut1>;<name2>,<mut2> (default: )

################################################################

reference
    D. Desirò, E. Barth, B. Hardin, M. Schwemmle and M. Marz.
    "vRNAsite: Prediction of packaging signals in segmented RNA viruses at the example of influenza A virus indicating flexible RNA-RNA interactions between segments."
    In Preparation, 2022.
    https://github.com/desiro/vRNAsite
"""

import argparse as ap
import sys
import os
import re
import traceback
import time
import pickle
import operator
try:
    from RNA import fold_compound, cvar, CONSTRAINT_DB, CONSTRAINT_DB_DEFAULT, bp_distance
except:
    print("Error: The mandatory library ViennaRNA 2.4 is missing, please read the README.md for more information!")
    exit()
from itertools import combinations, product, combinations_with_replacement
from numpy import arange, mean, median, prod, array, matrix, zeros, ones, nonzero, ones, linspace, transpose, concatenate, flip, full, nan_to_num
from math import ceil, floor, log, isnan
from subprocess import Popen, PIPE, call
from random import seed, randint

try:
    import multiprocessingPickle4
    import multiprocessing as mp
    ctx = mp.get_context()
    ctx.reducer = multiprocessingPickle4.MultiprocessingPickle4()
except:
    import multiprocessing as mp

# plotting
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from pandas import DataFrame
from bokeh.plotting import figure, output_file, show, save, ColumnDataSource
from bokeh.models import HoverTool, Legend, Plot, Range1d, Label, Slider, CustomJS, ColorBar, LinearColorMapper, BasicTicker
from bokeh.models.tools import *
from bokeh.models.glyphs import Text
from bokeh.models.widgets import RadioButtonGroup
from bokeh.layouts import column, row, gridplot, layout, widgetbox
from bokeh.io import export_svgs


################################################################################
## main
################################################################################

def main(opt):
    time_s = getTime()
    ############################################################################
    ## create output folder
    print(f"Status: Create output directory ...")
    opt = makeDir(**opt)
    makeTemp(**opt)
    time_s = getTime(time_s, f"Create output directory")
    ########################################################################
    ## read fasta file
    print(f"Status: Read fasta file ...")
    data_dict, mut_dict = readFasta(**opt)
    time_s = getTime(time_s, f"Read fasta file")
    ############################################################################
    ## create WT / mutant combinations
    print(f"Status: Create WT / mutant combinations ...")
    comb_list = defineCombinations(data_dict, **opt)
    time_s, opt["var_tem"] = getTime(time_s, f"Create WT / mutant combinations"), 50
    ############################################################################
    if opt["var_mat"] or opt["var_opl"]:
        ## load matrices
        print(f"Status: Load matrices ...")
        heat_dict, comb_list, splash_dict = loadMatrix(comb_list, **opt)
        time_s = getTime(time_s, f"Load matrices")
    else:
        heat_dict = dict()
        splash_dict = dict()
    ############################################################################
    if opt["var_spd"] and not opt["var_mat"] and not opt["var_opl"]:
        ## read SPLASH matrices and transform to heat_dict
        print(f"Status: Read SPLASH matrices ...")
        splash_dict = loadSPLASH(comb_list, data_dict, **opt)
        time_s = getTime(time_s, f"Read SPLASH matrices")
    ############################################################################
    if opt["var_cdl"] or opt["var_opl"]:
        ## load candidates
        print(f"Status: Load candidates ...")
        cand_list = loadCandidates(**opt)
        time_s = getTime(time_s, f"Load candidates")
    else: cand_list = list()
    ############################################################################
    if comb_list and not opt["var_opl"]:
        ## create sequence snippets
        print(f"Status: Create sequence snippets ...")
        snip_dict = createSnippets(data_dict, **opt)
        time_s = getTime(time_s, f"Create sequence snippets")
        ########################################################################
        ## combine sequence snippets
        print(f"Status: Combine sequence snippets ...")
        fold_list = combineSnippets(comb_list, data_dict, snip_dict, mut_dict, **opt)
        time_s = getTime(time_s, f"Combine sequence snippets")
        ########################################################################
        ## fold combinations and create matrices
        print(f"Status: Fold snippets and create matrices ...")
        if opt["var_big"]:
            res_list = foldCombinationsBig(fold_list, data_dict, **opt)
        else:
            res_list = foldCombinations(fold_list, data_dict, **opt)
        time_s = getTime(time_s, f"Fold snippets and create matrices")
        ########################################################################
        ## create matrices
        print(f"Status: Create matrices ...")
        heat_dict = createMatrices(res_list, heat_dict, mut_dict, comb_list, **opt)
        time_s = getTime(time_s, f"Create matrices")
    ############################################################################
    if opt["var_wgm"] and not opt["var_opl"]:
        ## orientate matrices
        print(f"Status: Orientate matrices ...")
        heat_dict = orientateMatrices(heat_dict, comb_list, **opt)
        time_s = getTime(time_s, f"Orientate matrices")
    ############################################################################
    ## save matrices
    print(f"Status: Save matrices ...")
    saveMatrix(heat_dict, splash_dict, **opt)
    time_s = getTime(time_s, f"Save matrices")
    ############################################################################
    if opt["var_plh"]:
        ## plot matrices
        print(f"Status: Plot matrices ...")
        plotMatrices(heat_dict, "clust_plots", **opt)
        opt["var_plm"] = 1.0
        plotMatrices(heat_dict, "heat_plots", **opt)
        if splash_dict:
            plm = opt["var_plm"]
            opt["var_plm"] = 0
            plotMatrices(splash_dict, "splash_plots", **opt)
            opt["var_plm"] = plm
        time_s = getTime(time_s, f"Plot matrices")
    ############################################################################
    if opt["var_bkh"]:
        ## plot bokeh
        print(f"Status: Plot bokeh ...")
        plotBokeh(heat_dict, data_dict, splash_dict, **opt)
        time_s = getTime(time_s, f"Plot bokeh")
    ############################################################################
    if not opt["var_opl"]:
        ## extract candidates
        print(f"Status: Extract candidates  ...")
        cand_list = extractCandidates(heat_dict, data_dict, cand_list, **opt)
        time_s = getTime(time_s, f"Extract candidates")
        ########################################################################
        ## extract candidates
        if opt["var_spd"]: print(f"Status: Add SPLASH data  ...")
        cand_list = addSplash(cand_list, splash_dict, **opt)
        if opt["var_spd"]: time_s = getTime(time_s, f"Add SPLASH data")
        ########################################################################
        ## save matrices
        print(f"Status: Save interactions ...")
        saveCandidates(cand_list, **opt)
        time_s = getTime(time_s, f"Save interactions")
    ############################################################################
    ## create circos plots
    if opt["var_dsh"]:
        print(f"Status: Create distribution plots ...")
        createDistributionPlot(cand_list, data_dict, **opt)
        time_s = getTime(time_s, f"Plot distributions")
    ############################################################################
    ## create circos plots
    if opt["var_cih"]:
        print(f"Status: Create circos plots ...")
        createCircosData(cand_list, data_dict, **opt)
        time_s = getTime(time_s, f"Plot circos")
    ############################################################################
    if opt["var_vrf"] or opt["var_vrh"]:
        ## plot structures
        print(f"Status: Plot structures ...")
        createStructures(cand_list, **opt)
        time_s = getTime(time_s, f"Plot structures")
    ############################################################################
    ## remove temp dir
    removeTemp(**opt)
    return opt["var_pfx"]




################################################################################
## functions
################################################################################

def getTime(time_s=0, name=""):
    if time_s and name:
        time_e = time.time()-time_s
        time_e = time.strftime("%H:%M:%S", time.gmtime(time_e))
        time_c = time.strftime('%x %X')
        print(f"Status: {name} finished at {time_c} in {time_e}")
    return time.time()

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

def makeTemp(**opt):
    ## create temp dir
    if not os.path.isdir(os.path.join(opt["var_pfx"],"vRNAsite_temp")):
        os.mkdir(os.path.join(opt["var_pfx"],"vRNAsite_temp"))

def removeTemp(**opt):
    ## remove temp dir
    if os.path.isdir(os.path.join(opt["var_pfx"],"vRNAsite_temp")):
        os.rmdir(os.path.join(opt["var_pfx"],"vRNAsite_temp"))

def readFasta(**opt):
    ## read fasta file
    ## mut_dict[(name,strain)] = (start,end)
    ## data_dict[(name,strain)] = RNA
    data_dict, mut_dict = dict(), dict()
    RNA = ""
    with open(opt["var_fsa"], "r") as infa:
        for line in infa:
            line = line.strip()
            if re.match(r">", line):
                if RNA:
                    if strain != "WT": RNA, strain, mut_dict = createMut(data_dict, mut_dict, name, RNA, strain, **opt)
                    data_dict[(name, strain)] = revComp(RNA, **opt)
                try:
                    name, strain = line[1:].split()[0], line[1:].split()[1]
                except ValueError: 
                    print(f"Error: \"{line}\" is not in shape \">vRNA# WT\" or \">vRNA# mut_name:start-end\"!")
                    sys.exit()
                RNA = ""
            else:
                RNA += line
        if strain != "WT": RNA, strain, mut_dict = createMut(data_dict, mut_dict, name, RNA, strain, **opt)
        data_dict[(name, strain)] = revComp(RNA, **opt)
    createFasta(data_dict, **opt)
    return data_dict, mut_dict

def createFasta(data_dict, **opt):
    ################################################################################
    ## create reverse complement fasta file
    ftype = ""
    if opt["var_rev"] or opt["var_cmp"]:
        if opt["var_rev"]: ftype += "_rev"
        if opt["var_cmp"]: ftype += "_cmp"
        with open(os.path.join(opt["var_pfx"], f"fasta{ftype}.fa"), "w") as outfa:
            for (name,strain),RNA in sorted(data_dict.items()):
                outfa.write(f">{name}-{strain}\n{RNA}\n")

def createMut(data_dict, mut_dict, name, RNA, strain, **opt):
    ## create mutant
    ## mut_dict[(name,strain)] = (start,end)
    strain, mrange = strain.split(":")
    start, end = mrange.split("-")
    WT_RNA = data_dict.get((name, "WT"), "")
    if not WT_RNA:
        print("Error: mutant defined in fasta before WT!")
        sys.exit()
    WT_RNA = revComp(WT_RNA, **opt)
    WTlen = len(WT_RNA)
    RNA = WT_RNA[:int(start)-1]+RNA+WT_RNA[int(end):]
    if opt["var_rev"]: mut_range = (WTlen-int(end), WTlen-int(start)+1)  # reverse position if var_rev added +1
    else: mut_range = (int(start)-1, int(end)) # removed -1
    mut_dict[(name, strain)] = mut_range
    return RNA, strain, mut_dict

def revComp(RNA, **opt):
    ## complement dictionary, or transform DNA to RNA
    RNA = RNA.upper()
    D2Rc = {"A":"U","T":"A","U":"A","C":"G","G":"C","R":"Y","Y":"R","M":"K",\
            "K":"M","S":"W","W":"S","B":"V","V":"B","D":"H","H":"D","N":"N"}
    if opt["var_cmp"]: RNA = "".join(D2Rc[i] for i in RNA)
    else:              RNA = RNA.replace("T","U")
    if opt["var_rev"]: RNA = RNA[::-1]
    return RNA

def defineCombinations(data_dict, **opt):
    ## create WT / mutant combinations
    ## comb_list.append(((name,strain),(name,strain)))
    comb_list = list()
    c_list = combinations(data_dict.keys(),2) if not opt["var_tra"] else combinations_with_replacement(data_dict.keys(),2)
    for nA,nB in c_list:
        if not opt["var_tra"] and nA[0] == nB[0]: continue
        la, lb = len(data_dict[nA]), len(data_dict[nB])
        if la < lb:   comb_list.append((nA,nB))
        elif la > lb: comb_list.append((nB,nA))
        else:         comb_list.append(tuple(sorted((nA,nB), reverse=False)))
    return comb_list

def loadMatrix(new_comb_list, **opt):
    ## load matrices from file and compare with fasta
    heat_dict, splash_dict = loadData(**opt)
    old_comb_list = heat_dict.keys()
    comb_list = list(set(new_comb_list).difference(set(old_comb_list)))
    return heat_dict, comb_list, splash_dict

def loadData(**opt):
    ## load data with pickle
    file = os.path.join(opt["var_pfx"], f"matrices.pcl")
    with open(file, "r+b") as fin:
        heat_dict, splash_dict = pickle.load(fin)
    return heat_dict, splash_dict

def saveMatrix(heat_dict, splash_dict, **opt):
    ## save data with pickle
    new_heats = dict()
    for ((nA,sA),(nB,sB)),lmat in heat_dict.items():
        if (len(sA.split('-')) > 1 or len(sB.split('-')) > 1) and not opt["var_svo"]: continue
        new_heats[((nA,sA),(nB,sB))] = lmat
    file = os.path.join(opt["var_pfx"], f"matrices.pcl")
    with open(file, "w+b") as pdat:
        pickle.dump((new_heats, splash_dict), pdat , protocol=4)

def createSnippets(data_dict, **opt):
    ## create sequence snippets
    t, l = opt["var_stp"], opt["var_slc"]
    snip_dict = dict()
    for name,RNA in data_dict.items():
        r = len(RNA)
        snip_dict[name] = [(i,i+l) for i in range(0,r-l+1,t)] + bool((r-l)%t) * [(r-l+t-(r-l)%t,r)]
    return snip_dict

def combineSnippets(comb_list, data_dict, snip_dict, mut_dict, **opt):
    ## combine snippets
    fold_list, l = list(), opt["var_slc"]
    for nA,nB in comb_list:
        a_snips, b_snips = snip_dict[nA], snip_dict[nB]
        if nA[1] != "WT":
            m1, m2 = mut_dict[nA]
            #a_snips = [(ai,aj) for ai,aj in a_snips if aj > m1-l and ai <= m2+l] old
            a_snips = [(ai,aj) for ai,aj in a_snips if aj >= m1 and ai < m2]
        if nB[1] != "WT":
            m1, m2 = mut_dict[nB]
            #b_snips = [(bi,bj) for bi,bj in b_snips if bj > m1-l and bi <= m2+l] old
            b_snips = [(bi,bj) for bi,bj in b_snips if bj >= m1 and bi < m2]
        if nA[1] != "WT" and nB[1] != "WT":
            tuplesA = [(nA,nB,a,b) for a,b in product(a_snips, snip_dict[nB])]
            tuplesB = [(nA,nB,a,b) for a,b in product(snip_dict[nA], b_snips)]
            fold_list += list(set(tuplesA).union(set(tuplesB)))
        else:
            fold_list += [(nA,nB,a,b) for a,b in product(a_snips, b_snips)]
    return fold_list

def foldCombinations(fold_list, data_dict, **opt):
    ## fold combinations
    res_map, res_list = list(), mp.Manager().list()
    thr, total = opt["var_thr"], len(fold_list)
    split = total//thr
    if split:
        in_split = [(s,s+split) for s in range(0,total,split)]
        in_split[-1] = (in_split[-1][0],total)
    else:
        in_split = [(0,total)]
    ## start multiprocessing
    for run,(i,j) in enumerate(in_split,1):
        p = mp.Process(target=windowFolding, args=(fold_list[i:j], data_dict, run, split, res_list, opt))
        res_map.append(p)
        p.start()
    for res in res_map:
        res.join()
    return res_list

def windowFolding(fold_list, data_dict, run, split, res_list, opt):
    ## create matrices for each combination 
    last_comb = ("","")
    for it,(nA, nB, (ai,aj), (bi,bj)) in enumerate(fold_list):
        if int(it*10000/split) % 10 == 0 and int(it*100/split) != 0:
            print(f"Status: Folding {run:>2d} ... {it*100/split:>4.1f} %             ", end="\r")
        ## do cofold
        if last_comb != (nA, nB):
            if last_comb != ("",""):
                res_list.append((last_comb[0], last_comb[1], lmat, cmat))
            last_comb = (nA,nB)
            RNAa, RNAb = data_dict[nA], data_dict[nB]
            m, n = len(RNAa), len(RNAb)
            lmat = zeros(shape=(m,n))
            cmat = zeros(shape=(m,n))
        RNA = f"{RNAa[ai:aj]}&{RNAb[bi:bj]}"
        constraint = f"{'<'*(aj-ai)}{'>'*(bj-bi)}"
        mfe, pattern = doCofold(RNA, constraint, **opt)
        for i,j in zip(range(ai,aj),range(bi,bj)[::-1]):
            lmat[i,j] += mfe
            cmat[i,j] += 1
    res_list.append((last_comb[0], last_comb[1], lmat, cmat))

def doCofold(RNA, constraint, **opt):
    ## do Cofold
    cvar.dangles = opt["var_dng"]
    cvar.noLP = int(opt["var_nlp"])
    cvar.temperature = opt["var_tmp"]
    fc = fold_compound(RNA)
    fc.constraints_add(constraint, CONSTRAINT_DB | CONSTRAINT_DB_DEFAULT)
    pattern, mfe = fc.mfe_dimer()
    return mfe, pattern




################################################################################
## handle big genomes
################################################################################

def foldCombinationsBig(fold_list, data_dict, **opt):
    ## fold combinations
    res_map, res_list = list(), mp.Manager().list()
    thr, total = opt["var_thr"], len(fold_list)
    split = total//thr
    if split:
        in_split = [(s,s+split) for s in range(0,total,split)]
        in_split[-1] = (in_split[-1][0],total)
    else:
        in_split = [(0,total)]
    ## start multiprocessing
    for run,(i,j) in enumerate(in_split,1):
        p = mp.Process(target=windowFoldingBig, args=(fold_list[i:j], data_dict, run, split, res_list, opt))
        res_map.append(p)
        p.start()
    for res in res_map:
        res.join()
    if opt["var_big"]:
        res_list = loadResList(res_list, **opt)
    return res_list

def windowFoldingBig(fold_list, data_dict, run, split, res_list, opt):
    ## create matrices for each combination 
    last_comb = ("","")
    for it,(nA, nB, (ai,aj), (bi,bj)) in enumerate(fold_list):
        if int(it*10000/split) % 10 == 0 and int(it*100/split) != 0:
            print(f"Status: Folding {run:>2d} ... {it*100/split:>4.1f} %             ", end="\r")
        ## do cofold
        if last_comb != (nA, nB):
            if last_comb != ("",""):
                if opt["var_big"]:
                    saveResList(last_comb[0], last_comb[1], run, lmat, "lmat", it, opt)
                    saveResList(last_comb[0], last_comb[1], run, cmat, "cmat", it, opt)
                    res_list.append((last_comb[0], last_comb[1],run,it))
                else:
                    res_list.append((last_comb[0], last_comb[1], lmat, cmat))
            last_comb = (nA,nB)
            RNAa, RNAb = data_dict[nA], data_dict[nB]
            m, n = len(RNAa), len(RNAb)
            lmat = zeros(shape=(m,n))
            cmat = zeros(shape=(m,n))
        RNA = f"{RNAa[ai:aj]}&{RNAb[bi:bj]}"
        constraint = f"{'<'*(aj-ai)}{'>'*(bj-bi)}"
        mfe, pattern = doCofold(RNA, constraint, **opt)
        for i,j in zip(range(ai,aj),range(bi,bj)[::-1]):
            lmat[i,j] += mfe
            cmat[i,j] += 1
    if opt["var_big"]:
        saveResList(last_comb[0], last_comb[1], run, lmat, "lmat", it, opt)
        saveResList(last_comb[0], last_comb[1], run, cmat, "cmat", it, opt)
        res_list.append((last_comb[0], last_comb[1],run,it))
    else:
        res_list.append((last_comb[0], last_comb[1], lmat, cmat))

def saveResList(nA, nB, run, mat, mtype, it, opt):
    ## saves matrix to file
    res_file = os.path.join(opt["var_pfx"], "vRNAsite_temp", f"{run}_{nA}_{nB}_{mtype}_{it}.pcl")
    with open(res_file, "w+b") as pout:
        pickle.dump(mat, pout, protocol=4)

def loadResList(res_files, **opt):
    ## load and combine matrices
    res_list, total = list(), len(res_files)
    last_comb = ("","")
    lmat, cmat = zeros(shape=(1,1)), zeros(shape=(1,1))
    for f,(nA,nB,run,it) in enumerate(sorted(res_files)):
        print(f"Status: Reading and combining {nA}x{nB} file number {f} of {total} ...                   ", end="\r")
        lres_file = os.path.join(opt["var_pfx"], "vRNAsite_temp", f"{run}_{nA}_{nB}_lmat_{it}.pcl")
        with open(lres_file, "r+b") as lin:
            nlmat = pickle.load(lin)
        cres_file = os.path.join(opt["var_pfx"], "vRNAsite_temp", f"{run}_{nA}_{nB}_cmat_{it}.pcl")
        with open(cres_file, "r+b") as cin:
            ncmat = pickle.load(cin)
        os.remove(lres_file)
        os.remove(cres_file)
        if last_comb != (nA, nB):
            if last_comb != ("",""):
                res_list.append((last_comb[0], last_comb[1], lmat, cmat))
            last_comb = (nA,nB)
            lmat, cmat = nlmat, ncmat
        else:
            lmat = nlmat + lmat
            cmat = ncmat + cmat
    res_list.append((last_comb[0], last_comb[1], lmat, cmat))
    return res_list




################################################################################
## handle matrices
################################################################################

def createMatrices(res_list, heat_dict, mut_dict, comb_list, **opt):
    ## combine matrices
    mat_dict = combineMatrices(res_list, **opt)
    ## average matrices
    heat_dict = averageMatrices(mat_dict, heat_dict, **opt)
    ## create mutant matrices
    heat_dict = mutateMatrices(heat_dict, mut_dict, comb_list, **opt)
    ## plot matrices
    return heat_dict

def combineMatrices(res_list, **opt):
    ## combine matrices
    mat_dict, last_comb, total = dict(), ("",""), len(res_list)
    for it,(nA, nB, lmat, cmat) in enumerate(res_list):
        if int(it*1000/total) % 1 == 0 and int(it*100/total) != 0:
            print(f"Status: Sorting ... {it*100/total:>4.1f} %             ", end="\r")
        if last_comb != (nA, nB):
            if last_comb != ("", ""):
                mat_dict[last_comb] = (nlmat, ncmat)
            last_comb = (nA, nB)
            nlmat, ncmat = mat_dict.get(last_comb, (zeros(shape=lmat.shape), zeros(shape=lmat.shape)))
        nlmat = nlmat + lmat
        ncmat = ncmat + cmat
    mat_dict[last_comb] = (nlmat, ncmat)
    return mat_dict

def averageMatrices(mat_dict, heat_dict, **opt):
    ## average matrices
    for (nA, nB),matrices in mat_dict.items():
        lmat, cmat = matrices
        cmat = cmat + (cmat == 0)
        lmat = lmat / cmat
        heat_dict[(nA, nB)] = lmat
    return heat_dict

def mutateMatrices(heat_dict, mut_dict, comb_list, **opt):
    ## create mutant matrices
    mut_list = [(nA, nB) for nA, nB in heat_dict.keys() if (nA[1] != "WT" or nB[1] != "WT") and (nA,nB) in comb_list]
    l = opt["var_slc"]
    for (nA,sA),(nB,sB) in mut_list:
        WT_mat, MT_mat = heat_dict[((nA,"WT"),(nB,"WT"))], heat_dict[((nA,sA),(nB,sB))]
        MW_mat, mask = zeros(shape=WT_mat.shape), zeros(shape=WT_mat.shape) != 0
        if sA != "WT":
            m1, m2 = mut_dict[(nA,sA)]
            mask[range(max(m1-l,0),min(m2+l,mask.shape[0]-1)+1),:] = True
        if sB != "WT":
            m1, m2 = mut_dict[(nB,sB)]
            mask[:,range(max(m1-l,0),min(m2+l,mask.shape[1]-1)+1)] = True
        MW_mat[mask] = MT_mat[mask]
        MW_mat[mask == False] = WT_mat[mask == False]
        heat_dict[((nA,sA),(nB,sB))] = MW_mat
    return heat_dict

def orientateMatrices(heat_dict, comb_list, **opt):
    ## create IAV orientation matrices
    var_pfx = opt["var_pfx"]
    opt["var_pfx"] = os.path.join(var_pfx, "weight_matrices")
    makeDir(**opt)
    orient_dict = dict() ## add pool
    with open("/home/vo54saz/projects/test_vsite.tsv", "w") as testx:
        testx.write(f"nA\tsA\tnB\tsB\tA_n\tA_m\tUU_n\tUU_m\tRR_n\tRR_m\tUR_n\tUR_m\tRU_n\tRU_m\n")
        for ((nA,sA),(nB,sB)),lmat in heat_dict.items():
            #testx.write(f"{sA}\t{len(sA.split('-'))}\t{sB}\t{len(sA.split('-'))}\n")
            if len(sA.split('-')) > 1 or len(sB.split('-')) > 1: continue
            wg_name = (((nA,f"{sA.split('-')[0]}-U"), (nB,f"{sB.split('-')[0]}-U")),
                       ((nA,f"{sA.split('-')[0]}-R"), (nB,f"{sB.split('-')[0]}-R")),
                       ((nA,f"{sA.split('-')[0]}-U"), (nB,f"{sB.split('-')[0]}-R")),
                       ((nA,f"{sA.split('-')[0]}-R"), (nB,f"{sB.split('-')[0]}-U")))
            wg_mats = createWeight(lmat, nA, sA, nB, sB, **opt)
            orient_dict[((nA,sA),(nB,sB))] = lmat
            n, m = lmat.shape
            testx.write(f"{nA}\t{sA}\t{nB}\t{sB}\t{n}\t{m}")
            for name,wgmat in zip(wg_name, wg_mats):
                n, m = wgmat.shape
                testx.write(f"\t{n}\t{m}")
                orient_dict[name] = lmat * wgmat
            testx.write(f"\n")
    opt["var_pfx"] = var_pfx
    return orient_dict

def createWeight(lmat, nA, sA, nB, sB, **opt):
    ################################################################################
    ## create weight matrix
    n, m = lmat.shape
    wgmat = zeros(shape=(m,n))
    m05, n05 = int(ceil(m/2)),  int(ceil(n/2))
    m04, n04 = int(floor(m/2)), int(floor(n/2))
    if opt["var_dst"] >= 1:
        nd = max(int(n05 - opt["var_dst"]),0)
        if opt["var_dsc"] and opt["var_dsc"] >= 1:
            nc = max(int(n05 - opt["var_dst"] - opt["var_dsc"]),0)
        else:
            nc = 0
    else:
        nd = int(ceil(n05 * (1.0 - opt["var_dst"])))
        if opt["var_dsc"] and opt["var_dsc"] < 1:
            nc = int(ceil(n05 * (1.0 - opt["var_dst"] - opt["var_dsc"])))
        else:
            nc = 0
    warray = [0.0] * nc + list(linspace(0,1,nd-nc)) + [1.0] * (n05-nd)
    for x,wgt in enumerate(warray):
        for i,j in enumerate(range(n05-x,n05+1)):
            xl = [i,   j-1, i,   j-1,   m-i-1, m-j, m-i-1, m-j  ]
            yl = [j-1, i,   n-j, n-i-1, j-1,   i,   n-j,   n-i-1]
            wgmat[xl,yl] = wgt
    for x,wgt in enumerate(warray[::-1]):
        for i,j in enumerate(range(n05-x,n05+1)):
            if n05+i < m05:
                wgmat[[n05+i, n05+i], [j-1, n-j]] = wgt
            if m-1-i-n05 >= m05:
                wgmat[[m-1-i-n05, m-1-i-n05], [n-j, j-1]] = wgt
    abcd = {"A":(0,n05,0,m05), "B":(0,n05,m04,m), "C":(n04,n,0,m05), "D":(n04,n,m04,m)}
    ABCD = {k:array([w[a:b] for w in wgmat[c:d]]) for k,(a,b,c,d) in abcd.items()}
    abcd = [("A","B"), ("B","A"), ("C","D"), ("D","C")]
    ABCD = {f"{i}{j}":concatenate((ABCD[i],ABCD[j][m05-m04:]), axis=0) for i,j in abcd}
    abcd = [("AB","CD"), ("DC","BA"), ("BA","DC"), ("CD","AB")]
    ABCD = [transpose(concatenate((ABCD[i],array([w[n05-n04:] for w in ABCD[j]])), axis=1)) for i,j in abcd]
    if opt["var_plh"]:
        plotMatricesMulti(ABCD[0], os.path.join(opt["var_pfx"], f"wgm_{nA}-{sA}_{nB}-{sB}_1"), f"{nA}-{sA}", f"{nB}-{sB}", ABCD[0].max(), opt) #UU
        plotMatricesMulti(ABCD[1], os.path.join(opt["var_pfx"], f"wgm_{nA}-{sA}_{nB}-{sB}_2"), f"{nA}-{sA}", f"{nB}-{sB}", ABCD[1].max(), opt) #RR
        plotMatricesMulti(ABCD[2], os.path.join(opt["var_pfx"], f"wgm_{nA}-{sA}_{nB}-{sB}_3"), f"{nA}-{sA}", f"{nB}-{sB}", ABCD[2].max(), opt) #sUlR
        plotMatricesMulti(ABCD[3], os.path.join(opt["var_pfx"], f"wgm_{nA}-{sA}_{nB}-{sB}_4"), f"{nA}-{sA}", f"{nB}-{sB}", ABCD[3].max(), opt) #sRlU
    return ABCD

#from numpy import concatenate, zeros, linspace, transpose, array, savetxt
#from math import ceil, floor
def weightmatrix(m, n, d):
    from numpy import concatenate, zeros, linspace, transpose, array, savetxt
    from math import ceil, floor
    wgmat = zeros(shape=(m,n))
    m05, n05 = int(ceil(m/2)),  int(ceil(n/2))
    m04, n04 = int(floor(m/2)), int(floor(n/2))
    nd = int(ceil(n05 * (1.0 - d)))
    warray = list(linspace(0,1,nd)) + [1.0] * (n05-nd)
    for x,wgt in enumerate(warray):
        for i,j in enumerate(range(n05-x,n05+1)):
            xl = [i,   j-1, i,   j-1,   m-i-1, m-j, m-i-1, m-j  ]
            yl = [j-1, i,   n-j, n-i-1, j-1,   i,   n-j,   n-i-1]
            wgmat[xl,yl] = wgt
    for x,wgt in enumerate(warray[::-1]):
        for i,j in enumerate(range(n05-x,n05+1)):
            if n05+i < m05:
                wgmat[[n05+i, n05+i], [j-1, n-j]] = wgt
            if m-1-i-n05 >= m05:
                wgmat[[m-1-i-n05, m-1-i-n05], [n-j, j-1]] = wgt
    abcd = {"A":(0,n05,0,m05), "B":(0,n05,m04,m), "C":(n04,n,0,m05), "D":(n04,n,m04,m)}
    ABCD = {k:array([w[a:b] for w in wgmat[c:d]]) for k,(a,b,c,d) in abcd.items()}
    abcd = [("A","B"), ("B","A"), ("C","D"), ("D","C")]
    ABCD = {f"{i}{j}":concatenate((ABCD[i],ABCD[j][m05-m04:]), axis=0) for i,j in abcd}
    abcd = [("AB","CD"), ("DC","BA"), ("BA","DC"), ("CD","AB")]
    ABCD = [transpose(concatenate((ABCD[i],array([w[n05-n04:] for w in ABCD[j]])), axis=1)) for i,j in abcd]
    savetxt("weight_matrix_1.csv",ABCD[0],delimiter=",")
    savetxt("weight_matrix_2.csv",ABCD[1],delimiter=",")
    savetxt("weight_matrix_3.csv",ABCD[2],delimiter=",")
    savetxt("weight_matrix_4.csv",ABCD[3],delimiter=",")
    return ABCD




################################################################################
## extract interactions
################################################################################

class links(object):
    def __init__(self, aSeq, aType, bSeq, bType, ai, aj, bi, bj, peak, mfe, structure, RNA, a_mfe, a_structure, b_mfe, b_structure, free_mfe, free_structure, alen, blen, distance):
        self.aSeq           = aSeq
        self.aType          = aType
        self.ai             = int(ai)
        self.aj             = int(aj)
        self.bSeq           = bSeq
        self.bType          = bType
        self.bi             = int(bi)
        self.bj             = int(bj)
        self.peak           = round(peak,2)
        self.distance       = int(distance)
        self.sp_mean        = float('nan')
        self.sp_min         = float('nan')
        self.sp_max         = float('nan')
        self.sp_median      = float('nan')
        self.mfe            = round(mfe,2)
        self.RNA            = RNA
        self.structure      = structure
        self.a_mfe          = round(a_mfe,2)
        self.a_structure    = a_structure
        self.a_mean         = float('nan')
        self.a_min          = float('nan')
        self.a_max          = float('nan')
        self.a_median       = float('nan')
        self.b_mfe          = round(b_mfe,2)
        self.b_structure    = b_structure
        self.b_mean         = float('nan')
        self.b_min          = float('nan')
        self.b_max          = float('nan')
        self.b_median       = float('nan')
        self.free_mfe       = round(free_mfe,2)
        self.free_structure = free_structure
        self.alen           = int(alen)
        self.blen           = int(blen)
    def plot(self, sep, extended=False):
        bio = ["sp_mean","sp_min","sp_max","sp_median","a_mean","a_min","a_max","a_median","b_mean","b_min","b_max","b_median",]
        ldat = sep.join([f"{var}" for key,var in vars(self).items() if not (not extended and key in bio)])
        return ldat

def extractCandidates(heat_dict, data_dict, cand_list, **opt):
    ## create matrix and extract candidates
    if opt["var_plh"]:
        var_pfx = opt["var_pfx"]
        opt["var_pfx"] = os.path.join(var_pfx, "cand_plots")
        makeDir(**opt)
    c_list, mplot_list, cnd_list = list(), list(), list(set([((lk.aSeq,lk.aType),(lk.bSeq,lk.bType)) for lk in cand_list]))
    pool_list = [(nA, sA, data_dict[(nA,sA.split("-")[0])],
                  nB, sB, data_dict[(nB,sB.split("-")[0])],
                  cbmat, opt) 
                  for ((nA,sA),(nB,sB)),cbmat in heat_dict.items()
                  if ((nA,sA),(nB,sB)) not in cnd_list 
                      and not (not opt["var_tra"] and (nA,sA) == (nB,sB))]
    if pool_list:
        if opt["var_cdp"]:
            with mp.Pool(processes=opt["var_thr"]) as p:
                p_list = p.starmap(extractCandidatesMulti, pool_list)
            for clist,mplot in p_list:
                mplot_list.append(mplot)
                for c in clist:
                    c_list.append(c)
        else:
            for tup in pool_list:
                clist,mplot = extractCandidatesMulti(*tup)
                mplot_list.append(mplot)
                for c in clist:
                    c_list.append(c)
        cand_list.extend(c_list)
    ## plot candidates
    if opt["var_plh"] and not opt["var_plp"]:
        for mplot in mplot_list:
            plotMatricesMulti(*mplot)
    ## remove duplicates
    if opt["var_tra"]:
        new_cand, cont_list = list(), list()
        for lk in cand_list:
            if (lk.aSeq,lk.aType) == (lk.bSeq,lk.bType):
                if (lk.bSeq,lk.bType,lk.aSeq,lk.aType,lk.bi,lk.bj,lk.ai,lk.aj) not in cont_list:
                    new_cand.append(lk)
                    cont_list.append((lk.aSeq,lk.aType,lk.bSeq,lk.bType,lk.ai,lk.aj,lk.bi,lk.bj))
            else:
                new_cand.append(lk)
        cand_list = new_cand
    if opt["var_plh"]: opt["var_pfx"] = var_pfx
    return cand_list

def extractCandidatesMulti(nA, sA, RNAa, nB, sB, RNAb, cbmat, opt):
    ## create matrix and extract candidates
    c_list = list()
    m, n = cbmat.shape
    cen_df = [cen for cblist in cbmat for cen in cblist]
    RNAa_df = [i for i in range(m) for j in range(n)]
    RNAb_df = [j for j in list(range(n))*m]
    df = [(cen,ai,bj) for cen,ai,bj in zip(cen_df,RNAa_df,RNAb_df) if cen <= opt["var_clb"]]
    clusts, cmat = makeClusters(sorted(df), cbmat, **opt)
    plm = opt["var_plm"]
    for it,(cen,ai,aj,bi,bj) in enumerate(clusts):
        aj, bj = aj+1, bj+1 # added +1
        RNA = f"{RNAa[ai:aj]}&{RNAb[bi:bj]}"
        aRNA, bRNA = RNA.split("&")
        if len(aRNA) < opt["var_cdm"] or len(bRNA) < opt["var_cdm"]: continue
        fmfe, fpat = doCofold(RNA, f"{'.'*(aj-ai)}{'.'*(bj-bi)}", **opt)
        fpat = fpat[:len(RNA.split("&")[0])]+"&"+fpat[len(RNA.split("&")[0]):]
        RNA, mfe, pattern, ai, aj, bi, bj = trimFold(RNA, ai, aj, bi, bj, nA, nB, sA, sB, **opt)
        if mfe > opt["var_cdc"]: continue
        aRNA, bRNA = RNA.split("&")
        if len(aRNA) < opt["var_cdm"] or len(bRNA) < opt["var_cdm"]: continue
        amfe, apat = doCofold(aRNA, len(aRNA)*".", **opt)
        bmfe, bpat = doCofold(bRNA, len(bRNA)*".", **opt)
        if bi >= ai: distance = bi - aj + 1
        else: distance = ai - bj # removed +1
        lk = links(nA, sA, nB, sB, ai, aj, bi, bj, cen, mfe, pattern, RNA, amfe, apat, bmfe, bpat, fmfe, fpat, m, n, distance)
        c_list.append(lk)
    if opt["var_plh"] and opt["var_plp"]:
        plotMatricesMulti(cmat, os.path.join(opt["var_pfx"], f"{nA}-{sA}_{nB}-{sB}"), f"{nA} {sA}", f"{nB} {sB}", plm, opt)
    return c_list, (cmat, os.path.join(opt["var_pfx"], f"{nA}-{sA}_{nB}-{sB}"), f"{nA} {sA}", f"{nB} {sB}", plm, opt)

def makeClusters(df, cbmat, **opt):
    ## cluster and split mfe points according to mfe
    cdf_list, total = list(), len(df)
    for it,(cen,ai,bj) in enumerate(df):
        print(f"Status: Clustering ... {it*100/total:>4.1f} %             ", end="\r")
        if cbmat[ai,bj] == 0: continue
        if cen > opt["var_clp"]: break
        cbmat, cluster = getNeighbors(cbmat, cen, ai, bj, **opt)
        cdf_list.append(cluster)
    return cdf_list, cbmat

def getNeighbors(cbmat, cen, ai, bj, **opt):
    ## return neighbors from i
    m, n = cbmat.shape
    cbmat[ai,bj] = opt["var_clv"]
    m_l, x_l, y_l, last_m = [cen], [ai], [bj], m
    queue = list(product([-1,0,+1],[-1,0,+1]))
    queue.remove((0,0))
    queue = [(cen,i,j) for i,j in queue]
    while queue:
        mq,i,j = queue.pop()
        if ai+i >= m or bj+j >= n or ai+i < 0 or bj+j < 0: continue
        mx = cbmat[ai+i,bj+j]
        if mx < mq or mx == 0 or mx > opt["var_clb"]: continue
        m_l.append(mx)
        x_l.append(ai+i)
        y_l.append(bj+j)
        cbmat[ai+i,bj+j] = opt["var_clv"]
        queue.extend([(mx, i+1, j+1), (mx, i-1, j-1), (mx, i+1, j-1), (mx, i-1, j+1),
                      (mx, i+1, j  ), (mx, i-1, j  ), (mx, i  , j+1), (mx, i  , j-1)])
    return cbmat, (min(m_l), min(x_l), max(x_l), min(y_l), max(y_l))

def trimFold(RNA, ai, aj, bi, bj, nA, nB, sA, sB, **opt):
    ## trim RNA fold pattern
    mfe, pattern = doCofold(RNA, f"{'<'*(aj-ai)}{'>'*(bj-bi)}", **opt)
    r1, r2 = RNA.split("&")
    l1, l2 = len(r1), len(r2)
    p1, p2 = pattern[:l1], pattern[l1:]
    sl1, sl2 = l1-len(p1.lstrip(".")), l2-len(p2.lstrip("."))
    sr1, sr2 = l1-len(p1.rstrip(".")), l2-len(p2.rstrip("."))
    if sl1 > opt["var_cdd"]: sl1 -= opt["var_cdd"]
    else: sl1 = 0
    if sl2 > opt["var_cdd"]: sl2 -= opt["var_cdd"]
    else: sl2 = 0
    if sr1 > opt["var_cdd"]: sr1 -= opt["var_cdd"]
    else: sr1 = 0
    if sr2 > opt["var_cdd"]: sr2 -= opt["var_cdd"]
    else: sr2 = 0
    ai, aj, bi, bj = ai+sl1, aj-sr1, bi+sl2, bj-sr2
    RNA = f"{r1[sl1:l1-sr1]}&{r2[sl2:l2-sr2]}"
    if len(RNA) > 1:
        mfe, pattern = doCofold(RNA, f"{'<'*(l1-sl1-sr1)}{'>'*(l2-sl2-sr2)}", **opt)
        pattern = pattern[:len(RNA.split("&")[0])]+"&"+pattern[len(RNA.split("&")[0]):]
    else:
        RNA, mfe, pattern = "", 0.0, ""
    return RNA, mfe, pattern, ai, aj, bi, bj

def intraStructure(RNA, pattern):
    ## searches for intra structure folds
    pattern = pattern[:len(RNA.split("&")[0])]+"&"+pattern[len(RNA.split("&")[0]):]
    aSeq = True
    for i in pattern:
        if   i == "(" and     aSeq: pass
        elif i == "(" and not aSeq: return True
        elif i == ")" and     aSeq: return True
        elif i == ")" and not aSeq: pass
        elif i == "&": aSeq = False
    return False

def saveCandidates(cand_list, **opt):
    ## save all candidates to file
    with open(os.path.join(opt["var_pfx"], f"candidates.pcl"), "w+b") as pin:
        pickle.dump(cand_list, pin , protocol=4)
    ext, end, ehd = "", "", ""
    if opt["var_nex"] == "peak": ext += f"_peak{opt['var_clp']}"
    if opt["var_nex"] == "mfe":  ext += f"_mfe{opt['var_cdc']}"
    if opt["var_tra"]:           ext += "_intra"
    with open(os.path.join(opt["var_pfx"], f"{os.path.basename(opt['var_pfx'])}{ext}.tsv"), "w") as wrt:
        for i,lk in enumerate(cand_list):
            if i == 0:
                bio = ["sp_mean","sp_min","sp_max","sp_median","a_mean","a_min","a_max","a_median","b_mean","b_min","b_max","b_median",]
                #ldat = sep.join([f"{var}" for key,var in vars(self).items() if not (not extended and key in bio)])
                if opt["var_spd"]:
                    wrt.write("\t".join(lk.__dict__.keys())+"\n")
                else:
                    wrt.write("\t".join([k for k in lk.__dict__.keys() if k not in bio])+"\n")
            if opt["var_rvp"]:
                lk.ai, lk.aj, lk.bi, lk.bj = lk.alen-lk.ai, lk.alen-lk.aj+1, lk.blen-lk.bi, lk.blen-lk.bj+1 # removed -1
            else:
                lk.ai, lk.aj, lk.bi, lk.bj = lk.ai+1, lk.aj, lk.bi+1, lk.bj
            wrt.write(lk.plot("\t", opt["var_spd"])+"\n")

def loadCandidates(**opt):
    with open(os.path.join(opt["var_pfx"], f"candidates.pcl"), "r+b") as pout:
        cand_list = pickle.load(pout)
    return cand_list




################################################################################
## handle SPLASH data
################################################################################

def loadSPLASH(comb_list, data_dict, **opt):
    ## load SPLASH data and transform to vRNAsite matrix
    if not os.path.isdir(opt["var_spd"]):
        print(f"Error: Not a valid directory!")
        sys.exit()
    splash_files = os.listdir(opt["var_spd"])
    total = len(splash_files)
    splash_dict = dict()
    for it,sfile in enumerate(splash_files):
        # initiate SPLASH matrix
        print(f"Status: Reading SPLASH matrix {it+1} of {total} ...             ", end="\r")
        An, As, Bn, Bs = os.path.splitext(sfile)[0].split("-")
        nA, nB = (An,As), (Bn,Bs)
        RNAa = data_dict.get(nA, "")
        if not RNAa:
            print(f"Warning: No fasta entry for \">{An} {As}\"!")
            continue
        RNAb = data_dict.get(nB, "")
        if not RNAb:
            print(f"Warning: No fasta entry for \">{Bn} {Bs}\"!")
            continue
        m, n = len(RNAa), len(RNAb)
        smat = zeros(shape=(m,n))
        with open(os.path.join(opt["var_spd"],sfile), "r") as sin:
            for line in sin:
                i,j,c = line.strip().split()
                c = float(c)
                if opt["var_rev"]: i,j = m-int(i)-1, n-int(j)-1
                else:              i,j = int(i), int(j)
                smat[i,j] += c
        # intra interactions
        if opt["var_tra"] and nA == nB:
            tmat = smat.transpose()
            if tmat.all() != smat.all():
                smat = smat + tmat
        elif not opt["var_tra"] and  nA == nB:
            continue
        # save interactions
        if (nA,nB) in comb_list:
            xmat = splash_dict.get((nA,nB), zeros(shape=(m,n)))
            xmat = xmat + smat
            splash_dict[(nA,nB)] = xmat
        elif (nB,nA) in comb_list:
            xmat = splash_dict.get((nB,nA), zeros(shape=(n,m)))
            xmat = xmat + smat.transpose()
            splash_dict[(nB,nA)] = xmat
    return splash_dict

def addSplash(cand_list, splash_dict, **opt):
    ## adds SPLASH data to candidate list
    new_cand = list()
    for lk in cand_list:
        #print(nA,sA,nB,sB) = SC35M_HA mutS4 SC35M_PB2 WT
        smat = splash_dict.get(((lk.aSeq,lk.aType.split("-")[0]),(lk.bSeq,lk.bType.split("-")[0])), zeros(shape=(1,1)))
        amat = splash_dict.get(((lk.aSeq,lk.aType.split("-")[0]),(lk.aSeq,lk.aType.split("-")[0])), zeros(shape=(1,1)))
        bmat = splash_dict.get(((lk.bSeq,lk.bType.split("-")[0]),(lk.bSeq,lk.bType.split("-")[0])), zeros(shape=(1,1)))
        data = ["sp","a","b"]
        pattern = [lk.structure,lk.a_structure,lk.b_structure]
        mat = [smat,amat,bmat]
        apos = [lk.ai,lk.ai,lk.bi]
        bpos = [lk.bi,lk.ai,lk.bi]
        for t,pat,xmat,ai,bi in zip(data, pattern, mat, apos, bpos):
            if xmat.any():
                s_list = list()
                for i,j in getBasePairs(pat, ai, bi, **opt):
                    sdat = xmat[i,j]
                    if not isnan(sdat):
                        s_list.append(sdat)
                if s_list:
                    lk.__dict__[f"{t}_min"] = round(min(s_list),3)
                    lk.__dict__[f"{t}_max"] = round(max(s_list),3)
                    lk.__dict__[f"{t}_mean"] = round(mean(s_list),3)
                    lk.__dict__[f"{t}_median"] = round(median(s_list),3)
        if opt["var_spr"] > lk.sp_mean: continue
        new_cand.append(lk)
    return new_cand

def getBasePairs(pat, ai, bi, **opt):
    ## determine base pairing positions of the structure
    st, bp, k, s = list(), list(), ai, 0
    for l,x in enumerate(pat):
        if x == "&":
            k = bi
            s = l+1
        elif x == "(":
            st.append(k+l-s)
        elif x == ")":
            i, j = st.pop(), k+l-s
            bp.append((i,j))
    return bp




################################################################################
## plot functions
################################################################################

def plotMatrices(heat_dict, name, **opt):
    ## plot all matrices
    var_pfx = opt["var_pfx"]
    opt["var_pfx"] = os.path.join(var_pfx, name)
    makeDir(**opt)
    plm = opt["var_plm"]
    print_list = [(cbmat, os.path.join(opt["var_pfx"], f"{nA}-{sA}_{nB}-{sB}"),
                   f"{nA} {sA}", f"{nB} {sB}", plm, opt) 
                   for ((nA,sA),(nB,sB)),cbmat in heat_dict.items()]
    if opt["var_plp"]:
        with mp.Pool(processes=opt["var_thr"]) as p:
            p.starmap(plotMatricesMulti, print_list)
    else:
        for (mat, f_name, xname, yname, plm, opt) in print_list:
            plotMatricesMulti(mat, f_name, xname, yname, plm, opt)
    opt["var_pfx"] = var_pfx

def plotMatricesMulti(mat, f_name, xname, yname, plm, opt):
    ## pyplot matrix plot
    xmat = nan_to_num(mat, copy=True)
    if not opt["var_pln"]:
        mat = xmat
    if not plm:
        plm = xmat.max()
    s = opt["var_pls"]
    n, m = mat.shape # default: 6.4, 4.8 -> 1.6 for legend
    mx = (4.8 / n) * m + 1.6
    pz = plt.figure(figsize=(mx*s,4.8*s))
    if opt["var_rvp"]: mat = flip(mat)
    if opt["var_plc"]:
        plt.imshow(mat, cmap='viridis_r', vmin=xmat.min(), vmax=plm, interpolation="none")
    else:
        plt.imshow(mat, cmap='viridis', vmin=xmat.min(), vmax=plm, interpolation="none")
    ax = plt.gca()
    xs, xe = ax.get_xlim()
    ys, ye = ax.get_ylim()
    ax.xaxis.set_tick_params(rotation=45)
    ax.xaxis.set_ticks(arange(xs, xe, opt["var_plx"]))
    ax.yaxis.set_ticks(arange(ye, ys, opt["var_ply"]))
    plt.gcf().subplots_adjust(bottom=0.15*s)
    plt.xlabel(yname)
    plt.ylabel(xname)
    plt.colorbar()
    if opt["var_plv"]: pz.savefig(f"{f_name}.svg", bbox_inches = 'tight', pad_inches = 0.1*s)
    pz.savefig(f"{f_name}.pdf", bbox_inches = 'tight', pad_inches = 0.1*s)
    plt.close(pz)




################################################################################
## distribution plot functions
################################################################################

def createDistributionPlot(cand_list, data_dict, **opt):
    ## pyplot matrix plot
    var_pfx = opt["var_pfx"]
    opt["var_pfx"] = os.path.join(var_pfx, "distribution_plots")
    makeDir(**opt)
    ## creat complete bar plot
    peaks = [lk.peak for lk in cand_list]
    pd = getDistribution(peaks)
    f_name = os.path.join(opt["var_pfx"], f"complete-distribution")
    plotBarPlot(pd, f_name, **opt)
    ## create bar plots for each segment
    for k1 in data_dict.keys():
        for k2 in data_dict.keys():
            if k1 == k2: continue
            cand_peaks = list()
            for lk in cand_list:
                c1, c2, peak = (lk.aSeq,lk.aType), (lk.bSeq,lk.bType), lk.peak
                if (k1==c1 and k2==c2) or (k1==c2 and k2==c1):
                    cand_peaks.append(peak)
            pd = getDistribution(cand_peaks)
            f_name = os.path.join(opt["var_pfx"], f"{k1[0]}-{k1[1]}_{k2[0]}-{k2[1]}")
            plotBarPlot(pd, f_name, **opt)
    opt["var_pfx"] = var_pfx

def getDistribution(peaks):
    ## calculate distribution
    pd = {m:0 for m in arange(int(round(max(peaks),0)),int(round(min(peaks),0))-1,-1)}
    for p in peaks:
        pd[int(round(p-0.001,0))] += 1
    return pd

def plotBarPlot(pd, f_name, **opt):
    ## create bar plot
    s = opt["var_dss"]
    pz = plt.figure(figsize=(len(pd.keys())*0.7*s,4*s))
    x = arange(len(pd.keys()))
    h = [p for k,p in pd.items()]
    plt.bar(x, height=h, color="#ef6548")
    l = [f"{k}" for k in pd.keys()]
    plt.xticks(x, l)
    label = [f"{i}" for i in h]
    for i in range(len(pd.keys())):
        plt.text(x=x[i] , y=h[i], s=label[i], ha="center", va="bottom")
    plt.xlabel("average free energy in kcal/mol")
    plt.ylabel("number of interacions")
    plt.subplots_adjust(bottom=0.2, top=1.2)
    if opt["var_dsv"]: pz.savefig(f"{f_name}.svg", bbox_inches = 'tight', pad_inches = 0.1*s)
    pz.savefig(f"{f_name}.pdf", bbox_inches = 'tight', pad_inches = 0.1*s)
    plt.close(pz)




################################################################################
## bokeh function
################################################################################

def plotBokeh(heat_dict, data_dict, splash_dict, **opt):
    ## create bokeh plots
    var_pfx = opt["var_pfx"]
    opt["var_pfx"] = os.path.join(var_pfx, "bokeh_plots")
    makeDir(**opt)
    print_list = [(nA,sA,data_dict[(nA,sA.split("-")[0])],
                   nB,sB,data_dict[(nB,sB.split("-")[0])],
                   cbmat,splash_dict.get(((nA,sA),(nB,sB)),zeros(shape=(1,1))),opt) 
                   for ((nA,sA),(nB,sB)),cbmat in heat_dict.items()]
    if opt["var_bkp"]:
        with mp.Pool(processes=opt["var_thr"]) as p:
            p.starmap(plotBokehMulti, print_list)
    else:
        for (nameB,strainB,RNAb,nameA,strainA,RNAa,cbmat,spmat,opt) in print_list:
            plotBokehMulti(nameB,strainB,RNAb,nameA,strainA,RNAa,cbmat,spmat,opt)
    opt["var_pfx"] = var_pfx

def plotBokehMulti(nameB, strainB, RNAb, nameA, strainA, RNAa, cbmat, spmat, opt):
    ## create bokeh plots
    m, n = len(RNAa), len(RNAb)
    ## create data frame items
    cbmat = cbmat.transpose() # a and b are reversed for plotting
    mfe_df = [mfe for cblist in cbmat for mfe in cblist]
    ## create SPLASH fram items
    if spmat.any():
        spmat = spmat.transpose() # a and b are reversed for plotting
        sps_df = [sps for splist in spmat for sps in splist]
    else:
        sps_df = ["-" for i in range(len(mfe_df))]
    ## RNA position
    pRNAa_df = [i+1 for i in range(m) for j in range(n)]
    pRNAb_df = [i+1 for i in list(range(n))*m]
    ## RNA string
    RNAa_df = [i for i in [s for s in RNAa] for j in range(n)]
    RNAb_df = [s for s in RNAb]*m
    csize = [8]*len(mfe_df)
    ## reverse the plot numbering
    if opt["var_rvp"]:
        pRNAa_df = pRNAa_df[::-1]
        pRNAb_df = pRNAb_df[::-1]
    ## create data dictionary
    df = DataFrame({"mfe":mfe_df, "posA":pRNAa_df, "posB":pRNAb_df, "nameA":RNAa_df, "nameB":RNAb_df, "csize":csize, "SPLASH":sps_df})
    ## shape dfs
    hoverA, hoverB = "",""
    #for i,shape_dict in enumerate(shape_list,1):
    #    typeaX, typebX = strainA.split(":")[0], strainB.split(":")[0]
    #    shapeA_list, shapeB_list = shape_dict[f"{nameA}:{typeaX}"], shape_dict[f"{nameB}:{typebX}"]
    #    shapeA_df = [i for i in [s for s in shapeA_list] for j in range(n)]
    #    shapeB_df = [s for s in shapeB_list]*m
    #    df[f"shapeA{i}"] = shapeA_df
    #    df[f"shapeB{i}"] = shapeB_df
    #    hoverA += f" - shape{i}: @shapeA{i}"
    #    hoverB += f" - shape{i}: @shapeB{i}"
    ## remove greater than var_cdc entries 
    df = df[df.mfe < opt["var_bkm"]]
    df = df.sort_values(by=['mfe'], ascending=False)
    ## set bokeh parameters
    hover = HoverTool(tooltips=[
        (f"{nameB} {strainB}", f"pos: @posB - base: @nameB{hoverB}"),
        (f"{nameA} {strainA}", f"pos: @posA - base: @nameA{hoverA}"),
        (f"mfe", "@mfe kcal/mol"),
        (f"SPLASH", "@SPLASH")])
    left, right, bottom, top = [0, m+1, 0, n+1]
    plt = figure(plot_width=int((m+1)*opt["var_bks"]), plot_height=int((n+1)*opt["var_bks"]),
                 x_range=(left, right), y_range=(top, bottom),
                 ## extract range
                 #x_range=(306, 328),           y_range=(67, 45),
                 x_axis_label=f"{nameA} {strainA}", y_axis_label=f"{nameB} {strainB}",
                 tools=[hover, ResetTool(), BoxZoomTool(), WheelZoomTool(),\
                        PanTool(), SaveTool()], toolbar_location="above")
    ## set boundaries
    plt.x_range.bounds = (left, right)
    plt.y_range.bounds = (bottom, top) # reversed
    ## add color bar
    mlow = min(df.mfe) if not df.mfe.empty else -1.0
    mhigh = max(df.mfe) if not df.mfe.empty else 0.0
    color_mapper = LinearColorMapper(palette="Viridis256", low=mlow, high=mhigh,
                                     low_color="blue", high_color="red")
    color_bar = ColorBar(color_mapper=color_mapper, ticker=BasicTicker(),
                     label_standoff=12, border_line_color=None, location=(0,0))
    plt.add_layout(color_bar, "right") 
    # create data source
    data_dict = {"posA":df["posA"], "posB":df["posB"], "nameA":df["nameA"],
                 "nameB":df["nameB"], "mfe":df["mfe"], "csize":df["csize"],
                 "SPLASH":df["SPLASH"]}
    #for i,shape_dict in enumerate(shape_list,1):
    #    data_dict[f"shapeA{i}"] = df[f"shapeA{i}"]
    #    data_dict[f"shapeB{i}"] = df[f"shapeB{i}"]
    source = ColumnDataSource(data=data_dict)
    plt.circle("posA", "posB", fill_color={"field":"mfe", "transform":color_mapper},
               size="csize", fill_alpha=1, source=source, line_color=None)
    ## create slider
    mlow = min(df.mfe) if not df.mfe.empty else -1.0
    mhigh = max(df.mfe) if not df.mfe.empty else 0.0
    mfe_slider = Slider(start=mlow, end=mhigh, value=mhigh, 
                        step=0.01, title="kcal/mol", direction="rtl")
                        #, 
                        #callback_policy="mouseup")
                        #, callback=callback)
        #var mfe_val = cb_obj.value;
    callback = CustomJS(args=dict(source=source, mfe=mfe_slider), code="""
        var csize_data = source.data["csize"];
        var mfe_data = source.data["mfe"]
        var mfe_val = mfe.value;
        for (var i = 0; i < mfe_data.length; i++) {
            if (mfe_data[i] <= mfe_val) {
                csize_data[i] = 8;
            } else {
                csize_data[i] = 0;
            }
        }
        source.change.emit();
        """)
    #source.change.emit()
    mfe_slider.js_on_change("value_throttled", callback)
        #source.trigger("change");""")
    ############################################################################
    ### create slider
    #callback = CustomJS(args=dict(source=source), code="""
    #    var mfe_val = cb_obj.value;
    #    var data = source.data;
    #    for (i = 0; i < data["posA"].length; i++) {
    #        if (data["mfe"][i] <= mfe_val) {
    #            data["csize"][i] = 8;
    #        } else {
    #            data["csize"][i] = 0;
    #        }
    #    }
    #    source.change.emit();""")
    #    #source.trigger("change");""")
    #mlow = min(df.mfe) if not df.mfe.empty else -1.0
    #mhigh = max(df.mfe) if not df.mfe.empty else 0.0
    #mfe_slider = Slider(start=mlow, end=mhigh, value=mhigh, 
    #                    step=0.01, title="kcal/mol", direction="rtl", 
    #                    js_event_callbacks=callback)
    #                    #callback_policy="mouseup", callback=callback)
    ############################################################################
    ## add weight plot option
    #weightOptions = RadioButtonGroup(labels=["no weight", f"{vRNAa} {typea}", f"{vRNAa} {typea}", f"{vRNAa} {typea}", f"{vRNAa} {typea}"], active=0)
    ## add RNA bases plot
    xa, ya = [[i for i in RNAa], [i for i,s in enumerate(RNAa)]]
    xb, yb = [[i for i in RNAb], [i for i,s in enumerate(RNAb)]]
    if opt["var_rvp"]:
        RNAlsta = list(RNAa)[::-1]
        RNAlstb = list(RNAb)
    else:
        RNAlsta = list(RNAa)
        RNAlstb = list(RNAb)[::-1]
    RNAa_source = ColumnDataSource(data=dict(
        posA=[i-0.2 for i in list(range(1,m+1))], posB=[0.3]*m, 
        text=RNAlsta))
    RNAb_source = ColumnDataSource(data=dict(
        posA=[0.5]*n, posB=[i+0.2 for i in list(range(1,n+1))[::-1]], 
        text=RNAlstb))
    RNAa_plt = Plot(plot_width=int((m+1)*opt["var_bks"]),   plot_height=30, 
                    x_range=plt.x_range, y_range=Range1d(0,2), 
                    tools=[hover, ResetTool(), BoxZoomTool(), WheelZoomTool(),\
                      PanTool(), SaveTool()], toolbar_location=None)
    RNAb_plt = Plot(plot_width=30,  plot_height=int((n+1)*opt["var_bks"]), 
                    x_range=Range1d(0,2),        y_range=plt.y_range, 
                    tools=[hover, ResetTool(), BoxZoomTool(), WheelZoomTool(),\
                      PanTool(), SaveTool()], toolbar_location=None)
    RNAa_glyph = Text(x="posA", y="posB", text="text", text_color="black", text_font_size="10pt")
    RNAb_glyph = Text(x="posA", y="posB", text="text", text_color="black", text_font_size="10pt")
    RNAa_plt.add_glyph(RNAa_source, RNAa_glyph)
    RNAb_plt.add_glyph(RNAb_source, RNAb_glyph)
    ## dummy
    dummy_source = ColumnDataSource(data=dict(posA=[], posB=[],text=[]))
    dummy_plt = Plot(plot_width=30, plot_height=30, x_range=Range1d(0,2),\
                     y_range=Range1d(0,2), tools=[hover, ResetTool(), BoxZoomTool(), WheelZoomTool(),\
                      PanTool(), SaveTool()], toolbar_location=None)
    dummy_glyph = Text(x="posA", y="posB", text="text", text_color="black", text_font_size="10pt")
    dummy_plt.add_glyph(dummy_source, dummy_glyph)
    dummy_plt.outline_line_alpha = 0.0
    ## arrange plots
    xA, xB = f"{nameB}-{strainB}", f"{nameA}-{strainA}"
    if opt["var_bkv"]:
        plt.output_backend="svg"
        RNAa_plt.output_backend="svg"
        RNAb_plt.output_backend="svg"
        dummy_plt.output_backend="svg"
        gplt = gridplot([[RNAb_plt, plt], [dummy_plt, RNAa_plt]], merge_tools=True)
        lay = layout([gplt,mfe_slider])
        output_file(os.path.join(opt["var_pfx"], f"{xA}_{xB}.html"), title=f"{xA} {xB}")
        save(lay)
        ## export svgs
        export_svgs(gplt, filename=os.path.join(opt["var_pfx"], f"{xA}_{xB}.svg"))
    else:
        gplt = gridplot([[RNAb_plt, plt], [dummy_plt, RNAa_plt]], merge_tools=False)
        lay = layout([gplt,mfe_slider])
        output_file(os.path.join(opt["var_pfx"], f"{xA}_{xB}.html"), title=f"{xA} {xB}")
        save(lay)




################################################################################
## structure function
################################################################################

def createStructures(cand_list, **opt):
    ## create varna structures
    var_pfx = opt["var_pfx"]
    ext = ""
    if opt["var_nex"] == "peak": ext += f"_peak{opt['var_vre']}"
    if opt["var_nex"] == "mfe":  ext += f"_mfe{opt['var_vrm']}"
    if opt["var_tra"]:           ext += "_intra"
    opt["var_pfx"] = os.path.join(var_pfx, f"varna_plots{ext}")
    makeDir(**opt)
    pool_list = [(lk, opt) for lk in cand_list if lk.peak <= opt["var_vre"]]
    with mp.Pool(processes=opt["var_thr"]) as p:
        p.starmap(createStructureMulti, pool_list)
    opt["var_pfx"] = var_pfx

def createStructureMulti(lk, opt):
    ## create 2. Structures with varna
    if opt["var_rvp"]: lk.ai, lk.aj, lk.bi, lk.bj = lk.alen-lk.ai, lk.alen-lk.aj+1, lk.blen-lk.bi, lk.blen-lk.bj+1
    Aname, Bname = f"{lk.aSeq}-{lk.aType}", f"{lk.bSeq}-{lk.bType}"
    RNA = lk.RNA
    Rlen = len(RNA)
    hstring, x, cstring = "", 1, {"U":"#ffff00","C":"#ff0000","A":"#0000ff","G":"#00ff00"}
    for i,n in enumerate(RNA):
        if n == "&": x = 0
        else: hstring += f"{i+x}-{i+x}:fill={cstring[n]},outline={cstring[n]},radius=1;"
    astring = f"{Aname}:type=B,anchor={1};{Bname}:type=B,anchor={Rlen-1};"
    ## base style
    ystring, pstring = list(), list()
    for b,s in cstring.items():
        #ystring.append(f"fill=#ff0000,outline={s},label=#ffffff")
        ystring.append(f"fill=#ffffff,outline={s},label=#000000")
        pstring.append(",".join([str(i) for i,r in enumerate(RNA.replace("&",""),1) if r == b]))
    ## do varna
    file_name = os.path.join(opt["var_pfx"], f"{Aname}_{lk.ai+1}-{lk.aj}_{Bname}_{lk.bi+1}_{lk.bj}_{lk.mfe:.2f}")
    title_name = f"{Aname} {lk.ai+1}-{lk.aj} {Bname} {lk.bi+1}-{lk.bj} ({lk.mfe:.2f} kcal/mol)"
    doVARNA(file_name, RNA, lk.structure, title_name, hstring, astring, ystring, pstring, 10, **opt)
    if opt["var_vrp"] and opt["var_vrd"] and opt["var_vrf"]: svg2pdf(file_name, **opt)
    ## do single varna
    if opt["var_vrs"]:
        file_name = os.path.join(opt["var_pfx"], f"{Aname}_{lk.ai+1}-{lk.aj}_{Bname}_{lk.bi+1}_{lk.bj}_single")
        title_name = f"{Aname} {lk.ai+1}-{lk.aj} ({lk.a_mfe:.2f} kcal/mol) {Bname} {lk.bi+1}-{lk.bj} ({lk.b_mfe:.2f} kcal/mol)"
        pattern = f"{lk.a_structure}&{lk.b_structure}"
        doVARNA(file_name, RNA, pattern, title_name, hstring, astring, ystring, pstring, 10, **opt)
        if opt["var_vrp"] and opt["var_vrd"] and opt["var_vrf"]: svg2pdf(file_name, **opt)

def doVARNA(file_name, RNA, constraint, title_name, hstring, astring, ystring, pstring, period, **opt):
    ## calls VARNA
    var_vrd, var_vrp, var_vra, var_vrf = opt["var_vrd"], opt["var_vrp"], opt["var_vra"], opt["var_vrf"]
    C_call=["java", "-cp", var_vrp, "fr.orsay.lri.varna.applications.VARNAcmd",
            "-sequenceDBN", RNA, "-structureDBN", constraint, "-o", f"{file_name}.svg",
            "-algorithm", var_vra, "-periodNum", f"{period}",
            "-spaceBetweenBases", "0.72", "-title", title_name,
            "-highlightRegion" , hstring, "-annotations", astring,
            "-backbone", "#000000", "-bp", "#000000", "-bpStyle", "-lw"]
    if ystring and pstring:
        for i,(ys,ps) in enumerate(zip(ystring,pstring),1):
            C_call.extend([f"-basesStyle{i}", ys, f"-applyBasesStyle{i}on", ps])
    if var_vrp and var_vrf:
        call(C_call, shell=False)
    with open(f"{file_name}.py", "w") as cvar:
        f_name = os.path.basename(file_name)
        C_call[9] = f"{f_name}.svg"
        if not var_vrp:
            C_call[2] = "VARNAv3-93.jar"
        x_call="', '".join(C_call)
        cvar.write(f"#!/usr/bin/env python3\n"+\
                   f"from subprocess import call\n"+\
                   f"C_call=['{x_call}']\n"+\
                   f"call(C_call, shell=False)\n"+\
                   f"C_call=['rsvg-convert', '-f', 'pdf', '-o', '{f_name}.pdf', '{f_name}.svg']\n"+\
                   f"call(C_call, shell=False)")
                   #f"C_call=['{var_vrd}', '-D', '{f_name}.svg', '--without-gui', '--export-pdf={f_name}.pdf']\n"+\
    #print(f"Status: VARNA call written to \"{file_name}.sh\"")
    #rsvg-convert -f pdf -o mygraph.pdf mygraph.svg

def svg2pdf(file_name, **opt):
    ## calls librsvg2
    C_call=[rsvg-convert, "-f", "pdf", "-o", f"{file_name}.pdf", f"{file_name}.svg"]
    #C_call=[opt["var_vrd"], "-D", f"{file_name}.svg", "--without-gui", f"--export-pdf={file_name}.pdf"]
    call(C_call, shell=False)




################################################################################
## circos plots
################################################################################

def createCircosData(cand_list, data_dict, **opt):
    ## create circos data
    var_pfx = opt["var_pfx"]
    ext = ""
    if opt["var_nex"] == "peak":  ext += f"_peak{opt['var_cie']}"
    if opt["var_nex"] == "mfe":   ext += f"_mfe{opt['var_cim']}"
    if opt["var_cit"] == "inter": ext += f"_inter"
    if opt["var_cit"] == "intra": ext += f"_intra"
    ## do mutant plots
    if opt["var_ciu"]:
        mutant_seq_circos = [(m.split(",")[0]) for m in opt["var_ciu"].split(";")]
        mutant_typ_circos = [(m.split(",")[1]) for m in opt["var_ciu"].split(";")]
        mutant_circos = {m.split(",")[0]:m.split(",")[1] for m in opt["var_ciu"].split(";")}
        mnameext = "-".join(mutant_typ_circos)
        ext += f"_{mnameext}"
        mutant_typ_circos.append("WT")
    if opt["var_cir"] and len(opt["var_cir"].split("-")) == 3:
        rSeq, rStart, rEnd = opt["var_cir"].split("-")
        ext += f"_{rSeq}-{rStart}-{rEnd}"
    opt["var_pfx"] = os.path.join(var_pfx, f"circos_plots{ext}")
    makeDir(**opt)
    karyo_file = os.path.join(var_pfx, "circos_karyo.tsv")
    bands_file = os.path.join(var_pfx, "circos_bands.tsv")
    # sort segments
    segs_dict = {k:v for k,v in sorted({(n,s):len(r) for (n,s),r in data_dict.items() if s == "WT"}.items(), key=operator.itemgetter(1), reverse=True)}
    seg_dict = dict()
    if opt["var_ciu"]:
        for (k,v),l in segs_dict.items():
            if k in mutant_circos.keys():
                seg_dict[(k,mutant_circos[k])] = l
            else:
                seg_dict[(k,v)] = l
    else:
        seg_dict = segs_dict
    # prepare data
    cover_dict = {(n,s):[0 for i in range(r)] for (n,s),r in seg_dict.items()}
    link_list = list()
    for lk in cand_list:
        if opt["var_ciu"]:
            if lk.aType.split("-")[0] == "WT" and lk.aSeq in mutant_seq_circos: continue
            if lk.bType.split("-")[0] == "WT" and lk.bSeq in mutant_seq_circos: continue
            if lk.aType.split("-")[0] not in mutant_typ_circos: continue
            if lk.bType.split("-")[0] not in mutant_typ_circos: continue
        else:
            if lk.aType.split("-")[0] != "WT" or lk.bType.split("-")[0] != "WT": continue
        if len(lk.aType.split("-")) > 1 or len(lk.bType.split("-")) > 1: continue
        if lk.peak > opt["var_cie"] or lk.mfe > opt["var_cim"]: continue
        if (opt["var_cit"] == "intra" and lk.aSeq != lk.bSeq): continue
        if (opt["var_cit"] == "inter" and lk.aSeq == lk.bSeq): continue
        if opt["var_cir"] and len(opt["var_cir"].split("-")) == 3:
            rSeq, rStart, rEnd = opt["var_cir"].split("-")
            rStart, rEnd = int(rStart), int(rEnd)
            if opt["var_rev"]:
                if rSeq == lk.aSeq: rStart, rEnd = lk.alen-rEnd, lk.alen-rStart
                if rSeq == lk.bSeq: rStart, rEnd = lk.blen-rEnd, lk.blen-rStart
            else:
                rStart, rEnd = rStart-1, rEnd-1
            if lk.aSeq == rSeq and (lk.aj < rStart or lk.ai > rEnd): continue
            if lk.bSeq == rSeq and (lk.bj < rStart or lk.bi > rEnd): continue
        # create histogram data
        # removed -1
        for a in range(lk.ai,lk.aj): cover_dict[(lk.aSeq,lk.aType.split("-")[0])][a] += 1
        for b in range(lk.bi,lk.bj): cover_dict[(lk.bSeq,lk.bType.split("-")[0])][b] += 1
        # create link list
        link_list.append(((lk.aSeq,lk.aType.split("-")[0]), lk.ai, lk.aj, (lk.bSeq,lk.bType.split("-")[0]), lk.bi, lk.bj))
    makeHistogram(cover_dict, **opt)
    # create colors for at least 12 segments, other will be random
    color_dict = dict()
    for i,((name,strain),seg_len) in enumerate(seg_dict.items()):
        color_dict[(name,strain)] = f"seg{i+1}"
    num_colors = len(color_dict.keys())
    # create karyo and bands
    makeKaryoBands(seg_dict, color_dict, **opt)
    # create link files
    makeCircosLinks(link_list, color_dict, **opt)
    # do circos
    circosConfig("all", num_colors, **opt)
    doCircos("all", **opt)
    if opt["var_cis"]:
        for i,((name,strain),seg_len) in enumerate(seg_dict.items()):
            circosConfig(f"{name}_{strain}", num_colors, **opt)
            doCircos(f"{name}_{strain}", **opt)
    opt["var_pfx"] = var_pfx

def doCofold(RNA, constraint, **opt):
    ## do Cofold
    cvar.dangles = opt["var_dng"]
    cvar.noLP = int(opt["var_nlp"])
    cvar.temperature = opt["var_tem"]
    fc = fold_compound(RNA)
    fc.constraints_add(constraint, CONSTRAINT_DB | CONSTRAINT_DB_DEFAULT)
    pattern, mfe = fc.mfe_dimer()
    return mfe, pattern

def makeCircosLinks(link_list, color_dict, **opt):
    ## create circos all links file
    links_out = os.path.join(opt["var_pfx"], f"circos_links_all.tsv")
    with open(links_out, "w") as linkout:
        for ((nA,sA),ai,aj,(nB,sB),bi,bj) in link_list:
            linkout.write(f"{nA}_{sA} {ai} {aj-1} {nB}_{sB} {bi} {bj-1} color={color_dict[(nA,sA)]}\n")
    ## create circos single links file
    if opt["var_cis"]:
        for (name,strain),seg in color_dict.items():
            links_out = os.path.join(opt["var_pfx"], f"circos_links_{name}_{strain}.tsv")
            with open(links_out, "w") as linkout:
                for ((nA,sA),ai,aj,(nB,sB),bi,bj) in link_list:
                    if opt["var_cit"] == "intra" and opt["var_cid"] > bi - aj + 1: continue
                    if (nA,sA) == (name,strain):
                        linkout.write(f"{nA}_{sA} {ai} {aj-1} {nB}_{sB} {bi} {bj-1} color={color_dict[(nB,sB)]}\n")
                    elif (nB,sB) == (name,strain):
                        linkout.write(f"{nA}_{sA} {ai} {aj-1} {nB}_{sB} {bi} {bj-1} color={color_dict[(nA,sA)]}\n")

def makeHistogram(cover_dict, **opt):
    ## creates a histogram coverage file for circos
    histo_file = os.path.join(opt["var_pfx"], f"circos_histo.tsv")
    with open(histo_file, "w") as histout:
        for (name,strain),cover_list in cover_dict.items():
            ist, cst = 0, cover_list[0]
            for i,c in enumerate(cover_list):
                if c == cst:
                    continue
                else:
                    histout.write(f"{name}_{strain} {ist} {i-1} {float(cst)}\n")
                    ist, cst = i, c
            histout.write(f"{name}_{strain} {ist} {i} {cst}\n")

def makeKaryoBands(seg_dict, color_dict, **opt):
    ## create circos karyo file
    karyo_file = os.path.join(opt["var_pfx"], f"circos_karyo.tsv")
    bands_file = os.path.join(opt["var_pfx"], f"circos_bands.tsv")
    with open(karyo_file, "w") as outkaryo, open(bands_file, "w") as outbands:
        for i,((name,strain),seg_len) in enumerate(seg_dict.items()):
            n,s = name.replace("_","-"), strain.replace("_","-")
            outbands.write(f"{name}_{strain} 0 {seg_len} {n}-{s}\n")
            outkaryo.write(f"chr - {name}_{strain} {n}-{s} 0 {seg_len} {color_dict[(name,strain)]}\n")
            outkaryo.write(f"band {name}_{strain} 1 1 0 {seg_len} {color_dict[(name,strain)]}\n")

def addColors(num_colors, colors):
    ## random color generator
    seed(1337)
    for i in range(13,num_colors+1):
        colors += f"""
        seg{i} = {randint(0,255)},{randint(0,255)},{randint(0,255)}"""
    return colors

def doCircos(fname, **opt):
    ## start circos and create plots
    circos_conf = os.path.abspath(os.path.join(opt["var_pfx"], f"circos_config_{fname}.conf"))
    wd = os.getcwd()
    os.chdir(opt["var_pfx"])
    C_call = [opt["var_cip"], "-conf", circos_conf]
    call(C_call, shell=False)
    os.rename("circos.png", f"circos_plot_{fname}.png")
    os.rename("circos.svg", f"circos_plot_{fname}.svg")
    os.chdir(wd)

def circosConfig(fname, num_colors, **opt):
    # need: circos_links, circos_histo
    ## tick config
    ticks = f"""
        radius = 1r
        color = black
        thickness = 2p
        multiplier = 1
        format = %d
        <tick>
            spacing = 2u
            size = 10p
        </tick>
        <tick>
            spacing = 10u
            size = 15p
            show_label = yes
            label_size = 14p
            label_offset = 7p
            format = %d
            label_parallel = yes
        </tick>"""
    ## ideogram config
    ideogram = f"""
        <spacing>
            default = 0.01r
        </spacing>
        radius = 0.80r
        thickness = 40p
        fill = yes
        stroke_color = dgrey
        stroke_thickness = 2p
        show_label = yes
        show_bands = yes
        fill_bands = yes
        band_transparency = 0
        band_stroke_color = black
        label_font = default
        label_radius = dims(image,radius) - 60p;
        label_size = 30
        label_parallel = yes"""
    ## color config
    colors = f"""
        seg1 = 141,211,199
        seg2 = 255,255,179
        seg3 = 190,186,218
        seg4 = 251,128,114
        seg5 = 128,177,211
        seg6 = 253,180,98
        seg7 = 179,222,105
        seg8 = 252,205,229
        seg9 = 217,217,217
        seg10 = 188,128,189
        seg11 = 204,235,197
        seg12 = 255,237,111"""
    if num_colors > 12:
        colors = addColors(num_colors, colors)
    ## link config
    links = f"""
        <link>
            file = circos_links_{fname}.tsv
            color = vlgrey
            radius = 0.88r
            bezier_radius = 0.1r
            thickness = 4
        </link>"""
    ## plot config
    plots = f"""
        <plot>
            type = text
            color = red
            file = circos_bands.tsv
            r0 = 1.008r
            r1 = 1.3r
            label_size = 24
            label_font = condensed
        </plot>
        <plot>
            type = histogram
            file = circos_histo.tsv
            r1 = 0.98r
            r0 = 0.88r
            max = 20
            min = 0
            stroke_type = outline
            thickness = 1
            <backgrounds> <background>
                color = vvlgrey
            </background> </backgrounds>
            <axes> <axis>
                spacing = 0.2r
                color = lgrey
                thickness = 1
            </axis> </axes>
        </plot>"""
    ## main config file
    config = f"""# main circos config file
        show_ticks = yes
        show_tick_labels = yes
        
        <ticks> {ticks}
        </ticks>
        <ideogram> {ideogram}
        </ideogram>
        <colors> {colors}
        </colors>
        
        karyotype = circos_karyo.tsv
        chromosomes_units = 100
        
        <links> {links}
        </links>
        <plots> {plots}
        </plots>
        
        # basic settings files
        <image> 
            <<include etc/image.conf>>
            radius* = 600p
        </image>
        <<include etc/colors_fonts_patterns.conf>> 
        <<include etc/housekeeping.conf>> """
    ## write circos config file
    circos_conf = os.path.join(opt["var_pfx"], f"circos_config_{fname}.conf")
    with open(circos_conf, "w") as circonf:
        circonf.write(config)




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
