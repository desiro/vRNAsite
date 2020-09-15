#!/usr/bin/env python3
# script: compare.py
# author: Daniel Desiro'
# dependencies: numpy, RNAfold, VARNAv3-93.jar, inkscape
script_usage="""
usage
    compare.py -f <in_fasta> -p <out_prefix> [options]

version
    compare.py 0.0.1 (alpha)

dependencies
    x

description
    Validates significance of flexible regions in a sequence.

--prefix,-pfx
    output directory and prefix for result files

--summary,-sum
    output summary table (default: )

--vRNAsite1,-vs1
    first vRNAsite table

--vRNAsite2,-vs2
    second vRNAsite table

--SPLASH1,-sp1
    first SPLASH table

--SPLASH2,-sp2
    second SPLASH table

--candidateRPM,-cdr
    define minimum RPM value (default: 5000)

--candidateMFE,-cdm
    define minimum MFE value (default: -10.0)

--candidateOverlap,-cdo
    define the number of overlapping bases (default: 5)

--orientation,-ort
    only allow specific orientations of the vRNPs; A = arbitrary, U = upright,
    R = reverse (default: A) (choices: A,U,R,UR,AU,AR,AUR)

--overwrite,-ovr
    overwrite data with folder named in prefix (default: False)

--minlength,-mle
    minimum length for an interaction (default: 5)

--mfeSplash,-msp
    use SPLASH mfe values for comparison (default: False)

reference
    Reference.
"""

import argparse as ap
import sys
import os
import re
import time
from operator import attrgetter, itemgetter
from itertools import combinations, product, combinations_with_replacement
from statistics import mean, median


### old


import random
try:
    from RNA import fold_compound, cvar, CONSTRAINT_DB, CONSTRAINT_DB_DEFAULT, bp_distance
except:
    print("Error: The mandatory library ViennaRNA 2.4 is missing, please read the README.md for more information!")
    exit()
from itertools import combinations, product, combinations_with_replacement
from multiprocessing import Process, Manager
from math import ceil, floor
from numpy import zeros
from pandas import DataFrame, set_option
# plotting
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from numpy import arange, mean, prod, array, matrix, zeros, nonzero
from numpy import sum as npsum
import pickle
from numpy import arange, mean, prod, array, asarray
import itertools
from subprocess import Popen, PIPE, call
from multiprocessing import Pool, Process, Manager, Lock
from copy import deepcopy
from itertools import permutations, product, combinations_with_replacement, combinations, chain
from collections import Counter

import networkx as nx
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt




################################################################################
## main
################################################################################

def main(opt):
    ############################################################################
    ## load fasta file
    time_s = getTime()
    print(f"Status: Read tables ...")
    tab1_list = readTable(opt["var_vs1"], os.path.splitext(os.path.basename(opt["var_vs1"]))[0], **opt)
    tab2_list = readTable(opt["var_vs2"], os.path.splitext(os.path.basename(opt["var_vs2"]))[0], **opt)
    tab3_list = readTable(opt["var_sp1"], os.path.splitext(os.path.basename(opt["var_sp1"]))[0], **opt)
    tab4_list = readTable(opt["var_sp2"], os.path.splitext(os.path.basename(opt["var_sp2"]))[0], **opt)
    ############################################################################
    ## create output folder
    opt = makeDir(**opt)
    ############################################################################
    ## compare data sets 103 69 215
    print(f"Status: Compare data sets ...")
    name = os.path.basename(opt["var_pfx"])
    compareData(tab1_list, tab2_list, tab3_list, tab4_list, f"{name}_vRNAsite-Rep1+Rep2_SPLASH-Rep1+Rep2_{opt['var_ort']}.tsv", **opt)
    time_s = getTime(time_s, f"Compare data sets")
    return opt["var_pfx"]




################################################################################
## functions
################################################################################

def readTable(table, tname, **opt):
    ## read table file
    data_type, tab_list, conv = F"vRNAsite", list(), {"PB2":1, "PB1":2, "PA":3, "HA":4, "NP":5, "NA":6, "M":7, "NS":8}
    with open(table, "r") as tabin:
        header = next(tabin)
        if "RPM" in header.strip().split(): data_type  = "SPLASH"
        lnk = dict()
        for i,line in enumerate(tabin):
            line = line.strip().split()
            if "NA" in line: continue
            if data_type == "SPLASH":
                if int(line[0]) > int(line[3]):
                    lnk["aSeg"], lnk["ai"], lnk["aj"], lnk["aOrt"] = line[3:6]+["A"]
                    lnk["bSeg"], lnk["bi"], lnk["bj"], lnk["bOrt"] = line[0:3]+["A"]
                    lnk["bk"], lnk["bl"], lnk["ak"], lnk["al"] = line[9:13]
                else:
                    lnk["aSeg"], lnk["ai"], lnk["aj"], lnk["aOrt"] = line[0:3]+["A"]
                    lnk["bSeg"], lnk["bi"], lnk["bj"], lnk["bOrt"] = line[3:6]+["A"]
                    lnk["ak"], lnk["al"], lnk["bk"], lnk["bl"] = line[9:13]
                lnk["read"], lnk["RPM"] = line[6:8]
                lnk["mfe"] = line[8]
                lnk["c_energy"] = 0.0
                lnk["RNA"], lnk["structure"] = line[13:15]
                lnk["a_mfe"], lnk["a_structure"], lnk["b_mfe"], lnk["b_structure"], lnk["free_mfe"], lnk["free_structure"] = 0.0, ".", 0.0, ".", 0.0, "."
            else:
                e1, e2 = 6, 7
                if len(line[0].split("-")) == 2: line[0] = f"{line[0]}-A"
                if len(line[3].split("-")) == 2: line[3] = f"{line[3]}-A"
                nline = [line[1], line[2], line[4], line[5]] + line[0].split("-") + line[3].split("-")
                nline[4], nline[7] = conv[nline[4]], conv[nline[7]]
                if opt["var_ort"]:
                    if nline[6] not in list(opt["var_ort"]) or nline[9] not in list(opt["var_ort"]): continue
                if nline[4] > nline[7]:
                    lnk["aSeg"], lnk["ai"], lnk["aj"], lnk["aOrt"] = nline[7], nline[2], nline[3], nline[9]
                    lnk["bSeg"], lnk["bi"], lnk["bj"], lnk["bOrt"] = nline[4], nline[0], nline[1], nline[6]
                else:
                    lnk["aSeg"], lnk["ai"], lnk["aj"], lnk["aOrt"] = nline[4], nline[0], nline[1], nline[6]
                    lnk["bSeg"], lnk["bi"], lnk["bj"], lnk["bOrt"] = nline[7], nline[2], nline[3], nline[9]
                lnk["ak"], lnk["al"], lnk["bk"], lnk["bl"] = 0, 0, 0, 0
                lnk["read"], lnk["RPM"] = 0, 0
                lnk["mfe"] = line[7]
                lnk["c_energy"] = line[6]
                lnk["RNA"], lnk["structure"], lnk["a_mfe"], lnk["a_structure"], lnk["b_mfe"], lnk["b_structure"], lnk["free_mfe"], lnk["free_structure"] = line[8:]
            lnk["rank"], lnk["type"] = 0, data_type
            lk = links(**lnk)
            if lk.aj - lk.ai + 1 < opt["var_mle"]: continue
            if lk.bj - lk.bi + 1 < opt["var_mle"]: continue
            tab_list.append(lk)
    if data_type == "SPLASH" and not opt["var_msp"]: tab_list = sorted(tab_list, key=attrgetter("RPM"), reverse=True)
    elif data_type == "SPLASH" and opt["var_msp"]: tab_list = sorted(tab_list, key=attrgetter("mfe"), reverse=False)
    else: tab_list = sorted(tab_list, key=attrgetter("c_energy"))
    for r,t in enumerate(tab_list, 1):
        t.rank = r
        # fix different lengths for better comparison
        # Udorn Rep2 HA 22 AAGGGTGTTTTTCCTTATATTTAATTACTAATAC CCTTATATTT  -> __________ AAGGGTGTTTTT__________AATTACTAATAC
        # Udorn Rep1 HA 1775 CCCTGCTTTTGC_ _ -> T CCCTGCTTTTGCT
        # PR8 Rep2 HA 1306 TGCCGTTACTCCTTTGGTTGTGTTGTG TTT -> ___ TGCCGTTACTCC___GTTTGTGTTGTG
        if tname == "Udorn-Rep2":
            if t.aSeg == 4 and t.ai >= 22: t.ai += 10
            if t.aSeg == 4 and t.aj >= 22: t.aj += 10
            if t.bSeg == 4 and t.bi >= 22: t.bi += 10
            if t.bSeg == 4 and t.bj >= 22: t.bj += 10
        if tname == "Udorn-Rep1":
            if t.aSeg == 4 and t.ai >= 1775: t.ai += 1
            if t.aSeg == 4 and t.aj >= 1775: t.aj += 1
            if t.bSeg == 4 and t.bi >= 1775: t.bi += 1
            if t.bSeg == 4 and t.bj >= 1775: t.bj += 1
        if tname == "PR8-Rep2":
            if t.aSeg == 4 and t.ai >= 1306: t.ai += 3
            if t.aSeg == 4 and t.aj >= 1306: t.aj += 3
            if t.bSeg == 4 and t.bi >= 1306: t.bi += 3
            if t.bSeg == 4 and t.bj >= 1306: t.bj += 3
    return tab_list

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

def comp(tx, x):
    ## compare positions
    if min([t.aj for t in tx]) - max([t.ai for t in tx]) >= x and\
       min([t.bj for t in tx]) - max([t.bi for t in tx]) >= x: return True
    else: return False

def compareData(tabx1_list, tabx2_list, tab3_list, tab4_list, name, **opt):
    ## compare data sets
    ## remove duplicate vRNAsite positions
    tab1_list, tab2_list = list(), list()
    for tx in tabx1_list:
        skip = False
        for ty in tabx1_list:
            if tx.aSeg == ty.aSeg and tx.bSeg == ty.bSeg and comp([tx,ty], 0) and tx.c_energy > ty.c_energy:
                skip = True
                continue
        if not skip:
            tab1_list.append(tx)
    for tx in tabx2_list:
        skip = False
        for ty in tabx2_list:
            if tx.aSeg == ty.aSeg and tx.bSeg == ty.bSeg and comp([tx,ty], 0) and tx.c_energy > ty.c_energy:
                skip = True
                continue
        if not skip:
            tab2_list.append(tx)
    ## create segment combination dictionaries
    write_list, comp_list = list(), list()
    tab_list, td = [tab1_list, tab2_list, tab3_list, tab4_list], list()
    for tl in tab_list:
        dl = {k:list() for k in combinations(range(1,9),2)}
        for t in tl:
            t_list = dl[(t.aSeg,t.bSeg)]
            t_list.append(t)
            dl[(t.aSeg,t.bSeg)] = t_list
        td.append(dl)
    ## compare positions
    for k in combinations(range(1,9),2):
        for t1 in td[0][k]:
            for t2 in td[1][k]:
                if comp([t1,t2], opt["var_cdo"]):
                    for t3 in td[2][k]:
                        if comp([t1,t2,t3], opt["var_cdo"]):
                            for t4 in td[3][k]:
                                if comp([t1,t2,t3,t4], opt["var_cdo"]):
                                    comp_list.append((t1,t2,t3,t4))
    ## compare RPM and MFE values
    for (v1,v2,s1,s2) in comp_list:
        if (v1.c_energy <= opt["var_cdm"] or v2.c_energy <= opt["var_cdm"]) and\
           (s1.RPM >= opt["var_cdr"] or s2.RPM >= opt["var_cdr"]):
            write_list.append((v1,v2,s1,s2,v1.rank,v2.rank,s1.rank,s2.rank))
    write_list = sorted(write_list, key=itemgetter(4,5,6,7))
    out_comp = os.path.join(opt["var_pfx"], name)
    ## write files
    with open(out_comp, "w") as outcomp:
        outcomp.write(f"candidate\tsource\taS\tai\taj\tbS\tbi\tbj\tpeak\tRPM\trank\tmfe\tRNA\tstructure\taOrt\tbOrt\ta_mfe\ta_structure\tb_mfe\tb_structure\tfree_structure\tfree_mfe\tread\tak\tal\tbk\tbl\n")
        for i,(v1,v2,s1,s2,v1.rank,v2.rank,s1.rank,s2.rank) in enumerate(write_list,1):
                outcomp.write(f"c{i}\tvRNAsite-Rep1\t{v1.aSeg}\t{v1.ai}\t{v1.aj}\t{v1.bSeg}\t{v1.bi}\t{v1.bj}\t"+\
                              f"{v1.c_energy}\t{v1.RPM}\t{v1.rank}\t{v1.mfe}\t{v1.RNA}\t{v1.structure}\t{v1.aOrt}\t{v1.bOrt}\t"+\
                              f"{v1.a_mfe}\t{v1.a_structure}\t{v1.b_mfe}\t{v1.b_structure}\t{v1.free_structure}\t"+\
                              f"{v1.free_mfe}\t{v1.read}\t{v1.ak}\t{v1.al}\t{v1.bk}\t{v1.bl}\n")
                outcomp.write(f"\tvRNAsite-Rep2\t{v2.aSeg}\t{v2.ai}\t{v2.aj}\t{v2.bSeg}\t{v2.bi}\t{v2.bj}\t"+\
                              f"{v2.c_energy}\t{v2.RPM}\t{v2.rank}\t{v2.mfe}\t{v2.RNA}\t{v2.structure}\t{v2.aOrt}\t{v2.bOrt}\t"+\
                              f"{v2.a_mfe}\t{v2.a_structure}\t{v2.b_mfe}\t{v2.b_structure}\t{v2.free_structure}\t"+\
                              f"{v2.free_mfe}\t{v2.read}\t{v2.ak}\t{v2.al}\t{v2.bk}\t{v2.bl}\n")
                outcomp.write(f"\tSPLASH-Rep1\t{s1.aSeg}\t{s1.ai}\t{s1.aj}\t{s1.bSeg}\t{s1.bi}\t{s1.bj}\t"+\
                              f"{s1.c_energy}\t{s1.RPM}\t{s1.rank}\t{s1.mfe}\t{s1.RNA}\t{s1.structure}\t{s1.aOrt}\t{s1.bOrt}\t"+\
                              f"{s1.a_mfe}\t{s1.a_structure}\t{s1.b_mfe}\t{s1.b_structure}\t{s1.free_structure}\t"+\
                              f"{s1.free_mfe}\t{s1.read}\t{s1.ak}\t{s1.al}\t{s1.bk}\t{s1.bl}\n")
                outcomp.write(f"\tSPLASH-Rep2\t{s2.aSeg}\t{s2.ai}\t{s2.aj}\t{s2.bSeg}\t{s2.bi}\t{s2.bj}\t"+\
                              f"{s2.c_energy}\t{s2.RPM}\t{s2.rank}\t{s2.mfe}\t{s2.RNA}\t{s2.structure}\t{s2.aOrt}\t{s2.bOrt}\t"+\
                              f"{s2.a_mfe}\t{s2.a_structure}\t{s2.b_mfe}\t{s2.b_structure}\t{s2.free_structure}\t"+\
                              f"{s2.free_mfe}\t{s2.read}\t{s2.ak}\t{s2.al}\t{s2.bk}\t{s2.bl}\n")

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

def getTime(time_s=0, name=""):
    if time_s and name:
        time_e = time.time()-time_s
        time_e = time.strftime("%H:%M:%S", time.gmtime(time_e))
        time_c = time.strftime('%x %X')
        print(f"Status: {name} finished at {time_c} in {time_e}")
    return time.time()




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
