# vRNAsite - manual

## Mandatory Prerequisites

* [python 3.7.1](https://www.python.org/downloads/release/python-385/)
* [viennaRNA 2.4.13](https://www.tbi.univie.ac.at/RNA/documentation.html#install)
* [numpy 1.16.4](http://www.numpy.org/)

## Optional Prerequisites

* [Miniconda3](https://docs.conda.io/en/latest/miniconda)
* [VARNA 3.93](http://varna.lri.fr/)
* [Bokeh 1.3.1](https://docs.bokeh.org/en/2.2.1/docs/installation.html)
* [Circos 0.69.8](http://circos.ca/software/download/)
* [matplotlib 3.3.1](https://matplotlib.org/users/installing.html)
* [pandas 1.1.2](https://pandas.pydata.org/getting_started.html)
* [Inkscape 0.92](https://inkscape.org/en/)



conda install -c bioconda viennarna
conda install -c conda-forge numpy=1.19.1
conda install -c conda-forge bokeh=1.3.1
conda install -c conda-forge matplotlib=3.3.1
conda install -c conda-forge pandas=1.1.2
conda install -c bioconda circos=0.69.8 
```

## usage
```
vRNAsite.py -fsa <in_vRNAsite_fasta> -pfx <out_prefix> [options]
```

## version
```
vRNAsite.py 0.0.1 (alpha)
```

## dependencies
```python v3.7.1```, ```numpy v1.16.4```, ```pandas v1.1.2```, ```bokeh v1.3.1```, ```ViennaRNA v2.4.13```, ```matplotlib v3.3.1```, ```VARNA v3.93```, ```circos v0.69.8```
```

## description
```vRNAsite``` can predict long-range RNA-RNA interaction between any two or more RNA sequences. The tool has been written in ```Python 3.7.1``` and relies heavily on the ```RNAcofold``` python site-package of the ```ViennaRNA Package 2.4.13```.

## options
```
################################################################

--prefix,-pfx
    output directory and prefix for result files

--fasta,-fsa
    fasta file with all mutants; one entry for each RNA from wild-type and all 
    mutants; will create all possible combinations between mutants and 
    wild-type; naming: \">name WT\" or \">name mut_name:start-end\"; start and 
    end defines which part in the WT should be replaced with the mutant string

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

--noLP,-nlp
    disable lonely pairs for RNAcofold (default: False)

--bigGenome,-big
    use this option for an alternative multiprocessing algorithm for big 
    genomes (default: False)

--onlyPlotting,-opl
    only to plotting, requires the candidates.pcl and the matrices.pcl to be
    in the prefix directory (default: False)

--reversePositions,-rvp
    reverse all positions, useful for negative stranded RNA (default: False)

--sequenceLength,-sql
    include the sequence length at the output table (default: False)

--namingExtension,-nex
    use the this parameter as an extension for naming (default: none) 
    (choices: none,peak,mfe,)

################################################################

--splashData,-spd
    optionally read SPLASH data from SPLASH tables in directory, names have to 
    match the fasta names (default: )

--splashLn,-spn
    use natural logarithm for SPLASH data (default: False)

--splashLog,-spl
    set mean minimum read log count threshold for SPLASH (default: 0.0)

--splashReads,-spr
    set mean minimum read count threshold for SPLASH (default: 0)

################################################################

--weightMatrix,-wgm
    use weight matrix to consider the relative positioning (default: False)

--weightDistance,-dst
    distance threshold modifier for links to be in the same area from the 
    center of their respective RNAs (default: 0.1)

--weightDescent,-dsc
    defines how the weight should descent (default: 0.0)

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

--candidateSingle,-cds
    disable single strain folds for each candidate (default: True)

--candidateLoad,-cdl
    load candidates from file to reduce extraction calculations (default: 
    False)

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
    defines the VARNA drawing algorithm (default: radiate) (choices: 
    line,circular,radiate,naview)

--varnaPeakEnergy,-vre
    maximum peak average free energy for a structure to be plotted with VARNA 
    (default: -10.0)

--varnaMinimumEnergy,-vrm
    maximum minimum free energy for a structure to be plotted with VARNA 
    (default: -10.0)

--varnaPath,-vrp
    use this VARNA path; example: VARNAv3-93.jar (default: VARNAv3-93.jar)

--varnaPdf,-vrd
    use this Inkscape path to create pdf files; example: inkscape
    (default: )

################################################################

--plotHeats,-plh
    plot interaction heat maps (default: False)

--plotProcesses,-plp
    turn off multi processing for matrix plotting (default: True)

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
    specify the interaction type which should be shown with circos (default: 
    all) (choices: all,intra,inter)

--circosRange,-cir
    set circos position plot range; start and end position has to be divided 
    by a minus symbol (default: )

################################################################
```