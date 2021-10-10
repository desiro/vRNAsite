# vRNAsite - manual

## Usage
```
vRNAsite.py -fsa <in_vRNAsite_fasta> -pfx <out_prefix> [options]
```

## Version
```
vRNAsite.py 0.0.1 (alpha)
```

## Dependencies
```python v3.7.1```, ```numpy v1.16.4```, ```pandas v1.3.3```, ```bokeh v2.4.0```, ```ViennaRNA v2.4.13```, ```matplotlib v3.4.3```, ```VARNA v3.93```, ```circos v0.69.8```

## Description
```vRNAsite``` can predict long-range RNA-RNA interaction between any two or more RNA sequences. The tool has been written in ```Python 3.7.1``` and relies heavily on the ```RNAcofold``` python site-package of the ```ViennaRNA Package 2.4.13```.

## Options

### Main
```
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
```

### SPLASH Data
```
--splashData,-spd
    optionally read SPLASH, SHAPE or similar data from matrices in directory, names have to match the fasta names (default: )

--splashReads,-spr
    set mean minimum threshold for SPLASH, SHAPE or similar data (default: 0.0)
```

### IAV Weight Matrix
```
--weightMatrix,-wgm
    use weight matrix to consider the relative positioning (default: False)

--weightDistance,-dst
    distance threshold modifier for links to be in the same area from the center
    of their respective RNAs (default: 0.1)

--weightDescent,-dsc
    defines how the weight should descent (default: 0.0)

--saveOrientation,-svo
    also save orientated matrices, this will use up more space (default: False)
```

### Watershed Segmentation
```
--clusterPeak,-clp
    maximum peak energy for clusters (default: -10.0)

--clusterBarrier,-clb
    maximum barrier energy for cluster (default: -10.0)

--clusterValue,-clv
    sets the cluster extraction value for plotting (default: 1.0)
```

### RNA Structure Prediction
```
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
```

### VARNA Plots
```
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
```

### Contact Matrix Plots
```
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
```

### bokeh Plots
```
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
```

### Peak Distribution Plots
```
--distributionPlots,-dsh
    enables peak energy distribution printing (default: False)

--distributionSize,-dss
    defines the distribution plot size modifier (default: 1.0)

--distributionSvg,-dsv
    also plot svg distribution plots (default: False)
```

### Circos Plots
```
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
    specify the interaction type which should be shown with circos (default: all) (choices: all,intra,inter)

-circosIntraDist,-cid
    define the minimum distance of intra circos intreactions (default: 10)

--circosRange,-cir
    set circos position plot range; sequence, start and end position has to be
    divided by a minus symbol (default: )

--circosMutants,-ciu
    create mutant circos plots; the mutant segments can be defined by this parameter;
    sequenes and type have to be separated by a comma, multiple mutants by semicolon (default: )
```
