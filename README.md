# vRNAsite

```vRNAsite``` can predict long-range RNA-RNA interaction between any two or more RNA sequences. The tool has been written in ```Python 3.7.1``` and relies heavily on the ```RNAcofold``` python site-package of the ```ViennaRNA Package 2.4.13```. For command-line options, please refer to the [manual](https://github.com/desiro/vRNAsite/blob/master/manual.md)

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

## Installation with Miniconda

Please download the [Miniconda3](https://docs.conda.io/en/latest/miniconda.html) application for your system. The following will demonstrate the installation and set up of miniconda on Linux, which should be similar on other platforms.

```
bash Miniconda3-latest-Linux-x86_64.sh -p ~/miniconda3
conda create --name vRNAsite
conda activate vRNAsite
conda install -c bioconda viennarna
conda install -c conda-forge numpy=1.19.1
conda install -c conda-forge bokeh=1.3.1
conda install -c conda-forge matplotlib=3.3.1
conda install -c conda-forge pandas=1.1.2
conda install -c bioconda circos=0.69.8 
```

## Examples

The basic input for ```vRNAsite.py``` includes the following parameters:
* ```-pfx``` - name of the output folder
* ```-fsa``` - input fasta file

### Fasta Header Formation

Wild-type sequences should have the format ```>segment_name WT```. The ```segment_name``` should not contain spaces, minus signs, colons or underscores. After a wild-type sequence there can be any number of mutants in the format ```>segment_name mutant_name:start-end``` with the same naming restrictions for ```segment_name``` and ```mutant_name```. The ```segment_name``` must match the wild type segment of the mutant. The ```start``` and ```end``` should be the start and end nucleotides in the wild-type sequence, which should be replaced by the following mutant sequence snippet

### Basic Example

```
python3 vRNAsite.py -pfx example -fsa example.fa -thr 4 -ovr -rev -cmp -rvp -nex peak -sql -clp -13.0
```

### Example with Plots

```
python3 vRNAsite.py -pfx example -fsa example.fa -thr 4 -ovr -rev -cmp -rvp -nex peak  -sql -clp -13.0 -cih -cis -cie -10.0 -vrh -vrs -vre -10.0 -plh -plv -pls 2.0 -bkh -bks 0.5
```

## Authors

* [Daniel Desirò](https://github.com/desiro)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Reference

Please cite [vRNAsite](https://doi.org/10.1101/424002) if you find our tool useful.

## Workflow overview

![workflow](https://github.com/desiro/vRNAsite/blob/master/methods_workflow.png "(a) creates all subsequences from two sequence (b) predicts structures with RNAcofold between subsequences (c) averages the free energy value for any two bases from involved substructure predictions and creates a contact matrix (d) extracts contac boundaries sequences from the contact matrix with a watershed segmentation algorithm (e) predicts structures with RNAcofold from the extracted contact boundaries (f) creates heat-plots, clust-plots and  cand-plots")
