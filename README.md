# [<samp>vRNAsite</samp>](https://github.com/desiro/vRNAsite)
[![License: GPL v3](https://img.shields.io/badge/License-GPL_v3-bd0000.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python v3.9.7](https://img.shields.io/badge/Language-Python_v3-75a8d3.svg)](https://www.python.org/)
[![Conda v4.11.0](https://img.shields.io/badge/Uses-Conda-43b02a.svg)](https://docs.conda.io/en/latest/miniconda.html)

***

## Description

This tool predicts potential intermolecular long-range RNA-RNA interactions between two or more RNA sequences. It can also be applied to predict intramolecular long-range RNA-RNA interactions. The tool predicts short consecutive stable interactions.

### Mandatory Prerequisites

* [![Python v3.9.7](https://img.shields.io/badge/Python_v3.9.7-75a8d3.svg)](https://www.python.org/downloads/release/python-397/)
* [![NumPy v1.22.2](https://img.shields.io/badge/NumPy_v1.22.2-013243.svg)](http://www.numpy.org/)
* [![ViennaRNA v2.5.0](https://img.shields.io/badge/ViennaRNA_v2.5.0-006795.svg)](https://www.tbi.univie.ac.at/RNA/)

### Optional Prerequisites

* [![Conda v4.11.0](https://img.shields.io/badge/Conda_v4.11.0-43b02a.svg)](https://docs.conda.io/en/latest/miniconda.html)
* [![VARNA v3.93](https://img.shields.io/badge/VARNA_v3.93-ffba27.svg)](http://varna.lri.fr/)
* [![Pandas v1.4.0](https://img.shields.io/badge/Pandas_v1.4.0-130654.svg)](https://pandas.pydata.org/)
* [![Bokeh v2.4.2](https://img.shields.io/badge/Bokeh_v2.4.2-542437.svg)](https://docs.bokeh.org/)
* [![Matplotlib v3.5.1](https://img.shields.io/badge/Matplotlib_v3.5.1-11557c.svg)](https://matplotlib.org/)
* [![Circos v0.69.8](https://img.shields.io/badge/Circos_v0.69.8-ec1c24.svg)](http://circos.ca/)

***

## Installation

To run <samp>vRNAsite</samp>, I recommend using Miniconda and following the steps below. If this is the first time using conda, you should probably restart your shell after the installation of Miniconda. The following will demonstrate the installation and set up of Miniconda on Linux, which should be similar on other platforms. For Windows 10 users, I advise using the [Ubuntu 20.04 LTS](https://www.microsoft.com/en-us/p/ubuntu-2004-lts/9n6svws3rx71?cid=msft_web_chart) subsystem. More information can be found on the [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and [Bioconda](https://bioconda.github.io/user/install.html) pages.

### Conda Installation

Installing Miniconda:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Updating Miniconda and setting channels:
```
conda update conda
conda update python
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Installing Conda packages:
```
conda create --name vRNAsite python=3.9.7
conda activate vRNAsite
conda install -c bioconda viennarna=2.5.0
conda install -c conda-forge numpy=1.22.2
conda install -c lb_arrakistx varna=3.93
conda install -c conda-forge pandas=1.4.0
conda install -c conda-forge bokeh=2.4.2
conda install -c conda-forge matplotlib=3.5.1
conda install -c bioconda circos=0.69.8
git clone https://github.com/desiro/vRNAsite.git
cd vRNAsite
```

### Installation without Conda

Installing the ViennaRNA package on Linux:
```
tar -zxvf ViennaRNA-2.5.0.tar.gz
cd ViennaRNA-2.5.0
./configure --with-python3
make
sudo make install
```

Installing the ViennaRNA package on MAC:
```
tar -zxvf ViennaRNA-2.5.0.tar.gz
cd ViennaRNA-2.5.0
./configure --enable-universal-binary --with-python3
make
sudo make install
```

***

## Examples

The basic input for ```vRNAsite.py``` includes the following parameters:
* ```-pfx``` - name of the output folder
* ```-fsa``` - input fasta file, uses special fasta header

### Fasta Header Formation

Wild-type sequences should have the format ```>segment_name WT```. The ```segment_name``` should not contain spaces, minus signs, colons or underscores. After a wild-type sequence there can be any number of mutants in the format ```>segment_name mutant_name:start-end``` with the same naming restrictions for ```segment_name``` and ```mutant_name```. The ```segment_name``` must match the wild type segment of the mutant. The ```start``` and ```end``` should be the start and end nucleotides in the wild-type sequence, which should be replaced by the following mutant sequence snippet

### Basic Example

```
python vRNAsite.py -pfx example -fsa example.fa -thr 4 -ovr -rev -cmp -rvp -nex peak -clp -10.0
```

### Example with Plots

```
python vRNAsite.py -pfx example -fsa example.fa -thr 4 -ovr -rev -cmp -rvp -nex peak -clp -10.0 -cih -cis -cie -10.0 -vrh -vrs -vre -10.0 -plh -plv -pls 2.0 -bkh -bks 0.5
```

### Options

For more command line options, see the [manual](https://github.com/desiro/vRNAsite/blob/master/manual.md).

***

## Authors

* [Daniel Desirò](https://github.com/desiro)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Reference

Please cite <samp>vRNAsite</samp> if you find our tool useful.

```
D. Desirò, E. Barth, B. Hardin, M. Schwemmle and M. Marz.
"vRNAsite: Prediction of packaging signals in segmented RNA viruses at the example of influenza A virus indicating flexible RNA-RNA interactions between segments."
In Preparation, 2022.
```
