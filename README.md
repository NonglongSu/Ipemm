## Synopsis
Ioemm (Indel-phased evolutionary Markov model) is a markov model that designed for estimating the evolutionary parameters, including {6 $\sigma$s, $\omega$, $\tau$} from GTR+MG94 model, gap opening for three indel phases, gap extension weight (insertion/deletion).   


## Dependencies
* [R 4.0 +]  (https://www.r-project.org/) and libraies `Biostrings`, `stringr`, `seqinr`, `stringi`, `plyr`, `R.utils`, `Matrix`, `expm`, `SQUAREM`, `jsonlite`, `matlib`, `ggplot2`, `ggrepel`, `tidyverse`, `stats`, `dfoptim`, `purrr`, `profvis`
* [GNU Make] (https://www.gnu.org/software/make/)
* Download coati (Codon-Aware Multiple Sequence Alignments) (https://github.com/jgarciamesa/coati)

## Download
git clone https://github.com/NonglongSu/Ipemm
cd chapter4          

## Usage  
* Generate 100 gillespie-simulated folders(500 files each).   
* EM-training on aligned-simulated data.
* EM-training & importance sampling on raw-simulated data.

```
Usage: make help    show all useful commands 
```
