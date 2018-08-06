# predCO
prediction of meiotic crossovers in plants 

## Citation

Demirci,  S., Peters, S.A., Ridder, D. de and Dijk, A.D.J.  van. (2018) DNA sequence and shape are predictive for meiotic crossovers throughout the plant kingdom. Plant J., 95, 686-699.

https://doi.org/10.1111/tpj.13979

### The dependencies are:

* python v3.4+,
* python v2.7+,
* R v3.3.2,
* bedtools v2.25,
* fimo (MEME suite), scikit-learn, pandas, numpy, scipy, stats, matplotlib (python v3 modules),
* Bio (python module),
* DNAshapeR (R v3.3.2 package)

Tested on Ubuntu 16.04.4 LTS

### How to use:
``` bash
bash predCO.sh [-h] [-p <true|false>] [-g <true|false>] [-n <sample_gene|genome|pericentromere] -r <reference_genome_folder> -i <input_folder> -o <output_folder> -s <path/to/scripts> -d <path/to/R-3.3.2/Rscript>
```
