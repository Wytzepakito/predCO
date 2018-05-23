# predCO
prediction of meiotic crossovers in plants 

Tested on Ubuntu 16.04.4 LTS

The dependencies are:

python v3.4+,
python v2.7+,
R v3.3.2,
bedtools v2.25,
fimo (MEME suite), scikit-learn, pandas, numpy, scipy, stats, matplotlib (python v3 modules),
Bio (python module),
DNAshapeR (R v3.3.2 package)

How to use:

bash predCO.sh [-h] [-p <true|false>] [-g <true|false>] [-n <sample_gene|genome|pericentromere] -r <reference_genome_folder> -i <input_folder> -o <output_folder> -s <path/to/scripts> -d <path/to/R-3.3.2/Rscript>
