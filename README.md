# JOnTADS: a unified caller for TADs and stripes in Hi-C data
JOnTADS is a versatile tool for identifying Topologically Associating Domains (TADs) and stripes in various chromatin conformation capture data, including population Hi-C, single-cell Hi-C and micro-C. It allows for easy analysis of Hi-C data across multiple samples.

## Dependencies
```sh
pip install numba==0.56.4  
pip install numpy==1.23.5  
pip install scipy==1.12.0  
pip install qpsolvers==2.7.3
```

## Installation
```sh
pip install JOnTADS
```

## Test Run
Suppose you are at the folder of this README file.
#### Single Sample Run
To identify TADs from a single sample, use the following command:
```sh
JOnTADS -F ./data/chr18.csv -O ./results/chr18.csv.tad
```

#### Multiple Samples Run
To analyze multiple samples simultaneously, use:
```sh
JOnTADS -F ./data/ES_rep1.chr18 ./data/ES_rep2.chr18 ./data/ME_rep1.chr18 ./data/ME_rep2.chr18 ./data/MS_rep1.chr18 ./data/MS_rep2.chr18 ./data/NP_rep1.chr18 ./data/NP_rep2.chr18 ./data/TP_rep1.chr18 ./data/TP_rep2.chr18 -O ./results/ES_rep1.chr18.tad ./results/ES_rep2.chr18.tad ./results/ME_rep1.chr18.tad ./results/ME_rep2.chr18.tad ./results/MS_rep1.chr18.tad ./results/MS_rep2.chr18.tad ./results/NP_rep1.chr18.tad ./results/NP_rep2.chr18.tad ./results/TP_rep1.chr18.tad ./results/TP_rep2.chr18.tad
```

#### Stripe Calling
To call stripes in addition to TADs:
```sh
JOnTADS -F ./data/chr18.csv -O ./results/chr18.csv.tad --stripe True --stripe_output ./results/chr18.csv.stripe --chr 18
```

## Usage
You can run JOnTADS in your terminal via
```sh
JOnTADS -F ./data/chr18.csv -O ./results/chr18.csv.tad --stripe True --stripe_output ./results/chr18.csv.stripe --chr 18
```
or
```sh
python JOnTADS.py -F ./data/chr18.csv -O ./results/chr18.csv.tad --stripe True --stripe_output ./results/chr18.csv.stripe --chr 18
```
or in your python file via
```sh
from JOnTADS import JOnTADS, get_stripe
tads = JOnTADS(file_names, outputs)
stripes = get_stripe(file_name, chr_num)
```

## Input and Output Format

#### Input Format

The input files should be Hi-C contact matrices separated by spaces or commas.

_In progress_: we are working on supporting supporting additional input formats.

#### Output Format

The output files contain information about the identified TADs or stripes. 

For TAD calling, the output contains four columns:
```sh
start, end, TAD score, TAD size
```
For stripe calling, the output contains six columns
```sh
chr, x1, x2, chr, y1, y2
```
where the stripe extends from `(x1, y1)` to `(x2, y2)`.

## Parameters

- **`-F`**: Input file(s) with Hi-C data.
- **`-O`**: Output file(s) for the detected TADs.
- **`-MAXSZ`**: Maximum size of TADs allowed, default 200.
- **`-MINSZ`**: Minimum size of TADs allowed, default 7.
- **`--stripe`**: Set to `True` to enable stripe detection.
- **`--chr`**: (When `stripe' is set to True) Chromosome number for stripe calling, e.g. 18.
- **`--stripe_output`**: (When `stripe' is set to True) Output file for stripe calling results.

## Contact

Feel free to contribute to the project by opening issues or pull requests. Any feedback or suggestions are highly appreciated. Correspondence should be addressed to qunhua.li@psu.edu. You can also contact the maintainer qiuhai.stat@gmail.com.

## Happy analyzing your Hi-C data with JOnTADS!
