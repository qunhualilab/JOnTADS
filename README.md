# JOnTADS
JOnTADS is a versatile tool for identifying Topologically Associating Domains (TADs) and stripes in various chromatin conformation capture data, including population Hi-C, single-109
cell Hi-C and micro-C. It allows for easy analysis of Hi-C data across multiple samples and outputs results in a structured format.

## Dependencies
```sh
pip install numba==0.56.4  
pip install numpy==1.23.5  
pip install scipy==1.12.0  
pip install qpsolvers==2.7.3
```

## Usage
Go to the repository of this folder

### Test Run
#### Single Sample Run
python JOnTADS.py -F ./data/chr18.csv -O ./results/chr18.csv.tad  

#### Multiple Samples Run
python JOnTADS.py -F ./data/ES_rep1.chr18 ./data/ES_rep2.chr18 ./data/ME_rep1.chr18 ./data/ME_rep2.chr18 ./data/MS_rep1.chr18 ./data/MS_rep2.chr18 ./data/NP_rep1.chr18 ./data/NP_rep2.chr18 ./data/TP_rep1.chr18 ./data/TP_rep2.chr18 -O ./results/ES_rep1.chr18.tad ./results/ES_rep2.chr18.tad ./results/ME_rep1.chr18.tad ./results/ME_rep2.chr18.tad ./results/MS_rep1.chr18.tad ./results/MS_rep2.chr18.tad ./results/NP_rep1.chr18.tad ./results/NP_rep2.chr18.tad ./results/TP_rep1.chr18.tad ./results/TP_rep2.chr18.tad  

#### Stripe Calling
python JOnTADS.py -F ./data/chr18.csv -O ./results/chr18.csv.tad --stripe_output ./results/chr18.csv.stripe -C 18 --stripe True  
