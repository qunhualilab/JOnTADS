# JOnTADS
A unified caller for TADs and stripes in various Hi-C data

Dependencies:

pip install numba

pip install numpy

pip install scipy

pip install sklearn

pip install qpsolvers==1.10.0

Go to the repository of this folder

Test run:

Single sample run:

python JOnTAD.py -F ./data/chr18.csv -O ./results/chr18.csv.tad

Multiple samples run:

python JOnTAD.py -F ./data/ES_rep1.chr18 ./data/ES_rep2.chr18 ./data/ME_rep1.chr18 ./data/ME_rep2.chr18 ./data/MS_rep1.chr18 ./data/MS_rep2.chr18 ./data/NP_rep1.chr18 ./data/NP_rep2.chr18 ./data/TP_rep1.chr18 ./data/TP_rep2.chr18 -O ./results/ES_rep1.chr18.tad ./results/ES_rep2.chr18.tad ./results/ME_rep1.chr18.tad ./results/ME_rep2.chr18.tad ./results/MS_rep1.chr18.tad ./results/MS_rep2.chr18.tad ./results/NP_rep1.chr18.tad ./results/NP_rep2.chr18.tad ./results/TP_rep1.chr18.tad ./results/TP_rep2.chr18.tad

Stripe calling:

python get_stripe.py -F ./data/chr18.csv -O ./results/chr18.csv.stripe -C 18
