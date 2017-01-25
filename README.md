# variantDetection

This is a toolkit to perform population based filtering of variants called from NGS data

### Prerequisites
* python
* tqdm (python module)

### Installation
```
git clone --recursive
cd variantDetection/htslib
make
cd ../
```

### usage
Two variables are set to call alleles:
 - The minimum read depth to call alleles at
 - The minimum proportion of reads at a position with an allele

To run the pipeline you will need a subdirectory
