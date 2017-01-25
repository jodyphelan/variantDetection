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

To run the pipeline you will need VCFs and a special pileup file (created using htsbox).
They should be deposited in folders named fastq and pileup.

The final structure will look like the following:
```
base_dir/
├──pileup/
│  ├──sample1_1.pileup.gz
│  └──sample1_1.pileup.gz.tbi
└──vcf/
   └──sample1.vcf.gz
```

Two variables are set to call alleles:
 - The minimum read depth to call alleles at
 - The minimum proportion of reads at a position with an allele

To run the pipeline run the command:
```
python /path/to/filter_variants.py raw <samples> <base_dir> <min_depth> <read_prop>
python /path/to/filter_variants/py filter
```

