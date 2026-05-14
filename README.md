# fst-gl

A program to parse vcf/bcf files and estimate allele frequencies and fst using Phred-scaled genotype likelihoods (PL).

## Installation
We provide a binary that can be found in the [releases](https://github.com/AndyJasonowicz-IPHC/fst-gl/releases). This release bundles htslib and the compression libraries needed to run ```fst-gl```. It does require OpenMP to be installed.

### Build from source
```fst-gl``` can be built from source as well.
#### Prerequisites:
* clang
* OpenMP
* nim 

Install the prerequisites using your Linux distribution's package manager. For RPM based systems run.
```bash 
dnf install -y clang libgomp 
```


Install Nim following the instructions here https://nim-lang.org/install_unix.html.
```bash
curl https://nim-lang.org/choosenim/init.sh -sSf | sh
```

#### Required Nim Libraries:
* zip/gzipfiles
* hts
* Therapist


#### Build Instructions
The easiest way to build fst-gl is to run ```nimble bundle``` to compile an executable binary ```bin/fst-gl```. This will automatically download dependencies and statically link htslib and libdeflate.

Run ```nimble build``` to compile an executable binary ```bin/fst-gl``` that dynamically links htslib (you will need to install this on your own). You may need to set ```$LD_LIBRARY_PATH``` if htslib is in a non standard location (see [build notes](#build-notes)).  

You can run ```nimble wipe``` to clean the build environment.

#### Build Notes:
When building from source HTSlib may need to be accessible on ```$LD_LIBRARY_PATH```. If you get an error that says ```could not load: libhts.so``` you may need to set your ```$LD_LIBRARY_PATH``` like this.
```bash
export LD_LIBRARY_PATH=/lib64:<path to directory that contains libhts.so>
```
If HTSlib is installed in a conda environment, we can add relevant path from our conda environment to ```$LD_LIBRARY_PATH```.
```bash
export LD_LIBRARY_PATH=/lib64:$CONDA_PREFIX/lib/
```


## Program usage
```bash
fst-gl --help
```
```bash
fst-gl - a program to estimate FST from genotype likelihoods.

Usage:
  fst-gl <BCF>
  fst-gl (--version | -h | --help)

Arguments:
  <BCF>                          bcf/vcf.gz/vcf file to read

Options:
  -p, --pop-file=<pop-file>      File contining sample names and population
                                 labels
  -o, --output=<output>          Prefix for output files [default: out]
  -r, --region=<region>          Specify a region to process [chr]:[start-stop].
                                 [default: *]
  -m, --allow-missing            Allow samples with missing data? Requires DP
                                 field in BCF. [default: false]
  -s, --random-starts            Use random starting values for EM algorithm?
                                 [default: false]
  -b, --boots=<boots>            Number of bootstrap resamplings
      --replace                  Perform re-sampling with replacement? [default:
                                 false]
  -t, --threads=<threads>        Number of threads to use in parallel sections
                                 of code. [default: 1]
  -c, --chunk-size=<chunk-size>  Number of SNPs to read into memory at a time.
                                 [default: 50000]
      --em-tol=<em-tol>          Tolerance required for stopping EM estimation
                                 of allele frequencies. [default: 1e-8]
      --max-iter=<max-iter>      Maximum iterations allowed for EM estimation of
                                 allele frequencies. [default: 200]
  -v, --verbose...               Verbose. Can be specified multiple times for
                                 more detailed output.
      --version                  Prints version
  -h, --help                     Show help message
```

## Usage notes
If program fails to converge, try increasing --max-iter to 500 or even 1000, this may be helpful in cases where sequencing depth is very low.

## Program Output
```fst-gl``` outputs two separate output files. One contains the _F<sub>ST</sub>_ estimates and the other contains the allele and genotype frequency estimates for each population. By default, these files are named 'out-af.txt.gz' and 'out-fst.txt.gz', this can be changed by using  the ```--output``` option.  

If we run the command on the test data (provided in the test directory). Note we can provide the ```-v``` flag multiple times to see additional output that may be useful to us.
```bash
mkdir test
./bin/fst-gl -t 2 -vvv --output test/test-results --pop-file example-data/popfile.txt example-data/allpops.bcf
```

The fst output is shown here
```bash
zcat test/test-results-fst.txt.gz | head | column -t
```

```
chrom  pos    fst        a          b          c         all_converged
chr1   4905   -0.010161  -0.002154  0.015152   0.199020  true
chr1   6714   0.144033   0.020136   0.020481   0.099184  true
chr1   9617   0.009810   0.002391   -0.116669  0.358056  true
chr1   10339  0.134766   0.018761   -0.010084  0.130537  true
chr1   12816  0.292023   0.065551   -0.009672  0.168592  true
chr1   13784  -0.000948  -0.000232  -0.013211  0.258208  true
chr1   14320  -0.002548  -0.000594  0.057547   0.176088  true
chr1   16535  0.387910   0.103612   0.003925   0.159566  true
chr1   17748  0.049507   0.012100   0.012368   0.219931  true
```

|Column | Description| 
|---|---|
|chrom | chromosome  |
|pos | position on chromosome   |
|fst | fst estimate  a/(a+b+c) |
|a | a from wc 84  |
|b | b from wc 84  |
|c | c from wc 84  |
|all_converged | did the EM algorithm converge for all populations?  |
|p | p-value only if ```--boots``` is used|
|converged_boots | number of resamplings in which EM algorithm converged for all populations (this value is used for p-value calculation) |

The allele and heterozygosity estimates are shown here
```bash
zcat test/test-results-af.txt.gz | head | column -t
```
```
chrom  pos   pop  n   ngood  mleaf     mlegf                       converged  n_iter  starts
chr1   4905  p0   50  22     0.318102  0.681154,0.001489,0.317357  true       105     0.333333,0.333333,0.333333
chr1   4905  p1   50  19     0.301583  0.396922,0.602991,0.000087  true       66      0.333333,0.333333,0.333333
chr1   4905  p2   50  17     0.295282  0.409898,0.589639,0.000463  true       83      0.333333,0.333333,0.333333
chr1   6714  p0   50  24     0.290627  0.419259,0.580230,0.000512  true       82      0.333333,0.333333,0.333333
chr1   6714  p1   50  11     0.181471  0.811181,0.014696,0.174123  true       562     0.333333,0.333333,0.333333
chr1   6714  p2   50  21     0.000090  0.999820,0.000180,0.000000  true       29      0.333333,0.333333,0.333333
chr1   9617  p0   50  15     0.467622  0.106766,0.851223,0.042011  true       327     0.333333,0.333333,0.333333
chr1   9617  p1   50  18     0.425115  0.272434,0.604902,0.122664  true       125     0.333333,0.333333,0.333333
chr1   9617  p2   50  18     0.348591  0.305304,0.692210,0.002486  true       236     0.333333,0.333333,0.333333
```

|Column | Description| 
|---|---|
|chrom | chromosome  |
|pos | position on chromosome  | 
|pop | population   |
|n | number of individuals used for estimation   |
|ngood | number of individuals with data   |
|mleaf | estimated allele frequency (ALT)   |
|mlegf | estimated genotype frequencies HomRef,Het,HomAlt   |
|converged | did the EM algorithm converge?   |
|n_iter | number of iterations until convergence   |
|starts | starting values for EM algorithm   |