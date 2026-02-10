# fst-gl

A program to parse vcf/bcf files and estimate allele frequencies and fst using Phred-scaled genotype likelihoods (PL).

## Installation
### Prerequisites:
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


### Build Instructions
Run ```nim build``` to compile an executable binary ```bin/fst-gl```.  

You can run ```nim clean``` to clean the build environment.

#### Build Notes:
HTSlib needs to be accessible on ```$LD_LIBRARY_PATH```. This can cause problems if it is 
installed in a conda environment because it may cause conflicts with system gcc and libstdc++.  
Here is one solution to this. We simply add the relevant path from our conda environment to ```$LD_LIBRARY_PATH```.
```bash
export LD_LIBRARY_PATH=/lib64:~/micromamba/envs/$CONDA_DEFAULT_ENV/lib/
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
  -m, --remove-missing           Remove samples with missing data? Requires DP
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



## Program Output
```fst-gl``` outputs two separate output files. One contains the _F<sub>ST</sub>_ estimates and the other contains the allele and genotype frequency estimates for each population. By default, these files are named 'out-af.txt.gz' and 'out-fst.txt.gz', this can be changed by using  the ```--output``` option.  

If we run the command on the test data (provided in the test directory)
```./bin/fst-gl -t 2 --output test/test --pop-file test/popfile.txt test/allpops.bcf```

The fst output is shown here
```bash
zcat test/test-fst.txt.gz | head | column -t
```

```
chrom  pos    fst       a         b         all_converged
chr1   4905   0.016957  0.003974  0.234356  true
chr1   6714   0.210428  0.040930  0.194509  true
chr1   7643   0.085329  0.008863  0.103869  true
chr1   9617   0.016707  0.003882  0.232341  true
chr1   10339  0.152228  0.023040  0.151352  true
chr1   12816  0.195082  0.036294  0.186043  true
chr1   13784  0.043413  0.010670  0.245771  true
chr1   14320  0.044174  0.010873  0.246148  true
chr1   16535  0.345063  0.088364  0.256080  true
```
|Column | Description| 
|---|---|
|chrom | chromosome  |
|pos | position on chromosome   |
|fst | fst estimate  |
|a | numerator a eqn 2   |
|abc | denominator a+b+c eqn 2 - 4  |
|all_converged | did the EM algorithm converge for all populations?  |


```
zcat test/test-af.txt.gz | head | column -t
```
```
chrom  pos   pop  n   ngood  mleaf     mlegf                       converged  n_iter  starts
chr1   4905  p0   50  50     0.390518  0.358570,0.501824,0.139606  true       3       0.333333,0.333333,0.333333
chr1   4905  p1   50  50     0.429524  0.239875,0.661204,0.098922  true       3       0.333333,0.333333,0.333333
chr1   4905  p2   50  50     0.282088  0.475701,0.484422,0.039877  true       3       0.333333,0.333333,0.333333
chr1   6714  p0   50  50     0.368995  0.339859,0.582292,0.077849  true       3       0.333333,0.333333,0.333333
chr1   6714  p1   50  50     0.339983  0.380011,0.560013,0.059976  true       3       0.333333,0.333333,0.333333
chr1   6714  p2   50  50     0.000000  1.000000,0.000000,0.000000  true       3       0.333333,0.333333,0.333333
chr1   7643  p0   50  50     0.160178  0.679645,0.320355,0.000000  true       3       0.333333,0.333333,0.333333
chr1   7643  p1   50  50     0.180312  0.679655,0.280066,0.040279  true       3       0.333333,0.333333,0.333333
chr1   7643  p2   50  50     0.000000  1.000000,0.000000,0.000000  true       3       0.333333,0.333333,0.333333
```

|Column | Description| 
|---|---|
|chrom | chromosome  |
|pos | position on chromosome  | 
|pop | population   |
|n | number of individuals used for estimation   |
|ngood | number of individuals with data   |
|mleaf | estimated allele frequencies   |
|mlegf | estimated genotype frequencies HomRef,Het,HomAlt   |
|converged | did the EM algorithm converge?   |
|n_iter | number of iterations until convergence   |
|starts | starting values for EM algorithm   |