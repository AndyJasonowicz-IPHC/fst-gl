# fst-gl

A program to parse vcf/bcf files and estimate genotype and allele frequencies and fst using Phred-scaled genotype likelihoods (PL).

## Installation
We provide a binary that can be found in the [releases](https://github.com/AndyJasonowicz-IPHC/fst-gl/releases). You may need to make the binary file executable first by running ```chmod +x fst-gl``` before first use. 

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
The easiest way to build fst-gl is to run ```nimble bundle``` to compile an executable binary ```bin/fst-gl```. This will automatically download dependencies and _statically link_ HTSlib and libdeflate.

Running ```nimble build``` will compile an executable binary ```bin/fst-gl``` that _dynamically links_ HTSlib (you will need to install this on your own). You may need to set ```$LD_LIBRARY_PATH``` if HTSlib is in a non standard location (see [build notes](#build-notes)).  

You can run ```nimble wipe``` to clean the build environment.

#### Build Notes:
When building from source, HTSlib may need to be accessible on ```$LD_LIBRARY_PATH```. If you get an error that says ```could not load: libhts.so``` you may need to set your ```$LD_LIBRARY_PATH``` like this.
```bash
export LD_LIBRARY_PATH=/lib64:<path to directory that contains libhts.so>
```
If HTSlib is installed in a conda environment, we can add relevant path from our conda environment to ```$LD_LIBRARY_PATH```.
```bash
export LD_LIBRARY_PATH=/lib64:$CONDA_PREFIX/lib/
```


## Program usage
Program options are shown below. 
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


### Example 

For more usage examples, please see our [tutorial](https://github.com/AndyJasonowicz-IPHC/fst-gl/wiki/Tutorial-Using-Included-Example-Data) Wiki page.


We have included a small example dataset in the ```example-data``` directory.
```bash
ls example-data
```

We include a BCF file (```allpops.bcf```) that consists of 50 individuals randomly sampled from 3 simulated populations with a symmetric migration rate of 0.05 among all populations. Sequence data were simulated to 2.5x depth and genotype likelihoods were estimated using [angsd](https://www.popgen.dk/angsd/index.php/ANGSD). 

```fst-gl``` requires VCF/BCF files that contain Phred-scaled genotype likelihoods stored in the PL field. A tab delimited file that describes what population each individual belongs to is also required to be supplied using the ```--pop-file``` option. The first column of the pop file specifies the individual sample ID that is stored in the VCF file and the second is the population that individual belongs to.

```bash
head -n 3 example-data/example-popfile.txt
```
The file should be structured like this.
```
p0_20   p0
p0_56   p0
p0_68   p0
```

```fst-gl``` outputs two separate output files. One contains the _F<sub>ST</sub>_ estimates and the other contains the allele and genotype frequency estimates for each population. By default, these files are named 'out-af.txt.gz' and 'out-fst.txt.gz', this can be changed by using  the ```--output``` option. If we run the command on the example data we can generate these outputs. Note we can provide the ```-v``` flag multiple times to see additional output that may be useful to us.
```bash
mkdir test
./bin/fst-gl -t 2 -vvv --output test/test-results --pop-file example-data/example-popfile.txt example-data/allpops.bcf
```

The fst output is shown here
```bash
zcat test/test-results-fst.txt.gz | head | column -t
```

```
chrom  pos    fst        a          b          c         all_converged
chr1   2182   0.130856   0.017734   0.021237   0.096553  true
chr1   4326   0.161400   0.034293   0.023220   0.154960  true
chr1   5143   0.135323   0.019332   0.003203   0.120322  true
chr1   7729   0.447450   0.128511   0.015052   0.143645  true
chr1   21726  0.481393   0.141329   -0.011702  0.163956  true
chr1   26982  0.251965   0.055682   -0.027256  0.192566  true
chr1   27301  0.120912   0.014001   0.004658   0.097135  true
chr1   34790  0.419278   0.115286   -0.019269  0.178947  true
chr1   42432  -0.007220  -0.001374  0.023514   0.168183  true
```
<details>
<summary>Click to expand table describing columns of  <em>F<sub>ST</sub></em> output</summary>

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
</details>


The allele and heterozygosity estimates are shown here
```bash
zcat test/test-results-af.txt.gz | head | column -t
```
```
chrom  pos   pop  n   ngood  mleaf     mlegf                       converged  n_iter  starts
chr1   2182  p0   43  43     0.269699  0.583611,0.293379,0.123010  true       8       0.333333,0.333333,0.333333
chr1   2182  p1   42  42     0.187055  0.669980,0.285931,0.044089  true       9       0.333333,0.333333,0.333333
chr1   2182  p2   44  44     0.000005  0.999991,0.000009,0.000000  true       8       0.333333,0.333333,0.333333
chr1   4326  p0   46  46     0.424950  0.398814,0.352473,0.248714  true       6       0.333333,0.333333,0.333333
chr1   4326  p1   49  49     0.345213  0.397011,0.515551,0.087438  true       21      0.333333,0.333333,0.333333
chr1   4326  p2   48  48     0.061235  0.907896,0.061738,0.030366  true       11      0.333333,0.333333,0.333333
chr1   5143  p0   47  47     0.253449  0.576870,0.339363,0.083768  true       16      0.333333,0.333333,0.333333
chr1   5143  p1   41  41     0.233253  0.575463,0.382566,0.041970  true       17      0.333333,0.333333,0.333333
chr1   5143  p2   49  49     0.000003  0.999995,0.000005,0.000000  true       8       0.333333,0.333333,0.333333
```

<details>
<summary>Click to expand table describing columns of allele and heterozygosity frequency output</summary>

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
</details>

## Usage notes
If program fails to converge, try increasing --max-iter to 500 or even 1000, this may be helpful in cases where sequencing depth is very low.
