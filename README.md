# fst-gl

A program to parse vcf/bcf files and estimate genotype and allele frequencies and fst using Phred-scaled genotype likelihoods (PL).

## Installation
We provide a binary that can be found in the [releases](https://github.com/AndyJasonowicz-IPHC/fst-gl/releases).

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
The easiest way to build fst-gl is to run ```nimble bundle``` to compile an executable binary ```bin/fst-gl```. This will automatically download dependencies and statically link HTSlib and libdeflate.

Run ```nimble build``` to compile an executable binary ```bin/fst-gl``` that dynamically links HTSlib (you will need to install this on your own). You may need to set ```$LD_LIBRARY_PATH``` if HTSlib is in a non standard location (see [build notes](#build-notes)).  

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
Program options are shown below, please see our [tutorial](https://github.com/AndyJasonowicz-IPHC/fst-gl/wiki/Tutorial-Using-Included-Example-Data) for more information on usage.
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
