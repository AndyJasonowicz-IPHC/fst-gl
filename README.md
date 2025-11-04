A program to parse vcf/bcf files and estimate allele frequencies and fst using genotype likelihoods (PL).


USER REQUIREMENTS:
clang
OpenMP
GCC
nim 

Install required libraries using your package manager.
```bash 
dnf install -y clang libgomp 
```

Install nim 
```bash
curl https://nim-lang.org/choosenim/init.sh -sSf | sh
```

Nim Libraries:
zip/gzipfiles
hts
Therapist


to build 
```bash
nim build 
```

to remove library run ```nim clean```

HTSlib needs to be accesablle on $LD_LIBRARY_PATH but this can cause problems if it is 
installed in a conda env becasue it may cause conflicts with system gcc and libstdc++. 
Here is one solution to this.
```bash
export LD_LIBRARY_PATH=/lib64:~/micromamba/envs/$CONDA_DEFAULT_ENV/lib/
```


#### calc-fst usage
```bash
fst-gl --help
```
```bash
fst-gl - a program to estimate FST from genotype likelihoods.

Usage:
  fst-gl <BCF>
  fst-gl (--version | -h | --help)

Arguments:
  <BCF>                        BCF/VCF file to read

Options:
  -p, --popFile=<popFile>      File contining sample names and population labels
  -o, --output=<output>        Prefix for output files [default: out]
  -r, --region=<region>        Specify a region to process [chr]:[start-stop].
                               [default: *]
  -m, --removeMissing          Remove samples with missing data? Requires DP
                               field in BCF. [default: false]
  -s, --randomStarts           Use random starting values for EM algorithm?
                               [default: false]
  -b, --boots=<boots>          Number of bootstrap resamplings
      --replace                Perform re-sampling with replacement? [default:
                               false]
  -t, --threads=<threads>      Number of threads to use in parallel sections of
                               code. [default: 1]
  -c, --chunkSize=<chunkSize>  Number of SNPs to read into memory at a time.
                               [default: 50000]
      --emTol=<emTol>          Tolerance required for stopping EM estimation of
                               allele frequencies. [default: 1e-8]
      --maxIter=<maxIter>      Maximum iterations allowed for EM estimation of
                               allele frequencies. [default: 200]
  -v, --verbose...             Verbose. Can be specified multiple times for more
                               detailed output.
      --version                Prints version
  -h, --help                   Show help message
```



