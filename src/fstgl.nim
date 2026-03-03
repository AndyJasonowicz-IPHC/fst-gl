import std/[times, os, math, sequtils, strformat, tables, enumerate, random]
import options
import strutils
import zip/gzipfiles

import hts

when defined(openmp):
    ## use -d:openmp to complile multithreaded version 
    {.passC: "-fopenmp".}
    {.passL: "-fopenmp".}

    {. pragma: omp, header:"omp.h" .}

    proc omp_set_num_threads*(x: cint) {.omp.}
    proc omp_get_num_threads*(): cint {.omp.}
    proc omp_get_max_threads*(): cint {.omp.}
    proc omp_get_thread_num*(): cint {.omp.}

randomize()

type 
    Gl2Beagle = tuple[gls: seq[float32], undeflow: bool]
type
    MLE_AF = tuple[n: int, ngood:int, af: float, gf: array[3, float], converged: bool, iterations: int, starts: array[3, float]]
type
    FST_OUT = tuple[n_pops: float, 
                    n_ave: float, 
                    p_ave: float, 
                    h_ave: float, 
                    he: float, 
                    n_c: float, 
                    s2: float,
                    a: float,
                    b: float,
                    c: float, 
                    fst: float,
                    ]

type
    Query* = ref object
        chrom*: string
        start*: string
        stop*: string
    
proc `$`*(q:Query): string =
  return "Query({q.chrom}:{q.start}-{q.stop})".fmt


proc parse_query_str(x: string, v:VCF): Query =
    ## Parse contig name and length from vcf 
    var contig_table = initTable[string, int]()
    for contig in contigs(v):
        contig_table[contig.name] = contig.length
    
    var qsplit = x.split(':')
    if len(qsplit) != 2:
        stderr.writeLine("Failed to parse query: '{x}'".fmt)
        return nil
    var coords = qsplit[1].split('-')
    if 2 < len(coords):
        stderr.writeLine("Failed to parse query: '{x}'".fmt)
        return nil
    else:
        if qsplit[1] == "": coords = @["1","{contig_table[qsplit[0]]}".fmt ]
        if len(coords) == 1: coords = @[qsplit[1], qsplit[1]]
        if coords[0] == "": coords[0] = "1"  
        if coords[1] == "": coords[1] = "{contig_table[qsplit[0]]}".fmt    

    let q = Query(chrom:qsplit[0], start:coords[0], stop:coords[1])
    return q


# Check for missing values based on type used in vcf
proc is_missing[T:int32|int8|int16|int64](v:T): bool {.inline.} =
    v == T.low


proc reshape_fmt*[T](values: seq[T], n_samples:int): seq[seq[T]] =
    result = newSeq[seq[T]](n_samples)
    let n_per = int(values.len / n_samples)
    for i in 0..<n_samples:
        var off = i * n_per
        for j in off..<off+n_per:
            result[i].add(values[j])
        
        
proc pl2prob(pl: seq[int32], scale: bool=false): seq[float64] = 
    var prob = pl.mapIt(pow(10, -it/10))
    if scale:
        prob = prob.mapIt(it/sum(prob))
    return prob

proc gl_nan_check(x:seq[float64]): bool = 
    var pass = true
    for i in x:
        if i.classify == fcNaN:
            pass = false
        if i.classify == fcInf:
            pass = false
    return pass

proc parse_pop_file(popfile: string): Table[string, string] = 
    ## parse pop_file 
    var pop_map = initTable[string, string]()

    for line in lines popfile:
        var p = line.split('\t')
        if pop_map.hasKeyOrPut(p[0], p[1]):
            pop_map[p[0]] = p[1]
        
    return pop_map


proc maf_em(beagle_gl: seq[seq[float64]], tol=1e-20, max_iter=200, random_starts:bool=true): MLE_AF = 
    var n = beagle_gl.len.float
    var ngood = len(beagle_gl.filterIt(it != @[1/3, 1/3, 1/3]));

    var gf: array[3,float]
    if random_starts: 
        # random starting values
        for i in 0..<3: gf[i] = rand(10.0)
        let isum = sum(gf)
        for i in 0..<3: gf[i] = gf[i]/isum

    else:
        # uniform starting values
        for i in 0..<3: gf[i] = 1.0/3.0
    
    let starts = gf

    var gf_indiv: array[3, float]
    var mle_gf: array[3, float]
    
    var mse_arr: array[3, float]

    var mse = 1 + tol
    var diff = 0.0

    var iteration = 0
    var converged = false

    while (mse > tol) and (iteration < max_iter):
        mle_gf = [0.0,0.0,0.0]
        for samp_gl in beagle_gl:
            for i in 0..<3: gf_indiv[i] = gf[i] * samp_gl[i]
            var prob_data = sum(gf_indiv)
            for i in 0..<3: gf_indiv[i] = gf_indiv[i] / prob_data
            for i in 0..<3: mle_gf[i] += gf_indiv[i]
        
        for i in 0..<3: mle_gf[i] = mle_gf[i] / n
        
        for i in 0..<3: mse_arr[i] = (gf[i] - mle_gf[i])^2
        mse = sum(mse_arr)
        gf = mle_gf        
    
        iteration += 1
    
    if mse < tol:
        converged = true
        #echo "ITERATION: {iteration}, COVERGED: {converged}".fmt
    #else:
    #    echo "WARNING: did not converge after {iteration} iterations.".fmt
        
    var mle_af = [mle_gf[0] + 0.5 * mle_gf[1], mle_gf[2] + 0.5 * mle_gf[1]]
    
    var res: MLE_AF
    res = (n: int(n), ngood: ngood, af: mle_af[1], gf: mle_gf, converged: converged, iterations: iteration, starts: starts)
    return res


proc fst_wc(pfreqs: seq[float64], hfreqs: seq[float64], ns: seq[int]): FST_OUT =  
    var n_pops = float(len(ns))
    var n_ave = sum(ns)/len(ns)
    var p_ave = sum(pfreqs)/float(len(pfreqs))
    var h_ave = sum(hfreqs)/float(len(hfreqs))

    var he = 1-sum([p_ave^2, (1-p_ave)^2])
    #(n_pops*n_ave - sum(sample_sizes^2)/(n_pops*n_ave))/(n_pops-1)
    var n_c = (n_pops*n_ave - sum(ns.mapIt(float(it)^2))/(n_pops*n_ave))/(n_pops-1)
    #s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)
    
    var s2_num: float64
    for i  in 0..int(n_pops)-1:
        s2_num += float(ns[i])*(pfreqs[i] - p_ave)^2

    var s2 = s2_num / ((n_pops-1)*n_ave)
    
    var a=n_ave/n_c*(s2 - 1/(n_ave-1)*(p_ave*(1-p_ave)-((n_pops-1)/n_pops)*s2-(1/4)*h_ave))
    var b=n_ave/(n_ave-1)*(p_ave*(1-p_ave) - (n_pops-1)/n_pops*s2 - (2*n_ave - 1)/(4*n_ave)*h_ave)
    var c=1/2*h_ave
    
    var fst = a/(a+b+c) 

    var res: FST_OUT
    res = (n_pops: n_pops, 
                n_ave: n_ave, 
                p_ave: p_ave, 
                h_ave: h_ave, 
                he: he, 
                n_c: n_c, 
                s2: s2,
                a: a,
                b: b,
                c: c, 
                fst: fst,
                )
    return res




proc formatFstOut(v: Variant, fst: FST_OUT, converged: bool, p: float64, good: int, precision: int=6): string =
    var lineout = "{v.CHROM}\t{v.POS}\t{fst.fst.formatFloat(ffDecimal, precision)}\t{fst.a.formatFloat(ffDecimal, precision)}\t{fst.b.formatFloat(ffDecimal, precision)}\t{fst.c.formatFloat(ffDecimal, precision)}\t{converged}".fmt
    if -99.9 < p:
        lineout = lineout & "\t{p.formatFloat(ffDecimal, precision)}\t{good}".fmt

    lineout = lineout & "\n"
    return lineout

proc formatMleAfOut(v: Variant, af: MLE_AF, pop:string, precision: int=6): string =
    #(n: int(n), af: mle_af[1], gf: mle_gf, converged: converged, iterations: iteration, starts: starts)
    let gtfs ="{af.gf[0].formatFloat(ffDecimal, precision)},{af.gf[1].formatFloat(ffDecimal, precision)},{af.gf[2].formatFloat(ffDecimal, precision)}".fmt
    let starts ="{af.starts[0].formatFloat(ffDecimal, precision)},{af.starts[1].formatFloat(ffDecimal, precision)},{af.starts[2].formatFloat(ffDecimal, precision)}".fmt
    var lineout = "{v.CHROM}\t{v.POS}\t{pop}\t{af.n}\t{af.ngood}\t{af.af.formatFloat(ffDecimal, precision)}\t{gtfs}\t{af.converged}\t{af.iterations}\t{starts}\n".fmt
    return lineout

proc boot_fst(glsr_rec: seq[seq[float]], v:Variant, pop_vec: seq[string], pop_labels: seq[string], sample_labels: seq[string], random_starts: bool, replace: bool): float =
    
    # Iterater over pops and estimate AFs
    var fst: float

    # storage for fst calculations
    var afs = newSeq[float64](len(pop_labels))
    var ns = newSeq[int](len(pop_labels))
    var pfreqs = newSeq[float64](len(pop_labels))
    var hfreqs = newSeq[float64](len(pop_labels))
    #echo "len(pop_rec): {len(pop_rec)}\tlen(glsr_rec): {len(glsr_rec)}".fmt
    ## estimation
    var pop_vec_boot = pop_vec
    # initalize rng and shuffle GLs
    var r = initRand()
    if replace:
        pop_vec_boot = pop_vec_boot.mapIt(r.sample(pop_vec))
    else:
        r.shuffle(pop_vec_boot)


    var convergence=true
    for ix, pop in enumerate(pop_labels):
        #echo "{pop} n={count(pop_rec, pop)}".fmt
        var pop_samps = zip(pop_vec_boot, sample_labels).filterIt(it[0]==pop).unzip()[1]
        var pop_glsr = zip(pop_vec_boot, glsr_rec).filterIt(it[0]==pop).unzip()[1]
        #echo "{pop}: n={len(pop_glsr)}".fmt
        var mle_af = maf_em(pop_glsr, tol=1e-8, max_iter=1000, random_starts=random_starts)

        ## print out some info if EM fails right away becuase of NaN/Inf
        if (mle_af.converged==false):
            convergence=false # sets to false if any of the pop AF EM fails to converge
            stderr.writeLine("Warning: MLE EM did not converge (bootstrapping) [{v}, pop: {pop}, n_iter: {mle_af.iterations}]".fmt)
            if (mle_af.iterations==1):
                stderr.writeLine("possible issue with Nan/Inf in GLs".fmt)
                for i in glsr_rec:
                    if gl_nan_check(i) == false:
                        stderr.writeLine("sample:{i}, unscaled PL/GL {i}, using GL?".fmt)

        #echo "\tESTIMATION - N: {mle_af.n}, AF: {mle_af.af:0.4f}, GF: {mle_af.gf[0]:0.4f},{mle_af.gf[1]:0.4f},{mle_af.gf[2]:0.4f}".fmt
        #afs.add(mle_af.af)
        ns[ix]=mle_af.n
        pfreqs[ix]=(mle_af.gf[1]/2) + mle_af.gf[2]
        hfreqs[ix]=mle_af.gf[1]

    ## estimate fst
    var fst_out=fst_wc(pfreqs, hfreqs, ns)
    if convergence:
        fst = fst_out.fst
    else:
        fst = -999.0
    
    return fst

proc processRecord(v: Variant, 
                   n:int, 
                   pop_vec: seq[string], 
                   pop_labels: seq[string], 
                   sample_labels: seq[string], 
                   remove_missing:bool=false, 
                   random_starts:bool=false, 
                   nboots:int=0, 
                   replace:bool=false, 
                   em_tol: float=1e-8, 
                   max_iter:int=200, 
                   verbosity: int=0): seq[string] {.gcsafe.} = 

    var pop_rec: seq[string]
    var samples_rec: seq[string]

    var pls = newSeq[int32](n*3)
    var plsr = newSeqWith(n, newSeq[int32](3))
    var glsr = newSeqWith(n, newSeq[float64](3))
    var glsr_rec: seq[seq[float64]]

    var dps = new_seq[int32](n)
    var dp_mask: seq[bool]
    dp_mask = dps.mapIt(true) ## Default for use all samples even if missing

    doAssert v.format.get("PL", pls) == Status.OK
    plsr = pls.reshape_fmt(n)

    try: glsr = plsr.mapIt(pl2prob(it, scale=true)) ## unscale to beagle gls
    except OverflowError: 
        var error = getCurrentException()
        if 0 < verbosity:  
            stderr.writeLine("OverflowError:\n{error.msg} exception occured when coverting PLs to probs at {v.CHROM}:{v.POS}...skipping this locus....".fmt)
            for (ix, pl) in enumerate(plsr):
                    if pl.mapIt(it.is_missing()).anyIt(it):
                        stderr.writeLine("Offending sample index: {ix}, PL:{pl}".fmt)
            return @["NULL", "NULL"] #glsr


    if remove_missing:
        doAssert v.format.get("DP", dps) == Status.OK
        dp_mask = dps.mapIt(0 < it)

        pop_rec = zip(dp_mask, pop_vec).filterIt(it[0]).unzip()[1]
        samples_rec = zip(dp_mask, sample_labels).filterIt(it[0]).unzip()[1]
        glsr_rec = zip(dp_mask, glsr).filterIt(it[0]).unzip()[1]
    else:
        pop_rec = pop_vec
        samples_rec = sample_labels
        glsr_rec = glsr
        
    ## FST Estimation
    ## Iterate over pops and estimate AFs
    # storage for fst calculations
    var ns = newSeq[int](len(pop_labels))
    var pfreqs = newSeq[float64](len(pop_labels))
    var hfreqs = newSeq[float64](len(pop_labels))
    var af_lines= ""
    
    var convergence=true
    for ix, pop in enumerate(pop_labels):
        #echo "{pop} n={count(pop_rec, pop)}".fmt
        var pop_samps = zip(pop_rec, samples_rec).filterIt(it[0]==pop).unzip()[1]
        var pop_glsr = zip(pop_rec, glsr_rec).filterIt(it[0]==pop).unzip()[1]
        ## !!! Check for missing data here? !!!
        #echo "{pop}: n={len(pop_glsr)}".fmt

        var mle_af = maf_em(pop_glsr, tol=em_tol, max_iter=max_iter, random_starts=random_starts)
        if mle_af.ngood == 0:
            echo "Warning: No data for {pop} at {v}. Setting EM estimates to NaN.".fmt
            mle_af.af = NaN;
            mle_af.gf = [NaN, NaN, NaN]
            


        af_lines = af_lines & formatMleAfOut(v, mle_af, pop, precision=6);

        ## print out some info if EM fails right away becuase of NaN/Inf
        if (mle_af.converged==false):
            convergence=false # sets to false if any of the pop AF EM fails to converge
            if 0 < verbosity:
                stderr.writeLine("Warning: MLE EM did not converge [{v}, pop: {pop}, n_iter: {mle_af.iterations}]".fmt)
                if (mle_af.iterations==1):
                    stderr.writeLine("possible issue with Nan/Inf in GLs".fmt)
                    for i in glsr_rec:
                        if gl_nan_check(i) == false:
                            stderr.writeLine("sample index:{i}, unscaled PL/GL {i}, using GL?".fmt)

        ns[ix]=mle_af.n
        pfreqs[ix]=(mle_af.gf[1]/2) + mle_af.gf[2]
        hfreqs[ix]=mle_af.gf[1]
                                         
    ## estimate fst
    var fst_out=fst_wc(pfreqs, hfreqs, ns)
    var fst: float64                            
    if convergence:
        fst = fst_out.fst
    else:
        fst = -999.0
    
    
    ## --- Bootstrapping --- 
    var good_boots = 0
    var p = NaN 
    var boots = newSeq[float](nboots) 
    if 0 < nboots:
        boots = boots.mapIt(-999.0)
        boots = boots.mapIt(boot_fst(glsr_rec, v, pop_rec, pop_labels, samples_rec, random_starts, replace=replace))
        good_boots = boots.countIt(-999.0!=it)
        p = boots.countIt(fst<=it)/good_boots


    return @[formatFstOut(v, fst_out,convergence, p, good_boots, precision=6), af_lines]
    #return @[fst, p]



proc wgsFst(bcf: string, 
             popfile: string,  
             out_prefix: string="out",
             query: string="*",  
             remove_missing: bool=false, 
             random_starts:bool=false, 
             nboots: int = 0, 
             replace:bool=false, 
             n_snps: int=50000, 
             em_tol: float=1e-8, 
             max_iter:int=200, 
             verbosity: int=1): void =


    
    let pop_map = parse_pop_file(popfile)

    var v:VCF
    doAssert(open(v, bcf))

    let strm_fst = newGzFileStream("{out_prefix}-fst.txt.gz".fmt, fmWrite)
    var header= "chrom\tpos\tfst\ta\tb\tc\tall_converged"
    if 0 < nboots:
        header = header & "\tp\tconverged_boots"
    strm_fst.writeLine(header)

    let strm_af = newGzFileStream("{out_prefix}-af.txt.gz".fmt, fmWrite)
    var headeraf= "chrom\tpos\tpop\tn\tngood\tmleaf\tmlegf\tconverged\tn_iter\tstarts"
    strm_af.writeLine(headeraf)

    ## --- Parse population file and set samples
    var all_samples: seq[string]            
    for key, value in pop_map:
        #echo "{key}: {value}".fmt
        all_samples &= key
    let 
        pops = deduplicate(all_samples.mapIt(pop_map[it]))
    
    
    if 0 < verbosity: 
        echo "POPFILE POPS: {pops}".fmt
        echo "{len(v.samples)} samples found in vcf".fmt

    var samples = v.samples
    var not_found = all_samples.filterIt(it.notin(samples))
    if 0 < len(not_found):
        if 0 < verbosity: 
            stderr.writeLine("WARNING: could not find {len(not_found)} samples listed in the pop file in {bcf}:".fmt)
            for i in not_found:
                stderr.writeLine("{i}".fmt)



    set_samples(v, all_samples)
    samples = v.samples

    let pop_vec = samples.mapit(pop_map[it])

    if 0 < verbosity: 
        echo "{len(samples)} samples in vcf after filtering".fmt
    
    ## ---- bcf processing ---

    var start_time =now()

    #var (output, excode)=execCmdEx("bcftools index --nrecords {bcf}".fmt);
    #let n_snps = parseInt(output.strip())
    #echo "found {n_snps} records in {now()-start_time}".fmt
    var query = query
    var query_obj: Query
    if query != "*":
        query_obj = parse_query_str(query, v)
        query = "{query_obj.chrom}:{query_obj.start}-{query_obj.stop}".fmt
        if 1 < verbosity:
            echo "parsed query {query_obj}='{query}'".fmt
            

    var chunk = 0
    var keep_reading = true

    var variants = newSeq[Variant](n_snps)
    while keep_reading:
        chunk += 1
        start_time =now()
        
        var ix = 0
        for rec in v.query(query):
            variants[ix] = rec.copy()
            ix += 1
            if (ix == n_snps):
                if query != "*":
                    query = "{query_obj.chrom}:{rec.POS+1}-{query_obj.stop}".fmt
                    if 1 < verbosity:
                        echo "CHUNK {chunk}: Ends at {rec}  [next region = {query}]".fmt
                break

        ## truncate incase we cant fill up enitre chunk 
        ## (eg. fewer variants than chunk size)
        if ix < len(variants):
            variants = variants[0..(ix-1)]
            keep_reading = false

        if 1 < verbosity: 
            echo "CHUNK {chunk}: {len(variants)} BCF/VCF records loaded into memory".fmt



        start_time = now()

        var fst_output = newSeq[seq[string]](len(variants))
        let bufOut = cast[ptr UncheckedArray[seq[string]]](fst_output[0].unsafeAddr)
        let bufIn = cast[ptr UncheckedArray[Variant]](variants[0].unsafeAddr)

        {.push stackTrace:off.}
        for i in 0||(len(variants)-1):
            bufOut[i] = processRecord(bufIn[i], len(samples), pop_vec, pops, samples, 
                                        remove_missing=remove_missing, random_starts=random_starts, 
                                        nboots=nboots, em_tol=em_tol, max_iter=max_iter, verbosity=verbosity)
        {.pop.}

        if 1 < verbosity: 
            echo "estimated fst for {len(variants)} SNPs in {now()-start_time}".fmt
        
        start_time =now()
        let chunk_size = 512*1000

        ## dump fst to file
        var text = fst_output.mapIt(it[0]).join("") 
        var num_bytes = text.len


        var idx = 0
        while true:
            strm_fst.writeData(text[idx].unsafeAddr, min(num_bytes, chunk_size))
            if num_bytes < chunk_size:
                break
            dec(num_bytes, chunk_size)
            inc(idx, chunk_size)
        

        ## Dump AF to file
        text = fst_output.mapIt(it[1]).join("") 
        num_bytes = text.len
        
        idx = 0
        while true:
            strm_af.writeData(text[idx].unsafeAddr, min(num_bytes, chunk_size))
            if num_bytes < chunk_size:
                break
            dec(num_bytes, chunk_size)
            inc(idx, chunk_size)

  
        if 1 < verbosity: 
            echo "file writing took {now()-start_time}\n".fmt
         
        
        ## for debug prposes, stop file reading after 3 chunks are processed
        #if chunk == 3:
        #    keep_reading = false
        GC_runOrc()
    strm_fst.close()
    strm_af.close()



when isMainModule:
    import therapist

    # The parser is specified as a tuple
    let spec = (
        # Name is a positional argument, by virtue of being surrounded by < and >
        bcf: newStringArg(@["<BCF>"], help="bcf/vcf.gz/vcf file to read"),
        # --times is an optional argument, by virtue of starting with - and/or --
        popfile: newStringArg(@["-p", "--pop-file"], help="File contining sample names and population labels"),
        out_prefix: newStringArg(@["-o", "--output"], help="Prefix for output files", defaultVal="out"),
        region: newStringArg(@["-r", "--region"], help="Specify a region to process [chr]:[start-stop].", defaultVal="*"),
        remove_missing: newFlagArg(@["-m", "--remove-missing"], help="Remove samples with missing data? Requires DP field in BCF. [default: false]"),
        random_starts: newFlagArg(@["-s", "--random-starts"], help="Use random starting values for EM algorithm? [default: false]"),
        nboots: newIntArg(@["-b", "--boots"], help="Number of bootstrap resamplings", defaultVal=0),
        replace: newFlagArg(@["--replace"], help="Perform re-sampling with replacement? [default: false]"),
        threads: newIntArg(@["-t", "--threads"], help="Number of threads to use in parallel sections of code.", defaultVal=1),
        n_snps: newIntArg(@["-c", "--chunk-size"], help="Number of SNPs to read into memory at a time.", defaultVal=50000),
        em_tol: newFloatArg(@["--em-tol"], help="Tolerance required for stopping EM estimation of allele frequencies.", defaultVal=1e-8),
        max_iter: newIntArg(@["--max-iter"], help="Maximum iterations allowed for EM estimation of allele frequencies.", defaultVal=200),
        verbosity: newCountArg(@["-v", "--verbose"], help="Verbose. Can be specified multiple times for more detailed output."),
        # --version will cause 0.1.0 to be printed
        version: newMessageArg(@["--version"], "0.0.1", help="Prints version"),
        # --help will cause a help message to be printed
        help: newHelpArg(@["-h", "--help"], help="Show help message"),
    )

    # `args` and `command` would normally be picked up from the commandline
    spec.parseOrQuit(prolog = "fst-gl - a program to estimate FST from genotype likelihoods.")
    # If a help message or version was requested or a parse error generated it would be printed
    # and then the parser would call `quit`. Getting past `parseOrQuit` implies we're ok.
    var remove_missing = false
    var random_starts = false
    var replace = false

    if 0 < spec.remove_missing.count:
        remove_missing = true

    if 0 < spec.random_starts.count:
        random_starts = true

    if 0 < spec.replace.count:
        replace = true


    echo "Processing BCF {spec.bcf.value} with the following options:".fmt
    echo "\t--pop-file={spec.popfile.value}".fmt
    echo "\t--output={spec.out_prefix.value}".fmt
    echo "\t--region={spec.region.value}".fmt
    echo "\t--remove-missing={spec.remove_missing.count} ({remove_missing})".fmt
    echo "\t--random-starts={spec.random_starts.count} ({random_starts})".fmt
    echo "\t--boots={spec.nboots.value}".fmt
    echo "\t--replace={spec.replace.count} ({replace})".fmt
    echo "\t--threads={spec.threads.value}".fmt
    echo "\t--chunk-size={spec.n_snps.value}".fmt
    echo "\t--em-tol={spec.em_tol.value}".fmt
    echo "\t--max-iter={spec.max_iter.value}".fmt
    echo "\t--verbose={spec.verbosity.count}".fmt


    ## set number of threads
    omp_set_num_threads(cint(spec.threads.value))
    if 1 < spec.verbosity.count:
        echo "OMP_MAX_THREADS: {omp_get_max_threads()}".fmt

    wgsFst(bcf=spec.bcf.value, 
            popfile=spec.popfile.value,  
            out_prefix=spec.out_prefix.value,
            query=spec.region.value,
            remove_missing=remove_missing, 
            random_starts=random_starts, 
            nboots=spec.nboots.value, 
            replace=replace, 
            n_snps=spec.n_snps.value, 
            em_tol=spec.em_tol.value,
            max_iter=spec.max_iter.value, 
            verbosity=spec.verbosity.count)
