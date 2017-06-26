# ProxyFinder
This is a JAVA 8 program to calculate LD proxies for a set of of input variants, and is partly built using the [HTSJDK library](https://github.com/samtools/htsjdk). It provides possibilities to search for proxies, given arbitrary genomic distances, LD and MAF cutoffs, but also allows for pairwise calculations. It uses Tabix to *quickly* grep SNPs out of large reference dataset VCFs (such as 1000 genomes). It's nevertheless probably slow as hell.

[Download binaries here]()

## Requirements
This program will run on anything that can run Java SE 8. It is a command line program. I do recommend ample RAM though (>4Gb), although it's hard to say how much, because Java. I also recommend running this software on a 64-bit machine, with a 64-bit Java 8 runtime environment.

# General usage notes
ProxyFinder has three modes:
* --proxy - a mode to find proxy variants given a list of snps
* --pairwise - a mode to calculate LD between pairs of variants
* --locusld - a mode to calculate LD between all variants within a region

**This program currently cannot properly handle Chromosome X and Y.**

### Output format
Output is plain text. You want it gzipped? Just append .gz to the output filename. Generally the output will look something like this:
</pre>
ChromA	PosA	RsIdA	ChromB	PosB	RsIdB	Distance	RSquared	Dprime
Chr2	202149589	rs1045485	Chr2	202149589	rs1045485	0	0.9999999999999981	0.999999999999999
Chr2	202149589	rs1045485	Chr2	202150914	rs35550815	1325	0.9999999999999981	0.999999999999999
</pre>

### Input genotypes and Tabix chromosome template
ProxyFinder relies on user supplied VCF variant files. These may be indexed with Tabix (e.g. [1000 genomes](http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/individual_chromosomes/)), but may also be unsorted VCF files. Please note however, that using non-indexed or unsorted VCF files will seriously harm (heh, that's part of my name) performance. You have been warned.
Some people like to store their Tabix VCF genotypes in a single file, others split everything in multiple chromosome files. I feel your pain, so I enable you to do both. If your files are split per chromosome, e.g. ``1kg.phase3.v5a.chr1.vcf.gz``, simply replace the chromosome number with the text CHR, e.g.: ``1kg.phase3.v5a.chrCHR.vcf.gz``. Magic will ensue. 

### Calculating LD for a subset of individuals
You want to exclude some people from the party? Make a list of individual IDs, matching the IDs of your VCF header, store them in a text file (can be gzipped), and specify ``--samplefilter /path/to/your/list.txt``. One individual per line. e.g.:
<pre>
sample1
sample2
sample3
</pre>

### Threads
Computers have many cores nowadays. Make use of them with ``--threads n``. Set n to the number of CPUs in your machine, or lower, depending on your sysadmins generosity.

### Matching on rsID
ProxyFinder matches your input against the VCF files you provide. When it can't find one of the specified variants, it will return the variant as a proxy for itself. By default, ProxyFinder will match on chromosome+location. It does not care about alleles, and rsID. If you do want to make it care about rsIds (you know, because it's oldschool), use the ```--matchrsid flag```. Hopefully, ProxyFinder will still be able to find your SNPs. 

# Proxy mode
Finds proxies for a given set of SNPs and other parameters. Example:
``java -Xmx4g -jar ./ProxyFinder.jar --proxy -i ./snps.txt -o ./proxies.txt.gz --tabix 1kg.phase3.v5a.chrCHR.vcf.gz``

Specify an LD output cutoff (e.g. r2>0.8) using:  ``--threshold 0.8``

Specify the width around the input SNP, where we should be looking for proxies (e.g. 1000 bp): ``--windowsize 1000`` 

Specify the minimal MAF (e.g. 1%) for proxies: ``--maf 0.01``
 
When using an unsorted VCF file, use: ``--vcf path.to.vcf.gz`` in stead of ``--tabix 1kg.phase3.v5a.chrCHR.vcf.gz``

#### Input file format
The ``snps.txt`` file in the above example, can have several formats. It expects at least three columns (only the first three are read). You can use chr\tpos\tid or id\tchr\tpos (which is used in [GoShifter](https://github.com/immunogenomics/goshifter)). The former format doesn't have a header, but the latter does (for GoShifter compatibility reasons):
<pre>
chr5	1279790	rs10069690
chr2	202149589	rs1045485
</pre>

or

<pre>
SNP	Chrom	BP
rs10069690	chr5	1279790
rs1045485	chr2	202149589
</pre>

# Pairwise mode
Calculate LD between pairs of SNPs. Example: 
``java -Xmx4g -jar ./ProxyFinder.jar --pairwise -i ./snps.txt -o ./proxies.txt.gz --tabix 1kg.phase3.v5a.chrCHR.vcf.gz``
Wow. That was familiar! No need to specify windows, thresholds or anything like that. You specify what you want to calculate. If both SNPs are present in the reference, you'll get a value.

When using an unsorted VCF file, use: ``--vcf path.to.vcf.gz`` in stead of ``--tabix 1kg.phase3.v5a.chrCHR.vcf.gz``

#### Input file format
As input, you can use either the three column format specified above (note, not the GoShifter format), and ProxyFinder will make all possible combinations (fun!). Or, you can specify the pairs you want to calculate, using a six column input:
<pre>
chr5	1279790	rs10069690	rs10941679	chr5	44706498
chr2	202149589	rs1045485	rs16857609	chr2	218296508
</pre>

# LocusLD
You want to know the LD within a region? Look no further. Example:
``java -Xmx4g -jar ./ProxyFinder.jar --locusld -r ./regions.bed -o ./proxydir/ --tabix 1kg.phase3.v5a.chrCHR.vcf.gz``
Note that here, in stead of ```-i /snps.txt```, we're using ```-r regions.bed```. 

Please note: these calculations can take quite a while, depending on the size of the regions. Also, this mode can use quite a lot of memory, since all variants are loaded into memory.

This mode will write an output file for each region. Specify the output directory with ``-o /outdir/``. Depending on the number of variants in the region, this may use quite some disk space, and are therefore gzipped by default.

#### Input file format
For specifying regions, use the BED3 format:
<pre>
chr1	1000	1001000
chr3	400001	400003
</pre>

# All command line parameters
<pre>
usage:
 -i,--snps <arg>           SNP path (format: 3 or 6 columns, tab
                           separated, one or two snps per line: chr pos
                           rsid)
    --locusld              Perform Pairwise LD calculation within a region
                           (provide regions with --regions)
 -m,--maf <arg>            MAF threshold [default: 0.005]
    --matchrsid            Match variants on RS id
 -o,--out <arg>            Output path
    --pairwise             Perform Pairwise LD calculation (use 6 column
                           file for --snps)
    --proxy
 -r,--regions <arg>        Region bed file path
    --samplefilter <arg>   Limit samples to individuals in this list (one
                           sample per line)
 -t,--threshold <arg>      R-squared threshold [default: 0.8]
    --tabix <arg>          Prefix for tabix path [format
                           /path/to/chrCHR.vcf.gz]. Replace the chromosome
                           number with CHR (will be replaced by chr number
                           depending on input SNP or region).
    --threads <arg>        Nr of threads [default: 1]
    --vcf <arg>            Use non-indexed VCF as input
 -w,--windowsize <arg>     Window size [default 1000000]
</pre>







