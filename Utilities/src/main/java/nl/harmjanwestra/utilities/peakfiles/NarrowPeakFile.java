package nl.harmjanwestra.utilities.peakfiles;

/**
 * Created by hwestra on 9/1/15.
 */
public class NarrowPeakFile {

	/*
	This format is used to provide called peaks of signal enrichment based on pooled, normalized (interpreted) data. It is a BED6+4 format.

    chrom - Name of the chromosome (or contig, scaffold, etc.).
    chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
    chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
    name - Name given to a region (preferably unique). Use '.' if no name is assigned.
    score - Indicates how dark the peak will be displayed in the browser (0-1000). If all scores were '0' when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
    strand - +/- to denote strand or orientation (whenever applicable). Use '.' if no orientation is assigned.
    signalValue - Measurement of overall (usually, average) enrichment for the region.
    pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
    qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
    peak - Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.

Here is an example of narrowPeak format:

track type=narrowPeak visibility=3 db=hg19 name="nPk" description="ENCODE narrowPeak Example"
browser position chr1:9356000-9365000
chr1    9356548 9356648 .       0       .       182     5.0945  -1  50
chr1    9358722 9358822 .       0       .       91      4.6052  -1  40
chr1    9361082 9361182 .       0       .       182     9.2103  -1  75

	 */
}
