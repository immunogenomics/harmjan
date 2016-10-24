package nl.harmjanwestra.miscscripts;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.Correlation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 10/20/16.
 */
public class GenomeAutoCorrelation {

	public static void main(String[] args) {

		GenomeAutoCorrelation c = new GenomeAutoCorrelation();
		Feature region = new Feature(Chromosome.TWO, 190250000, 193475000);
		try {
			c.run("/Data/Ref/Ensembl/Genome/Homo_sapiens.GRCh37.75.dna.chromosome.2.fa.gz", "", region, 20);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void run(String filein, String fileout, Feature region, int resolution) throws IOException {

		FastaSequenceFile f = new FastaSequenceFile(new File(filein), true);

		ReferenceSequence seq = f.nextSequence();

		byte[] bases = seq.getBases();

		// 190250000-193475000

		int width = region.getStop() - region.getStart();


		int nrwindows = width / resolution;
		int bpperwindow = (int) Math.floor((double) width / nrwindows);
		ArrayList<SeqWindow> results = new ArrayList<SeqWindow>();

		int zzz = 0;
		for (int b = region.getStart(); b < region.getStop(); b += bpperwindow) {
			int windowstart = b - (resolution / 2);
			int windowstop = b + (resolution / 2);
			double[] seq1 = getseq(bases, windowstart, windowstop);
			for (int b2 = b + (bpperwindow / 2); b2 < region.getStop(); b2 += bpperwindow) {
				int window2start = b2 - (resolution / 2);
				int window2stop = b2 + (resolution / 2);
				double[] seq2 = getseq(bases, window2start, window2stop);

				Pair<double[], double[]> remainingbases = filterregions(seq1, seq2);
				double[] d1 = remainingbases.getLeft();
				double[] d2 = remainingbases.getRight();
				double r = Correlation.correlate(d1, d2);
				System.out.println(windowstart + "\t" + window2start + "\t" + r + "\t" + d1.length);
				zzz++;
				if (zzz == 10000) {
					System.exit(-1);
				}
			}
		}

	}

	private Pair<double[], double[]> filterregions(double[] seq1, double[] seq2) {

		int nrremain = 0;
		for (int b = 0; b < seq1.length; b++) {
			if (BaseAnnot.isDNABase(seq1[b]) && BaseAnnot.isDNABase(seq2[b])) {
				nrremain++;
			}
		}

		double[] seq1out = new double[nrremain];
		double[] seq2out = new double[nrremain];
		int q = 0;

		for (int b = 0; b < seq1.length; b++) {
			if (BaseAnnot.isDNABase(seq1[b]) && BaseAnnot.isDNABase(seq2[b])) {
				seq1out[q] = seq1[b];
				seq2out[q] = seq2[b];

				q++;
			}
		}

		return new Pair<double[], double[]>(seq1out, seq2out);
	}

	private double[] getseq(byte[] bases, int windowstart, int windowstop) {
		double[] out = new double[windowstop - windowstart];
		int z = 0;
		for (int q = windowstart; q < windowstop; q++) {
			out[z] = bases[q];
			z++;
		}
		return out;
	}
}

class SeqWindow {

	public int r1;
	public int r2;
	public double d;

}
