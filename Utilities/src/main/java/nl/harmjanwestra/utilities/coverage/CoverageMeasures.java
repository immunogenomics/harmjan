package nl.harmjanwestra.utilities.coverage;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.filter.AggregateFilter;
import nl.harmjanwestra.utilities.bamfile.BamFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by hwestra on 6/3/15.
 */
public class CoverageMeasures {


	public double calcRPKM(int nrReadsMappedInRegion, long totalReadsMappedOnGenome, int regionsize) {
		// RPKM = (10^9 * C)/(N * L)
		double rpkm = (Math.pow(10, 9) * nrReadsMappedInRegion) / ((long) regionsize * totalReadsMappedOnGenome);
		return rpkm;
	}

	// determine RPKM for a given region
	public Pair<Integer, Integer> countReadFeatures(Feature f, SAMRecordIterator it, AggregateFilter filter) {

		int start = f.getStart();
		int stop = f.getStop();
		int featuresize = Math.abs(stop - start);
		int nrReadsMappedInRegion = 0;
		int nrFragmentsMappedInRegion = 0;

		while (it.hasNext()) {
			SAMRecord record = it.next();

			if (filter == null || !filter.filterOut(record)) {

				if (record.getFirstOfPairFlag()) {
					Feature seqfeature = getSeqFeature(record, f);
					if (seqfeature.overlaps(f)) {
						nrFragmentsMappedInRegion++;
					}
				}
				nrReadsMappedInRegion++;

			}
		}

		return new Pair<Integer, Integer>(nrReadsMappedInRegion, nrFragmentsMappedInRegion);
	}

	public Pair<double[], int[]> determineRPKM(ArrayList<Feature> regions, String bamfile, AggregateFilter filter) throws Exception {

		BamFileReader reader = new BamFileReader(new File(bamfile));
		HashMap<String, String> chromosomeToSequence = reader.matchChromosomeNames(regions);

		int distributionsize = 1000;
		BasicStats stats = new BasicStats(reader, filter, distributionsize);

		int nrTotalReadsMapped = 0;


		int[] reads = new int[regions.size()];
		int[] fragments = new int[regions.size()];
		double[] rpkms = new double[regions.size()];
		double[] fpkms = new double[regions.size()];

		Pair<double[], int[]> output = new Pair<double[], int[]>(rpkms, reads);

		for (int i = 0; i < regions.size(); i++) {
			Feature f = regions.get(i);

			Chromosome c = f.getChromosome();
			int start = f.getStart();
			int stop = f.getStop();
			int regionSize = stop - start;
			String chrName = chromosomeToSequence.get(c.getName());

			if (chrName != null) {
				SAMRecordIterator it = reader.query(chrName, start, stop, true);

				Pair<Integer, Integer> readstats = countReadFeatures(f, it, filter);
				int nrReadsMappedToRegion = readstats.getLeft();
				int fragmentsForRegion = readstats.getRight();
				reads[i] = nrReadsMappedToRegion;
				fragments[i] = fragmentsForRegion;
				rpkms[i] = calcRPKM(nrReadsMappedToRegion, nrTotalReadsMapped, regionSize);

				it.close();
			}
		}

		reader.close();
		return output;
	}

	private Feature getSeqFeature(SAMRecord record, Feature region) {
		int start = record.getStart();
		int end = record.getMateAlignmentStart() + record.getReadLength();
		Feature f = new Feature();
		f.setChromosome(region.getChromosome());
		f.setStart(start);
		f.setStop(end);
		return f;
	}

	public Pair<double[][], ArrayList<String>> convertToTPM(double[][] dataMatrix, ArrayList<String> peakNames, double[] fraglen) {

		// filter out peaks smaller than fraglen


		boolean[] includePeak = new boolean[dataMatrix.length];
		int nrPeaksRemain = 0;
		ArrayList<String> newPeakList = new ArrayList<String>();
		double[][] effLen = new double[dataMatrix.length][dataMatrix[0].length];
		for (int i = 0; i < includePeak.length; i++) {
			String name = peakNames.get(i);
			String[] elems1 = name.split(":");
			String[] elems2 = elems1[1].split("-");

			int start = Integer.parseInt(elems2[0]);
			int end = Integer.parseInt(elems2[1]);
			int len = end - start;
			int ct = 0;
			for (int j = 0; j < fraglen.length; j++) {
				if (fraglen[j] > len) {
					ct++;
				}
				effLen[i][j] = len - fraglen[j] + 1;
			}
			if (ct > 0) {
				includePeak[i] = false;
			} else {
				includePeak[i] = true;
				nrPeaksRemain++;
				newPeakList.add(name);
			}

		}
		System.out.println(newPeakList.size() + " features remain for TPM");

		double[][] tpm = new double[nrPeaksRemain][fraglen.length];

		for (int col = 0; col < dataMatrix[0].length; col++) {
			double[] rate = new double[nrPeaksRemain];
			double ratesum = 0;
			int peakCtr = 0;
			for (int row = 0; row < dataMatrix.length; row++) {
				if (includePeak[row]) {
					rate[peakCtr] = Math.log(dataMatrix[row][col]) - Math.log(effLen[row][col]);
					ratesum += Math.exp(rate[peakCtr]);
					peakCtr++;
				}
			}

			double denom = Math.log(ratesum);
			for (int row = 0; row < rate.length; row++) {
				tpm[row][col] = Math.exp(rate[row] - denom + Math.log(1E6));
			}
		}

		return new Pair<double[][], ArrayList<String>>(tpm, newPeakList);

	}


	/*


	xi == counts or nr reads
	liexp = effective length (of region)
	liexp = li - mufld +1 ==> expregionlen = regionlen - mufraglen + 1
	effCounts = xi * li/liexp
	cpm = xi / N * 10^6  ==> N == number of fragments sequenced

	// within sample
	tpm = xi / liexp * 1/sum(xj/liexpj) * 10E6
	fpkm = xi / liexp * N * 10E9





	 */

}
