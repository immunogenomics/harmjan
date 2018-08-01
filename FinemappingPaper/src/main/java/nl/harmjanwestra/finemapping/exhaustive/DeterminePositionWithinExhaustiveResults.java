package nl.harmjanwestra.finemapping.exhaustive;

import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

/**
 * Created by Harm-Jan on 06/16/16.
 */
public class DeterminePositionWithinExhaustiveResults {

	public static void main(String[] args) {
		String regionfile = "";
		String assocfileprefix = "";
		String modelfileprefix = "";
		String outputfile = "";

		DeterminePositionWithinExhaustiveResults d = new DeterminePositionWithinExhaustiveResults();
		try {

//			regionfile = "D:\\Cloud\\Dropbox\\2016-03-RAT1D-Finemappng\\Data\\2016-07-25-SummaryStats\\Conditional\\RA-lociWithIndependentFX.txt";
//			assocfileprefix = "D:\\tmp\\2016-08-14-exhaustive\\RA-assoc0.3-COSMO-chr";
//			modelfileprefix = "D:\\tmp\\2016-08-14-exhaustive\\models\\RA-assoc0.3-COSMO-chr";
//			outputfile = "D:\\tmp\\2016-08-14-exhaustive\\output\\RA-exhaustiveout-PositionOfTopConditionalEffects.txt";
//			d.mergeCredibleSets(regionfile, assocfileprefix, modelfileprefix, outputfile);
//
//			regionfile = "D:\\Cloud\\Dropbox\\2016-03-RAT1D-Finemappng\\Data\\2016-07-25-SummaryStats\\Conditional\\T1D-lociWithIndependentFX.txt";
//			assocfileprefix = "D:\\tmp\\2016-08-14-exhaustive\\T1D-assoc0.3-COSMO-chr";
//			modelfileprefix = "D:\\tmp\\2016-08-14-exhaustive\\models\\T1D-assoc0.3-COSMO-chr";
//			outputfile = "D:\\tmp\\2016-08-14-exhaustive\\output\\T1D-exhaustive-PositionOfTopConditionalEffects.txt";
//			d.mergeCredibleSets(regionfile, assocfileprefix, modelfileprefix, outputfile);

//			regionfile = "/Data/tmp/tnfaip3/tnfaip3.txt";
//			assocfileprefix = "/Data/tmp/tnfaip3/RA-assoc0.3-COSMO-chr";
//			modelfileprefix = "/Data/tmp/tnfaip3/RA-assoc0.3-COSMO-chr";
//			outputfile = "/Data/tmp/tnfaip3/exhaustivecomp.txt";
//			d.compare(regionfile, assocfileprefix, modelfileprefix, outputfile);


			// CD28
//			regionfile = "/Data/Projects/2016-Finemapping/Exhaustive/cd28/cd28.bed";
//			assocfileprefix = "/Data/Projects/2016-Finemapping/Exhaustive/RA-assoc0.3-COSMO-chr";
//			modelfileprefix = "/Data/Projects/2016-Finemapping/Exhaustive/RA-assoc0.3-COSMO-chr";
//			outputfile = "/Data/Projects/2016-Finemapping/Exhaustive/cd28/exhaustivecompRA.txt";
//			d.compare(regionfile, assocfileprefix, modelfileprefix, outputfile);
//
//			regionfile = "/Data/Projects/2016-Finemapping/Exhaustive/cd28/cd28.bed";
//			assocfileprefix = "/Data/Projects/2016-Finemapping/Exhaustive/T1D-assoc0.3-COSMO-chr";
//			modelfileprefix = "/Data/Projects/2016-Finemapping/Exhaustive/RA-assoc0.3-COSMO-chr";
//			outputfile = "/Data/Projects/2016-Finemapping/Exhaustive/cd28/exhaustivecompT1D.txt";
//			d.compare(regionfile, assocfileprefix, modelfileprefix, outputfile);


			regionfile = "/Data/Projects/2016-Finemapping/Exhaustive/tnfaip3.bed";
			assocfileprefix = "/Data/Projects/2016-Finemapping/Exhaustive/RA-assoc0.3-COSMO-chr";
			modelfileprefix = "/Data/Projects/2016-Finemapping/Exhaustive/RA-assoc0.3-COSMO-chr";
			outputfile = "/Data/Projects/2016-Finemapping/Exhaustive/tnfaip3/exhaustivecompRA.txt";
			d.run(regionfile, assocfileprefix, modelfileprefix, outputfile);


//			regionfile = "/Data/tmp/2016-06-20/TNFAIP3.bed";
//			assocfileprefix = "/Data/tmp/2016-06-20/RA-assoc0.3-COSMO-TNFAIP3-chr";
//			modelfileprefix = "/Data/tmp/2016-06-20/RA-assoc0.3-COSMO-chr";
//			outputfile = "/Data/tmp/2016-06-20/RA-TNFAIP3";
//			d.determineRegionSignificanceThresholds(regionfile, assocfileprefix, modelfileprefix, outputfile);

		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void run(String regionfile,
					String assocfileprefix,
					String modelfileprefix,
					String outputfile) throws IOException {


		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(regionfile);


		TextFile outf = new TextFile(outputfile, TextFile.W);
		outf.writeln("Region\tConditionalVariant1\tConditionalVariant2\tP\tPresentInExhaustiveOutput\tNrLowerPvals\tNrTotalPvals\tExhaustiveTopVariant1\tExhaustiveTopVariant2\tTopP");
		for (int r = 0; r < regions.size(); r++) {
			Feature region = regions.get(r);
			String modelfile = modelfileprefix + region.getChromosome().getNumber() + "-gwas-conditional-models.txt";
			String assocfile = assocfileprefix + region.getChromosome().getNumber() + "-pairwise.txt.gz";
			String[] variants = loadConditionalVariants(region, modelfile);

			if (variants == null) {
				System.err.println("Model path is broken for some reason: " + modelfile);
				System.exit(-1);
			}


			Double pvalForCombo = 0d;
			ArrayList<Double> allPvals = new ArrayList<>();

			String regionStr = region.toString();

			String[] topEffect = null;
			Double maxP = 0d;
			System.out.println("Looking for variants: " + variants[0] + "\t" + variants[1]);
			boolean foundit = false;

			Pair<Double, String> lowestPair = null;

			ArrayList<Pair<Double, String>> outputBuffer = null;
			double highestPInBuffer = 0;
			int numberOfTopFx = 2000;
			ArrayList<Pair<Double, String>> workBuffer = new ArrayList<Pair<Double, String>>(numberOfTopFx);

			boolean lineIsWhatWereLookingFor = false;
			// get all pvalues in region
			TextFile tf2 = new TextFile(assocfile, TextFile.R);
			String headerln = tf2.readLine();
			String[] elems = tf2.readLineElems(TextFile.tab);
			while (elems != null) {
				lineIsWhatWereLookingFor = false;
				Feature s1 = new Feature();
				Chromosome chr = Chromosome.parseChr(elems[0]);
				Integer pos1 = Integer.parseInt(elems[1]);
				String snp1Id = elems[2];
				Integer pos2 = Integer.parseInt(elems[5]);
				String snp2Id = elems[6];

				SNPFeature snp1 = new SNPFeature();
				snp1.setChromosome(chr);
				snp1.setStart(pos1);
				snp1.setStop(pos1);
				snp1.setName(snp1Id);

				SNPFeature snp2 = new SNPFeature();
				snp2.setChromosome(chr);
				snp2.setStart(pos2);
				snp2.setStop(pos2);
				snp2.setName(snp2Id);


				if (region.overlaps(snp1) && region.overlaps(snp2)) {

					String snp1str = snp1.getChromosome().getNumber() + "_" + snp1.getStart() + "_" + snp1.getName();
					String snp2str = snp2.getChromosome().getNumber() + "_" + snp2.getStart() + "_" + snp2.getName();
//					System.out.println(snp1str);
//					System.exit(0);
					Double p = Double.parseDouble(elems[elems.length - 1]);

					if ((snp1str.equals(variants[0]) && snp2str.equals(variants[1]))
							|| (snp1str.equals(variants[1]) && snp2str.equals(variants[0]))) {
						System.out.println("found it.");
						pvalForCombo = p;
						foundit = true;
						lineIsWhatWereLookingFor = true;
					}

					if (p > maxP) {
						topEffect = elems;
						maxP = p;
					}
					if (!lineIsWhatWereLookingFor) {
						allPvals.add(p);
					}

					if (!Double.isNaN(p) && !Double.isInfinite(p)) {
						Pair<Double, String> pair = new Pair<Double, String>(p, Strings.concat(elems, Strings.tab), Pair.SORTBY.LEFT);
						if (lowestPair == null) {
							lowestPair = pair;
						} else {
							if (p > lowestPair.getLeft()) {
								lowestPair = pair;
							}
						}
						workBuffer.add(pair);
						if (workBuffer.size() == numberOfTopFx) {
							if (outputBuffer == null) {
								outputBuffer = workBuffer;
								Collections.sort(outputBuffer);
								System.out.println("Set outputbuffer: " + outputBuffer.size());
							} else {
//								System.out.println("Update outputbuffer");
								outputBuffer.addAll(workBuffer);
								Collections.sort(outputBuffer, new PairSorter());
								ArrayList<Pair<Double, String>> tmp = new ArrayList<>(numberOfTopFx);
								tmp.addAll(outputBuffer.subList(0, numberOfTopFx));
								outputBuffer = tmp;
							}
							workBuffer = new ArrayList<>(numberOfTopFx);
						}
					}

				}
				elems = tf2.readLineElems(TextFile.tab);
			}
			tf2.close();

			if (outputBuffer == null) {
				outputBuffer = workBuffer;
				Collections.sort(outputBuffer);
			} else {
				outputBuffer.addAll(workBuffer);
				Collections.sort(outputBuffer);
				ArrayList<Pair<Double, String>> tmp = new ArrayList<>(numberOfTopFx);
				tmp.addAll(outputBuffer.subList(0, numberOfTopFx));
				outputBuffer = tmp;
			}


			TextFile outtop = new TextFile(outputfile + "_" + region.toString() + "-top" + numberOfTopFx + ".txt.gz", TextFile.W);
			outtop.writeln(headerln);
			for (int d = 0; d < outputBuffer.size(); d++) {
				outtop.writeln(outputBuffer.get(d).getRight());
			}
			outtop.close();

			int nrLowerPvals = 0;
			for (Double d : allPvals) {
				if (d >= pvalForCombo) {
					nrLowerPvals++;
				}
			}


			Chromosome chr = Chromosome.parseChr(topEffect[0]);
			Integer pos1 = Integer.parseInt(topEffect[1]);
			String snp1Id = topEffect[2];
			Integer pos2 = Integer.parseInt(topEffect[5]);
			String snp2Id = topEffect[6];
			SNPFeature snp1 = new SNPFeature();
			snp1.setChromosome(chr);
			snp1.setStart(pos1);
			snp1.setStop(pos1);
			snp1.setName(snp1Id);

			SNPFeature snp2 = new SNPFeature();
			snp2.setChromosome(chr);
			snp2.setStart(pos2);
			snp2.setStop(pos2);
			snp2.setName(snp2Id);

			String output = regionStr + "\t" +
					variants[0] + "\t" +
					variants[1] + "\t" +
					pvalForCombo + "\t" +
					foundit + "\t" +
					nrLowerPvals + "\t" +
					allPvals.size() + "\t" +
					snp1.toString() + "\t" +
					snp2.toString() + "\t" +
					maxP;

			System.out.println(region.toString() + "\nConditional effect: " + (nrLowerPvals) + " lower pvals out of " + allPvals.size());
			System.out.println("lowest pair: " + lowestPair.getLeft() + "\t" + lowestPair.getRight());

			double max = Collections.max(allPvals);
			double min = Collections.min(allPvals);
			System.out.println();
			System.out.println("Pval of combo: " + pvalForCombo);
			System.out.println("Max: " + max);
			System.out.println("Min: " + min);

			int nrbins = 100;
			int[] bins = new int[nrbins];
			Range range = new Range(min, min, max, max);
			for (double d : allPvals) {
				double pos = range.getRelativePositionX(d);
				int bin = (int) Math.floor(pos * nrbins);
				if (bin >= nrbins) {
					bin = nrbins - 1;
				}
				if (bin < 0) {
					bin = 0;
				}
				bins[bin]++;
			}


			outf.writeln(output);

			TextFile binout = new TextFile(outputfile + "_" + region.toString() + "-bins.txt", TextFile.W);
			binout.writeln("PvalBin\tCt\tPerc");
			double stepsize = range.getRangeX() / nrbins;
			for (int i = 0; i < nrbins; i++) {
				binout.writeln(
						(min + (stepsize * i))
								+ "\t" + bins[i]
								+ "\t" + ((double) bins[i] / allPvals.size())
				);
			}
			binout.close();
		}

		outf.close();

	}

	private String[] loadConditionalVariants(Feature region, String modelfile) throws IOException {

		TextFile tf = new TextFile(modelfile, TextFile.R);
		String ln = tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String modelregion = elems[0];
			if (region.toString().equals(modelregion)) {
				String iter = elems[1];
				if (iter.equals("2")) {
					String snps = elems[2];
					return snps.split(";");
				}

			}
			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();

		return null;
	}

	private class PairSorter implements Comparator<Pair<Double, String>> {
		@Override
		public int compare(Pair<Double, String> o1, Pair<Double, String> o2) {
			if (o1.getLeft() > o2.getLeft()) {
				return 1;
			} else if (o1.getLeft() < o2.getLeft()) {
				return -1;
			} else {
				return 0;
			}
		}
	}
}
