package nl.harmjanwestra.finemapping.genotypesandimputation;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 6/30/16.
 */
public class MergeAccuracyAndInfoScoreFiles extends VariantCounter {


	public static void main(String[] args) {


		boolean windows = false;
		String seqpanelvcf = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/seqpanel/seqpanelfiltered-maf0005-cr0950-rd10-gq30-runNamesFixed-RASampleNamesFixed-badSamplesRemoved-mixupsFixed.vcf.gz";
		String variantsOnIC = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/RAAndT1D-recode-maf0005-ICRegionsW100kb-samplenamefix.vcf.gz-updatedRSId-stats.vcf.gz";

		String[] files = new String[]{
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-08-31-Accuracy/T1D-EUR.txt",
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-08-31-Accuracy/T1D-COSMO.txt",
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-08-31-Accuracy/T1D-HRC-EAGLE.txt",
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-08-31-Accuracy/T1D-HRC-SHAPEIT.txt",
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-08-31-Accuracy/T1D-HRC-EAGLE-Michigan.txt"
		};
		String[] labels = new String[]{
				"EUR",
				"COSMO",
				"HRC / HRC / EAGLE",
				"HRC / HRC / SHAPEIT",
				"HRC / HRC / EAGLE / MICHIGAN"
		};
		String diseaseprefix = "T1D";
		String samplelist = null;
//
//		files = new String[]{
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-08-31-Accuracy/RA-EUR.txt",
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-08-31-Accuracy/RA-COSMO.txt",
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-08-31-Accuracy/RA-HRC-w100kb.txt",
//		};
//		labels = new String[]{
//				"EUR",
//				"COSMO",
//				"HRC / COSMO / EAGLE"
//		};
//		diseaseprefix = "RA";

		seqpanelvcf = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/seqpanel/seqpanelfiltered-maf0005-cr0950-rd10-gq30-runNamesFixed-RASampleNamesFixed-badSamplesRemoved-mixupsFixed.vcf.gz-updatedRSId.vcf.gz";
		variantsOnIC = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/RAAndT1D-recode-maf0005-ICRegionsW100kb-samplenamefix.vcf.gz-updatedRSId-stats.vcf.gz";

		String bedregions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		String outfile = "/Data/tmp/2016-06-29-quals/T1D-INFOAndAccMerged.txt";

		boolean includeId = true;
		boolean includeIndels = true;
		boolean includeICVariants = false;
		double mafthreshold = 0.01;


		try {
			MergeAccuracyAndInfoScoreFiles c = new MergeAccuracyAndInfoScoreFiles();
			c.ttest(files, labels, seqpanelvcf, bedregions, mafthreshold, variantsOnIC, includeId, includeIndels, includeICVariants, samplelist);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void ttest(
			String[] files,
			String[] labels,
			String seqpanelvcf,
			String bedregions,
			double mafthreshold,
			String variantsOnIC,
			boolean includeId,
			boolean includeIndels,
			boolean includeICVariants,
			String samplelist) throws IOException {


		HashSet<String> variantsOnICHash = loadVariantHash(variantsOnIC, includeId);
		Triple<ArrayList<VCFVariant>, ArrayList<VCFVariant>, ArrayList<VCFVariant>> seqpanelvariants = loadSequencedVariants(
				seqpanelvcf, bedregions, mafthreshold, 1d, variantsOnICHash, includeId, includeIndels, samplelist
		);

		ArrayList<VCFVariant> seqpanel = seqpanelvariants.getLeft();
		ArrayList<VCFVariant> variantsOnImmunoChip = seqpanelvariants.getMiddle();
		ArrayList<VCFVariant> variantsNotOnImmunoChip = seqpanelvariants.getRight();

		System.out.println(seqpanel.size() + " variants in VCF");
		System.out.println(variantsNotOnImmunoChip.size() + " not on IC");
		System.out.println(variantsOnImmunoChip.size() + " on IC");

		HashSet<String> sequencedVariantsHash = null;
		if (includeICVariants) {
			ArrayList<VCFVariant> allvars = new ArrayList<>();
			allvars.addAll(variantsNotOnImmunoChip);
			allvars.addAll(variantsOnImmunoChip);
			sequencedVariantsHash = hashit(allvars, includeId);
		} else {
			sequencedVariantsHash = hashit(variantsNotOnImmunoChip, includeId);
		}

		System.out.println(sequencedVariantsHash.size() + " after hashing");
		System.out.println(files.length + " files");

		// load qual scores
		double[][] values = new double[files.length][sequencedVariantsHash.size()];

		for (int f = 0; f < files.length; f++) {
			// get the imputation accuracies for these variants
			TextFile tf2 = new TextFile(files[f], TextFile.R);

			String[] elems = tf2.readLineElems(TextFile.tab);
			int nrSequenced = 0;
			int nrSequencedPassingRSQ = 0;
			int nrSequencedPassingMaf = 0;
			int nrSequencdPassingMafAndRSQ = 0;
			int ctr = 0;
			while (elems != null) {


				if (!elems[rsqlcol].equals("null")) {
					double val = Double.parseDouble(elems[rsqlcol]);
					double maf = Double.parseDouble(elems[maf2col]);

					String[] varElems = elems[0].split("_");

					boolean sequenced = isVariantInHash(varElems, sequencedVariantsHash, includeId);

					if (sequenced) {
						nrSequenced++;
						if (maf > mafthreshold) {
							nrSequencedPassingMaf++;

							//
							values[f][ctr] = val;
							ctr++;
						}

					}
				}

				elems = tf2.readLineElems(TextFile.tab);
			}
			tf2.close();

//			System.out.println("----");
//			System.out.println();

//			System.out.println("nrSequenced\t" + nrSequenced + "\t" + ((double) nrSequenced / sequencedVariantsHash.size()));
//			System.out.println("nrSequencedPassingRSQ\t" + nrSequencedPassingRSQ + "\t" + ((double) nrSequencedPassingRSQ / sequencedVariantsHash.size()) + "\t" + ((double) nrSequencedPassingRSQ / nrSequenced));
//			System.out.println("nrSequencedPassingMaf\t" + nrSequencedPassingMaf + "\t" + ((double) nrSequencedPassingMaf / sequencedVariantsHash.size()));
//			System.out.println("nrSequencdPassingMafAndRSQ\t" + nrSequencdPassingMafAndRSQ + "\t" + ((double) nrSequencdPassingMafAndRSQ / sequencedVariantsHash.size()));
//			System.out.println();

			System.out.println(nrSequenced + "\t" + nrSequencedPassingMaf + "\t" + nrSequencedPassingRSQ + "\t" + nrSequencdPassingMafAndRSQ);
		}

		// now make the table
		double[][] table = new double[files.length][files.length];
		double[] medians = new double[files.length];
		double[] means = new double[files.length];
		double[] stdevs = new double[files.length];
		for (int f = 0; f < files.length; f++) {
			for (int f2 = f + 1; f2 < files.length; f2++) {
				double[] x = values[f];
				double[] y = values[f2];

				org.apache.commons.math3.stat.inference.TTest t = new org.apache.commons.math3.stat.inference.TTest();
				double p = t.pairedTTest(x, y);
				double medianx = JSci.maths.ArrayMath.median(x);
				double mediany = JSci.maths.ArrayMath.median(y);

				table[f][f2] = p;
				table[f2][f] = (medianx - mediany);

			}
			medians[f] = JSci.maths.ArrayMath.median(values[f]);
			means[f] = JSci.maths.ArrayMath.mean(values[f]);
			stdevs[f] = JSci.maths.ArrayMath.standardDeviation(values[f]);
		}

		String header = "-";
		for (int i = 0; i < files.length; i++) {
			header += "\t" + labels[i];
		}

		System.out.println(header + "\tmedian\tmean\tstdev");

		for (int i = 0; i < files.length; i++) {
			String line = labels[i];
			for (int j = 0; j < files.length; j++) {
				if (i == j) {
					line += "\t-";
				} else {
					line += "\t" + table[i][j];
				}
			}
			System.out.println(line + "\t" + medians[i] + "\t" + means[i] + "\t" + stdevs[i]);
		}

	}

	public void merge(String seqpanelvcf, String variantsOnIC, String[] accfiles, String[] infofiles, String[] labels, String outfile) throws IOException {

		// get a list of maf > 0.005 variants that on the sequencingpanel
		TextFile tf = new TextFile(seqpanelvcf, TextFile.R);
		ArrayList<VCFVariant> seqpanel = new ArrayList<>();
		String ln = tf.readLine();

		double mafthreshold = 0.01;
		double upperthreshold = 1;
//		double infothreshold = 0.5;
		boolean includeICVariants = false;

		boolean includeId = true;
		boolean includeIndels = true;

		VariantCounter counter = new VariantCounter();
		HashSet<String> variantsOnICHash = counter.loadVariantHash(variantsOnIC, includeId);


		System.out.println(variantsOnICHash.size() + " total on IC");
		ArrayList<VCFVariant> variantsNotOnImmunoChip = new ArrayList<>();
		ArrayList<VCFVariant> variantsOnImmunoChip = new ArrayList<>();
		while (ln != null) {
			if (!ln.startsWith("#")) {
				String[] elems = ln.split("\t");
				Chromosome chr = Chromosome.parseChr(elems[0]);
				if (chr.isAutosome()) {
					VCFVariant variant = new VCFVariant(ln, VCFVariant.PARSE.ALL);
					if (variant.getMAF() > mafthreshold && variant.getMAF() < upperthreshold) {
						seqpanel.add(variant);
						boolean varOnIc = counter.isVariantInHash(elems, variantsOnICHash, includeId);
						boolean indel = counter.isIndel(elems);

						if (includeIndels || (!includeIndels && !indel)) {
							if (!varOnIc) {
								variantsNotOnImmunoChip.add(variant);
							} else {
								variantsOnImmunoChip.add(variant);
							}
						}
					}
				}
			}
			ln = tf.readLine();
		}
		tf.close();

		System.out.println(seqpanel.size() + " variants in VCF");
		System.out.println(variantsNotOnImmunoChip.size() + " not on IC");
		System.out.println(variantsOnImmunoChip.size() + " on IC");
		System.out.println("MAF> " + mafthreshold);
//		System.out.println("INFO> " + infothreshold);


		HashSet<String> sequencedVariantsHash = null;
		ArrayList<VCFVariant> allvars = null;
		if (includeICVariants) {
			allvars = new ArrayList<>();
			allvars.addAll(variantsNotOnImmunoChip);
			allvars.addAll(variantsOnImmunoChip);
			sequencedVariantsHash = counter.hashit(allvars, includeId);
		} else {
			allvars = new ArrayList<>();
			allvars.addAll(variantsNotOnImmunoChip);
			sequencedVariantsHash = counter.hashit(variantsNotOnImmunoChip, includeId);
		}

		System.out.println(allvars.size() + " variants to evaluate");

		HashMap<String, Integer> variantToId = new HashMap<String, Integer>();
		for (int v = 0; v < allvars.size(); v++) {
			ln = allvars.get(v).toString();
			variantToId.put(ln, v);
		}

		double[][][] scores = new double[infofiles.length][sequencedVariantsHash.size()][2]; // [referencepanels][variants][info/acc]

		for (int q = 0; q < accfiles.length; q++) {
			TextFile tf2 = new TextFile(accfiles[q], TextFile.R);
			String[] elems = tf2.readLineElems(TextFile.tab);
			while (elems != null) {
				String variantName = elems[0];
				Integer varId = variantToId.get(variantName);
				if (varId != null) {
					scores[q][varId][0] = Double.parseDouble(elems[12]);
				}
				elems = tf2.readLineElems(TextFile.tab);
			}
			tf2.close();
		}

		for (int q = 0; q < infofiles.length; q++) {
			TextFile tf2 = new TextFile(infofiles[q], TextFile.R);
			System.out.println("parsing: " + accfiles[q]);
			ln = tf2.readLine();
			while (ln != null) {
				if (!ln.startsWith("#")) {
					VCFVariant var = new VCFVariant(ln);
					String varname = var.toString();
					Integer varId = variantToId.get(varname);
					if (varId != null) {
						scores[q][varId][1] = var.getImputationQualityScore();
					}
				}
				ln = tf2.readLine();
			}
			tf.close();
		}

		TextFile out = new TextFile(outfile, TextFile.W);

		String header = "variant\tmafInSeqPanel";
		for (int l = 0; l < labels.length; l++) {
			header += "\t" + labels[l] + "-Accuracy\t" + labels[l] + "-INFO";
		}
		out.writeln(header);

		for (int v = 0; v < allvars.size(); v++) {
			ln = allvars.get(v).toString() + "\t" + allvars.get(v).getMAF();
			for (int l = 0; l < labels.length; l++) {
				ln += "\t" + scores[l][v][0] + "\t" + scores[l][v][1];
			}
			out.writeln(ln);
		}
		out.close();
	}

}
