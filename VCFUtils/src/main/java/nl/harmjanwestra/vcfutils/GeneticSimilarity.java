package nl.harmjanwestra.vcfutils;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Created by Harm-Jan on 07/02/16.
 */
public class GeneticSimilarity {


	public static void main(String[] args) {
//		GeneticSimilarity sim = new GeneticSimilarity();
//		String listofvariants = "/Data/tmp/2016-07-22/plink.prune.out";
//		String vcf1 = "/Data/ImmunoChip/RA/2015-09-23-IncX/ES/genotypes-filtered-sorted.vcf.gz";
//		String vcf2 = "/Data/ImmunoChip/RA/2015-09-23-IncX/NL/genotypes-filtered-sorted.vcf.gz";
//		String outfile = "/Data/tmp/2016-07-22/debug";
//		try {
//			sim.determineGeneticSimilaritySingleDataset(listofvariants, vcf1, outfile);
////			sim.determineGeneticSimilarityBetweenDatasets(listofvariants, vcf1, vcf2, outfile);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
	}

	public void determineGeneticSimilaritySingleDataset(String listOfVariants, boolean exclude, String vcf1, String outfile) throws IOException {
		// read list of variants
		TextFile tf1 = new TextFile(listOfVariants, TextFile.R);
		ArrayList<String> listOfVariantsToExcludeArr = tf1.readAsArrayList();
		tf1.close();
		HashSet<String> hashSetVariants = new HashSet<String>();
		hashSetVariants.addAll(listOfVariantsToExcludeArr);

		// load the variants for ds1
		System.out.println("VCF input: " + vcf1);
		Pair<ArrayList<String>, ArrayList<VCFVariant>> data1 = loadData(vcf1, hashSetVariants, exclude);
		ArrayList<String> samples1 = data1.getLeft();
		ArrayList<VCFVariant> variants1 = data1.getRight();
		System.out.println(samples1.size() + " samples in vcf1");
		System.out.println(variants1.size() + " variants in vcf1");

		VCFVariant[] allvariants = variants1.toArray(new VCFVariant[0]);
		nl.harmjanwestra.utilities.vcf.GeneticSimilarity sim = new nl.harmjanwestra.utilities.vcf.GeneticSimilarity();
		int cores = Runtime.getRuntime().availableProcessors();
		System.out.println("Detected " + cores + " Processors ");
		ExecutorService threadPool = Executors.newFixedThreadPool(cores);

		sim.calculateWithinDataset(allvariants, threadPool);
		Triple<DoubleMatrix2D, DoubleMatrix2D, DoubleMatrix2D> similarity = sim.calculateWithinDataset(allvariants, threadPool);
		DoubleMatrix2D geneticSimilarity = similarity.getLeft();
		DoubleMatrix2D sharedgenotypes = similarity.getMiddle();
		DoubleMatrix2D correlationmatrix = similarity.getRight();

		writeMatrix(geneticSimilarity, samples1.toArray(new String[0]), samples1.toArray(new String[0]), outfile + "geneticsim.txt.gz");
		writeMatrix(correlationmatrix, samples1.toArray(new String[0]), samples1.toArray(new String[0]), outfile + "correlations.txt.gz");
		writeMatrix(sharedgenotypes, samples1.toArray(new String[0]), samples1.toArray(new String[0]), outfile + "sharedgenotypes.txt.gz");

		// make two dosage matrices
		System.out.println("Submitting correlation jobs..");
		// correlate the samples

		String header = "Sample1\tSample2\tCorrelation\tGeneticSimilarity\tPercSharedGenotypes";
		TextFile out = new TextFile(outfile + "pairsAboveCorrelationThreshold0.25.txt", TextFile.W);
		out.writeln(header);
		for (int i = 0; i < samples1.size(); i++) {
			for (int j = i; j < samples1.size(); j++) {
				if (correlationmatrix.get(i, j) > 0.25) {
					out.writeln(samples1.get(i) + "\t" + samples1.get(j) + "\t" + correlationmatrix.get(i, j) + "\t" + geneticSimilarity.get(i, j) + "\t" + sharedgenotypes.get(i, j));
				}
			}
		}
		out.close();


		TextFile out3 = new TextFile(outfile + "toppairpersample1.txt", TextFile.W);
		out3.writeln(header);
		for (int i = 0; i < samples1.size(); i++) {
			double maxx = 0;
			int maxj = 0;
			for (int j = 0; j < samples1.size(); j++) {
				if (correlationmatrix.get(i, j) > maxx) {
					maxx = correlationmatrix.get(i, j);
					maxj = j;
				}
			}
			out3.writeln(samples1.get(i) + "\t" + samples1.get(maxj) + "\t" + correlationmatrix.get(i, maxj) + "\t" + geneticSimilarity.get(i, maxj) + "\t" + sharedgenotypes.get(i, maxj));
		}
		out3.close();
	}

	public void determineGeneticSimilarityBetweenDatasets(String listOfVariantsFile, boolean exclude, String vcf1, String vcf2, String outfile) throws IOException {

		// read list of variants
		TextFile tf1 = new TextFile(listOfVariantsFile, TextFile.R);
		ArrayList<String> listOfVariantsToExcludeArr = tf1.readAsArrayList();
		tf1.close();
		HashSet<String> hashSetVariantsToExclude = new HashSet<String>();
		hashSetVariantsToExclude.addAll(listOfVariantsToExcludeArr);

		// load the variants for ds1
		System.out.println("VCF1: " + vcf1);
		System.out.println("VCF2: " + vcf2);
		Pair<ArrayList<String>, ArrayList<VCFVariant>> data1 = loadData(vcf1, hashSetVariantsToExclude, exclude);
		ArrayList<String> samples1 = data1.getLeft();
		ArrayList<VCFVariant> variants1 = data1.getRight();
		System.out.println(samples1.size() + " samples in vcf1");
		System.out.println(variants1.size() + " variants in vcf1");

		Pair<ArrayList<String>, ArrayList<VCFVariant>> data2 = loadData(vcf2, hashSetVariantsToExclude, exclude);
		ArrayList<String> samples2 = data2.getLeft();
		ArrayList<VCFVariant> variants2 = data2.getRight();
		System.out.println(samples2.size() + " samples in vcf2");
		System.out.println(variants2.size() + " variants in vcf2");

		// merge set of loaded variants
		HashSet<String> hashOfVariants = new HashSet<String>();
		for (VCFVariant v : variants1) {
			hashOfVariants.add(v.toString());
		}
		for (VCFVariant v : variants2) {
			hashOfVariants.add(v.toString());
		}

		ArrayList<String> listOfVariants = new ArrayList<>();
		listOfVariants.addAll(hashOfVariants);
		HashMap<String, Integer> variantToId = new HashMap<String, Integer>();
		int varctr = 0;
		for (String s : listOfVariants) {
			variantToId.put(s, varctr);
			varctr++;
		}
		System.out.println(variantToId.size() + " unique variants..");
		VCFVariant[][] allvariants = new VCFVariant[2][listOfVariants.size()];
		for (VCFVariant v : variants1) {
			Integer id = variantToId.get(v.toString());
			if (id != null) {
				allvariants[0][id] = v;
			}
		}
		for (VCFVariant v : variants2) {
			Integer id = variantToId.get(v.toString());
			if (id != null) {
				allvariants[1][id] = v;
			}
		}

		// prune the set of variants
		int missing = 0;
		for (int i = 0; i < listOfVariants.size(); i++) {
			if (allvariants[0][i] == null || allvariants[1][i] == null) {
				missing++;
			}
		}
		System.out.println(missing + " variants not present in VCF1 or VCF2.");
		VCFVariant[][] tmpvariants = new VCFVariant[2][listOfVariants.size() - missing];
		int present = 0;
		for (int i = 0; i < listOfVariants.size(); i++) {
			if (allvariants[0][i] != null && allvariants[1][i] != null) {
				tmpvariants[0][present] = allvariants[0][i];
				tmpvariants[1][present] = allvariants[1][i];
				present++;
			}
		}
		allvariants = tmpvariants;

		System.out.println(allvariants[0].length + " variants shared between datasets");
		if (allvariants[0].length == 0) {
			System.exit(-1);
		}

		// now see if the alleles match
		VCFMerger merger = new VCFMerger();
		for (int i = 0; i < allvariants[0].length; i++) {
			VCFVariant variant1 = allvariants[0][i];
			VCFVariant variant2 = allvariants[1][i];

			String[] alleles1 = variant1.getAlleles();
			String allele1minor = variant1.getMinorAllele();
			String[] alleles2 = variant2.getAlleles();
			String allele2minor = variant2.getMinorAllele();

			int nrIdenticalAlleles = merger.countIdenticalAlleles(alleles1, alleles2);

			if (nrIdenticalAlleles != 2 || !allele1minor.equals(allele2minor)) {
				allvariants[0][i] = null;
				allvariants[1][i] = null;
			}
		}

		// prune the set of variants
		missing = 0;
		for (int i = 0; i < allvariants[0].length; i++) {
			if (allvariants[0][i] == null || allvariants[1][i] == null) {
				missing++;
			}
		}
		tmpvariants = new VCFVariant[2][allvariants[0].length - missing];


		present = 0;
		for (int i = 0; i < allvariants[0].length; i++) {
			if (allvariants[0][i] != null && allvariants[1][i] != null) {
				tmpvariants[0][present] = allvariants[0][i];
				tmpvariants[1][present] = allvariants[1][i];
				present++;
			}
		}
		allvariants = tmpvariants;
		System.out.println(allvariants[0].length + " variants have identical alleles");
		if (allvariants.length == 0) {
			System.exit(-1);
		}

		if (allvariants[0].length < 100) {
			System.out.println("Error: not enough info...");
			System.out.println(allvariants[0].length + " variants remain after pruning");
			System.exit(-1);
		}


		// get the genetic similarity
		nl.harmjanwestra.utilities.vcf.GeneticSimilarity sim = new nl.harmjanwestra.utilities.vcf.GeneticSimilarity();

		int cores = Runtime.getRuntime().availableProcessors();
		System.out.println("Detected " + cores + " Processors ");

		ExecutorService threadPool = Executors.newFixedThreadPool(cores);
		Triple<DoubleMatrix2D, DoubleMatrix2D, DoubleMatrix2D> similarity = sim.calculateBetweenDatasets(allvariants, threadPool);
		threadPool.shutdown();

		DoubleMatrix2D geneticSimilarity = similarity.getLeft();
		DoubleMatrix2D sharedgenotypes = similarity.getMiddle();
		DoubleMatrix2D correlationmatrix = similarity.getRight();

		writeMatrix(geneticSimilarity, samples1.toArray(new String[0]), samples2.toArray(new String[0]), outfile + "geneticsim.txt.gz");
		writeMatrix(correlationmatrix, samples1.toArray(new String[0]), samples2.toArray(new String[0]), outfile + "correlations.txt.gz");
		writeMatrix(sharedgenotypes, samples1.toArray(new String[0]), samples2.toArray(new String[0]), outfile + "sharedgenotypes.txt.gz");


		// make two dosage matrices
		System.out.println("Submitting correlation jobs..");
		// correlate the samples

		String header = "Sample1\tSample2\tCorrelation\tGeneticSimilarity\tPercSharedGenotypes";
		TextFile out = new TextFile(outfile + "pairsAboveCorrelationThreshold0.25.txt", TextFile.W);
		out.writeln(header);
		for (int i = 0; i < samples1.size(); i++) {
			for (int j = 0; j < samples2.size(); j++) {
				if (correlationmatrix.get(i, j) > 0.25) {
					out.writeln(samples1.get(i) + "\t" + samples2.get(j) + "\t" + correlationmatrix.get(i, j) + "\t" + geneticSimilarity.get(i, j) + "\t" + sharedgenotypes.get(i, j));
				}
			}
		}
		out.close();


		TextFile out3 = new TextFile(outfile + "toppairpersample1.txt", TextFile.W);
		out3.writeln(header);
		for (int i = 0; i < samples1.size(); i++) {
			double maxx = 0;
			int maxj = 0;
			for (int j = 0; j < samples2.size(); j++) {
				if (correlationmatrix.get(i, j) > maxx) {
					maxx = correlationmatrix.get(i, j);
					maxj = j;
				}
			}
			out3.writeln(samples1.get(i) + "\t" + samples2.get(maxj) + "\t" + correlationmatrix.get(i, maxj) + "\t" + geneticSimilarity.get(i, maxj) + "\t" + sharedgenotypes.get(i, maxj));
		}
		out3.close();

		TextFile out4 = new TextFile(outfile + "toppairpersample2.txt", TextFile.W);
		out4.writeln(header);
		for (int j = 0; j < samples2.size(); j++) {
			double maxx = 0;
			int maxi = 0;
			for (int i = 0; i < samples1.size(); i++) {

				if (correlationmatrix.get(i, j) > maxx) {
					maxx = correlationmatrix.get(i, j);
					maxi = i;
				}
			}
			out4.writeln(samples1.get(maxi) + "\t" + samples2.get(j) + "\t" + correlationmatrix.get(maxi, j) + "\t" + geneticSimilarity.get(maxi, j) + "\t" + sharedgenotypes.get(maxi, j));
		}
		out4.close();

		threadPool.shutdown();
	}

	private void writeMatrix(DoubleMatrix2D geneticSimilarity, String[] samplesRows, String[] samplesCols, String s) throws IOException {
		TextFile out = new TextFile(s, TextFile.W);
		String header = "-";
		header += "\t" + Strings.concat(samplesCols, Strings.tab);
		out.writeln(header);
		DecimalFormat f = new DecimalFormat("#.###");
		ProgressBar pb = new ProgressBar(geneticSimilarity.rows(), "printing matrix");
		for (int r = 0; r < geneticSimilarity.rows(); r++) {
			DoubleMatrix1D row = geneticSimilarity.viewRow(r);
			double[] rowarr = row.toArray();
			out.writeln(samplesRows[r] + "\t" + Strings.concat(rowarr, Strings.tab, f));

			pb.iterate();
		}
		out.close();
		pb.close();
	}

	private Pair<ArrayList<String>, ArrayList<VCFVariant>> loadData(String vcffileloc, HashSet<String> hashsetVariants, boolean exclude) throws IOException {

		String[] files = vcffileloc.split(",");

		ArrayList<String> samples = null;
		ArrayList<VCFVariant> variants = new ArrayList<>();
		for (String file : files) {
			System.out.println("Parsing: " + file);
			VCFGenotypeData data1 = new VCFGenotypeData(file);
			if (samples == null) {
				samples = data1.getSamples();
			}
			data1.close();

			TextFile tf = new TextFile(file, TextFile.R);
			String ln = tf.readLine();

			int lnsparsed = 0;
			while (ln != null) {
				if (!ln.startsWith("#")) {

					String substr = ln;
					if (ln.length() > 150) {
						ln.substring(0, 150);
					}
					String[] elems = substr.split("\t");
					String varId = elems[2];
					if ((exclude && !hashsetVariants.contains(varId))
							||
							(!exclude && hashsetVariants.contains(varId))
							) {
						VCFVariant variant = new VCFVariant(ln);
						if (variant.isAutosomal() && variant.isBiallelic() && variant.getMAF() > 0.05 && variant.getHwep() > 0.001) {
							variants.add(variant);
						}
					}
				}
				lnsparsed++;
				if (lnsparsed % 1000 == 0) {
					System.out.print(lnsparsed + " lines parsed. " + variants.size() + " in memory\r");
				}
				ln = tf.readLine();
			}
			tf.close();
			System.out.println();
			System.out.println("Done parsing: " + file);
			System.out.println("");
		}
		return new Pair<>(samples, variants);
	}

	public class CorrelationTask implements Callable<Triple<Integer, Integer, Double>> {

		int i;
		int j;
		double[] dosage1;
		double[] dosage2;

		public CorrelationTask(int i, int j, double[] dosage1, double[] dosage2) {
			this.i = i;
			this.j = j;
			this.dosage1 = dosage1;
			this.dosage2 = dosage2;
		}

		@Override
		public Triple<Integer, Integer, Double> call() throws Exception {
			// prune missing
			int missinggt = 0;
			for (int k = 0; k < dosage1.length; k++) {
				if (dosage1[k] == -1 || dosage2[k] == -1) {
					missinggt++;
				}
			}

			double[] x = new double[dosage1.length - missinggt];
			double[] y = new double[dosage2.length - missinggt];

			int presentgt = 0;
			for (int k = 0; k < dosage1.length; k++) {
				if (dosage1[k] == -1 || dosage2[k] == -1) {

				} else {
					x[presentgt] = dosage1[k];
					y[presentgt] = dosage2[k];
					presentgt++;
				}
			}

			double corr = JSci.maths.ArrayMath.correlation(x, y);
			return new Triple<Integer, Integer, Double>(i, j, corr);
		}
	}

}
