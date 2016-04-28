package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.*;

/**
 * Created by hwestra on 12/8/15.
 */
public class VCFCorrelator {

	public void updateVCFInfoScore(String vcfin, String vcfOut) throws IOException {
		System.out.println("Will replace imputation quals.");
		System.out.println("In: " + vcfin);
		System.out.println("Out: " + vcfOut);
		TextFile tf = new TextFile(vcfin, TextFile.R);
		TextFile out = new TextFile(new File(vcfOut), TextFile.W, true);
		String ln = tf.readLine();
		int submitted = 0;

		int cores = Runtime.getRuntime().availableProcessors();
		System.out.println("Detected " + cores + " Processors ");
		ExecutorService threadPool = Executors.newFixedThreadPool(cores);
		CompletionService<String> jobHandler = new ExecutorCompletionService<>(threadPool);

		int nrRead = 0;

		while (ln != null) {
			if (ln.startsWith("#")) {
				out.writeln(ln);
			} else {
				VCFVariantInfoUpdater task = new VCFVariantInfoUpdater(ln);
				jobHandler.submit(task);
				submitted++;
				nrRead++;
				if (submitted % 10000 == 0) {
					int returned = 0;
					while (returned < submitted) {
						Future<String> future = null;
						try {
							future = jobHandler.take();
							if (future != null) {
								String outStr = future.get();
								out.writeln(outStr);
								returned++;
							}
						} catch (InterruptedException e) {
							e.printStackTrace();
						} catch (ExecutionException e) {
							e.printStackTrace();
						}
					}
					submitted = 0;
					System.out.println(nrRead + " variants read");
				}
			}
			ln = tf.readLine();
		}

		int returned = 0;
		while (returned < submitted) {
			Future<String> future = null;
			try {
				future = jobHandler.take();
				if (future != null) {
					String outStr = future.get();
					out.writeln(outStr);
					returned++;
				}
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}
		submitted = 0;

		System.out.println(nrRead + " variants total.");
		out.close();
		tf.close();
		threadPool.shutdown();
	}

	public void run(String vcf1, String vcf2, String variantsToTestFile, String out) throws IOException {

		// get samples in vcf1 and 2
		// index samples
		// load list of variants (if any)
		// load variants in one dataset
		// iterate the other dataset
		// correlate

		System.out.println("Output here: " + out);

		VCFGenotypeData data1 = new VCFGenotypeData(vcf1);
		VCFGenotypeData data2 = new VCFGenotypeData(vcf2);

		ArrayList<String> samples1 = data1.getSamples();
		System.out.println(samples1.size() + " samples in " + vcf1);
		ArrayList<String> samples2 = data2.getSamples();
		System.out.println(samples2.size() + " samples in " + vcf2);

		data2.close();

		HashSet<String> samples1Hash = new HashSet<String>();
		for (int i = 0; i < samples1.size(); i++) {
			samples1Hash.add(samples1.get(i));
		}

		ArrayList<String> sharedSamples = new ArrayList<String>();
		HashMap<String, Integer> sharedSamplesIndex = new HashMap<String, Integer>();

		int sctr = 0;
		for (String s : samples2) {
			if (samples1Hash.contains(s)) {
				sharedSamples.add(s);
				sharedSamplesIndex.put(s, sctr);
				sctr++;
			}
		}
		System.out.println(sharedSamples.size() + " samples shared between VCFs");

		HashSet<String> varsToTest = null;
		if (variantsToTestFile != null) {
			TextFile fin = new TextFile(variantsToTestFile, TextFile.R);
			ArrayList<String> str = fin.readAsArrayList();
			boolean splitweirdly = false;
			varsToTest = new HashSet<String>();
			for (String s : str) {
				String[] elems = s.split(";");
				if (elems.length > 1) {
					varsToTest.add(elems[0] + "-" + elems[1]);
					splitweirdly = true;
				} else {
					varsToTest.add(s);
				}

			}
			if (splitweirdly) {
				System.out.println("split weirdly ;)");
			}
			fin.close();
			System.out.println(varsToTest.size() + " variants to test from " + variantsToTestFile);
		}

		HashMap<String, VCFVariant> variantMap = new HashMap<String, VCFVariant>();
		int ctr1 = 0;
		while (data1.hasNext()) {
			VCFVariant var = data1.nextLoadHeader();
			Chromosome chr = Chromosome.parseChr(var.getChr());

			if (!chr.equals(Chromosome.X)) {
				String varStr = var.toString();
				if (varsToTest == null || varsToTest.contains(varStr)) {
					if (var.getTokens() != null) {
						var.parseGenotypes(var.getTokens(), VCFVariant.PARSE.ALL);
						var.cleartokens();
						variantMap.put(var.toString(), var);
					}
				}
			}
			ctr1++;
			if (ctr1 % 1000 == 0) {
				System.out.println(ctr1 + " variants parsed from vcf1");
			}
		}
		data1.close();

		System.out.println(variantMap.size() + " variants loaded from " + vcf1);
		TextFile tfot = new TextFile(out, TextFile.W);


		HashMap<String, VCFVariant> variantMap2 = new HashMap<>();
		HashSet<String> writtenVariants = new HashSet<String>();
		int ctr2 = 0;
		TextFile tfVCF2 = new TextFile(vcf2, TextFile.R);
		String ln = tfVCF2.readLine();

		while (ln != null) {
			if (!ln.startsWith("#")) {

				int strlen = ln.length();
				int substrlen = 500;
				if (strlen < substrlen) {
					substrlen = strlen;
				}
				String lnheader = ln.substring(0, substrlen);
				String[] lnheaderelems = lnheader.split("\t");

				// return this.chr + "_" + this.pos + "_" + this.id;
				String varStr = lnheaderelems[0] + "_" + lnheaderelems[1] + "_" + lnheaderelems[2];

				VCFVariant var1 = variantMap.get(varStr);
				if (var1 != null) {
					VCFVariant var2 = new VCFVariant(ln);

//					if (var2.getTokens() != null) {
//						var2.parseGenotypes(var2.getTokens(), VCFVariant.PARSE.ALL);
//						var2.cleartokens();
//					} else {
//						System.out.println(var2.toString());
//					}

					double[][] gprobs1 = null;
					HashMap<String, Double> info1 = var1.getInfo();
					HashMap<String, Double> info2 = var2.getInfo();
					if (var1.hasImputationProbs()) {
						gprobs1 = var1.getImputedDosages();
					} else {
						gprobs1 = var1.getDosages();
					}

					double[][] gprobs2 = null;
					if (var2.hasImputationProbs()) {
						gprobs2 = var2.getImputedDosages(); // format [samples][alleles]
					} else {
						gprobs2 = var2.getDosages();
					}


					// check if variants have the same number of alleles
					if (gprobs1[0].length == gprobs2[0].length) {

						// recode: make sure sample ordering is identical

						gprobs1 = reorder(gprobs1, samples1, sharedSamplesIndex, sharedSamples.size());
						gprobs2 = reorder(gprobs2, samples2, sharedSamplesIndex, sharedSamples.size());

						// remove missing genotypes
						Pair<double[][], double[][]> data = removeNulls(gprobs1, gprobs2);

						for (int a = 0; a < gprobs1[0].length; a++) {
							// do some remapping here
							double[] x = toArr(data.getLeft(), a);
							double[] y = toArr(data.getRight(), a);

							double r = JSci.maths.ArrayMath.correlation(x, y);
							var1.recalculateMAFAndCallRate();
							var2.recalculateMAFAndCallRate();

							String var1Str = var1.getMinorAllele() + "\t" + Strings.concat(var1.getAlleles(), Strings.comma) + "\t" + var1.getMAF() + "\t" + var1.getCallrate();
							String var2Str = var2.getMinorAllele() + "\t" + Strings.concat(var2.getAlleles(), Strings.comma) + "\t" + var2.getMAF() + "\t" + var2.getCallrate();

							if (Double.isNaN(r)) {
								ln = var1.toString() + "\t" + var1Str + "\t" + var2Str + "\t" + (a + 1) + "\t" + data.getLeft().length + "\t" + 0 + "\t" + 0 + "\t" + var1.getImputationQualityScore() + "\t" + var2.getImputationQualityScore();
							} else {
								double rsq = r * r;
								ln = var1.toString() + "\t" + var1Str + "\t" + var2Str + "\t" + (a + 1) + "\t" + data.getLeft().length + "\t" + r + "\t" + rsq + "\t" + var1.getImputationQualityScore() + "\t" + var2.getImputationQualityScore();
							}
							tfot.writeln(ln);
							writtenVariants.add(varStr);
						}
					} else {
						// ?
					}
				} else {
					if (varsToTest != null && varsToTest.contains(varStr)) {
						VCFVariant var2 = new VCFVariant(ln);
						variantMap2.put(varStr, var2);
					}
				}
				ctr2++;
				if (ctr2 % 1000 == 0) {
					System.out.println(ctr2 + " variants parsed from vcf2");
				}
			}
			ln = tfVCF2.readLine();
		}
		tfVCF2.close();

		// write variants that are not in variantlist
		if (varsToTest != null) {
			for (String variant : varsToTest) {

				if (!writtenVariants.contains(variant)) {
					VCFVariant var1 = variantMap.get(variant);
					VCFVariant var2 = variantMap2.get(variant);

					if ((var1 == null || var2 == null) || (var1 == null && var2 == null)) {

						String var1Str = "";
						String var2Str = "";
						String infoStr1 = null;
						String infoStr2 = null;
						String varstr = variant;


						if (var1 == null && var2 == null) {
							var1Str = null + "\t" + null + "\t" + 0 + "\t" + 0;
							var2Str = null + "\t" + null + "\t" + 0 + "\t" + 0;
						} else if (var1 != null && var2 == null) {
							varstr = var1.toString();
							infoStr1 = "" + var1.getInfo().get("AR2");
							var1Str = var1.getMinorAllele() + "\t" + Strings.concat(var1.getAlleles(), Strings.comma) + "\t" + var1.getMAF() + "\t" + var1.getCallrate();
							var2Str = null + "\t" + null + "\t" + 0;
						} else if (var1 == null && var2 != null) {
							varstr = var2.toString();
							infoStr2 = "" + var2.getInfo().get("AR2");
							var1Str = null + "\t" + null + "\t" + 0;
							var2Str = var2.getMinorAllele() + "\t" + Strings.concat(var2.getAlleles(), Strings.comma) + "\t" + var2.getMAF() + "\t" + var2.getCallrate();
						}
						ln = varstr + "\t" + var1Str + "\t" + var2Str + "\t" + 0 + "\t" + 0 + "\t" + 0 + "\t" + 0 + "\t" + infoStr1 + "\t" + infoStr2;
						tfot.writeln(ln);
					}
				}

			}
		}


		tfot.close();
	}

	private Pair<double[][], double[][]> removeNulls(double[][] gprobs1, double[][] gprobs2) {

		int nrNull = 0;
		for (int i = 0; i < gprobs1.length; i++) {
			double d = gprobs1[i][0];
			double d2 = gprobs2[i][0];
			if (Double.isNaN(d) || Double.isNaN(d2)) {
				nrNull++;
			}
		}

		double[][] out1 = new double[gprobs1.length - nrNull][];
		double[][] out2 = new double[gprobs1.length - nrNull][];
		int ctr = 0;
		for (int i = 0; i < gprobs1.length; i++) {
			double d = gprobs1[i][0];
			double d2 = gprobs2[i][0];
			if (!Double.isNaN(d) && !Double.isNaN(d2)) {
				out1[ctr] = gprobs1[i];
				out2[ctr] = gprobs2[i];
				ctr++;
			}
		}

		return new Pair<double[][], double[][]>(out1, out2);
	}

	private double[][] reorder(double[][] gprobs2, ArrayList<String> samples, HashMap<String, Integer> sharedSamplesIndex, int size) {
		double[][] output = new double[size][gprobs2[0].length];

		// initialize with NaNs
		for (int i = 0; i < output.length; i++) {
			for (int j = 0; j < output[i].length; j++) {
				output[i][j] = Double.NaN;
			}
		}

		for (int i = 0; i < gprobs2.length; i++) {
			String sample = samples.get(i);
			Integer index = sharedSamplesIndex.get(sample);
			if (index != null) {
				for (int j = 0; j < gprobs2[i].length; j++) {
					output[index][j] = gprobs2[i][j];
				}
			}
		}
		return output;
	}

	private double[][] recode1(double[][] gprobs1, boolean[] includeSample1, int size) {
		double[][] output = new double[size][gprobs1[0].length];
		int ctr = 0;
		for (int i = 0; i < gprobs1.length; i++) {
			if (includeSample1[i]) {
				for (int j = 0; j < gprobs1[i].length; j++) {
					output[ctr][j] = gprobs1[i][j];
				}
				ctr++;
			}
		}
		return output;
	}

	private double[] toArr(double[][] x, int q) {
		double[] arr = new double[x.length];
		for (int i = 0; i < arr.length; i++) {
			arr[i] = x[i][q];
		}
		return arr;
	}

	public class VCFVariantInfoUpdater implements Callable<String> {

		private String in = null;

		public VCFVariantInfoUpdater(String in) {
			this.in = in;
		}

		@Override
		public String call() throws Exception {

			VCFVariant variant = new VCFVariant(in);
			int nrAlleles = variant.getAlleles().length;
			double rsq = 0;
			if (nrAlleles == 2) {
				double[] probabilies = variant.getGenotypeDosages();
				double[] bestguess = convertGenotypesToDouble(variant);
				double r = JSci.maths.ArrayMath.correlation(bestguess, probabilies);
				rsq = r * r;
			}

			String[] elems = in.split("\t");
			elems[7] = "INFO=" + rsq;

			return Strings.concat(elems, Strings.tab);
		}

		private double[] convertToProbsToDouble(VCFVariant vcfVariant) {
			double[][] dosages = vcfVariant.getImputedDosages(); // [samples][alleles];
			double[] output = new double[dosages.length];
			for (int i = 0; i < output.length; i++) {
				output[i] = dosages[i][0];
			}
			return output;
		}

		private double[] convertGenotypesToDouble(VCFVariant vcfVariant) {
			byte[][] alleles = vcfVariant.getGenotypeAlleles();
			double[] output = new double[alleles[0].length];
			for (int i = 0; i < alleles[0].length; i++) {
				if (alleles[0][i] == -1) {
					output[i] = -1;
				} else {
					output[i] = (alleles[0][i] + alleles[1][i]);
				}
			}

			return output;
		}


	}

}
