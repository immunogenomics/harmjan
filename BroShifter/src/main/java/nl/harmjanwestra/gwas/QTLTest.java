package nl.harmjanwestra.gwas;

import cern.jet.random.tdouble.StudentT;
import nl.harmjanwestra.gwas.CLI.QTLTestOptions;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.enums.Strand;
import nl.harmjanwestra.utilities.features.*;
import nl.harmjanwestra.utilities.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;

/**
 * Created by hwestra on 3/14/16.
 */
public class QTLTest {

	private final QTLTestOptions options;

	public QTLTest(QTLTestOptions options) throws IOException {
		this.options = options;
		this.run();
	}


	// TODO: this code does not assume missing genotypes!!!

	public void run() throws IOException {

		String outfile = options.out;
		String gtfAnnotationFile = options.annotation;
		String expDataFile = options.expression;
		String covDataFile = options.covariates;
		String genotypeDataFile = options.genotype;
		String genotypeToExpressionCouplingFile = options.genotypeToExpression;
		String geneLimitFile = options.genelimit;
		String variantLimitFile = options.variantlimit;
		String covariateLimitFile = options.covariatelimit;
		double imputationqualthreshold = options.impqualthreshold;

		int cisWindow = options.ciswindow;
		Random random = new Random();
		long[] permutationSeedNumber = new long[options.nrpermutationspergene];
		for (int i = 0; i < permutationSeedNumber.length; i++) {
			permutationSeedNumber[i] = random.nextLong(); // To think about: do we want a fixed order in our permutations? if so, each permutation should get a unique seed, shared between each gene
		}


		double mafThreshold = options.mafthreshold;
		double callrateThreshold = options.callratethreshold;
		int nrPermutationsPerVariant = options.nrpermutationspergene;

		HashSet<String> limitGenes = null;
		HashSet<String> limitVariants = null;
		HashSet<String> limitCovarariates = null;


		// load links between genotype and gene expression if any
		HashMap<String, String> genotypeToExpression = null;
		if (genotypeToExpressionCouplingFile != null) {
			System.out.println("Loading GTE: " + genotypeToExpressionCouplingFile);
			genotypeToExpression = new HashMap<String, String>();
			TextFile tf = new TextFile(genotypeToExpressionCouplingFile, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				if (elems.length >= 2) {
					genotypeToExpression.put(elems[0], elems[1]);
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			System.out.println(genotypeToExpression.size() + " sample couplings loaded");
		}

		// initialize genotype reader

		System.out.println("Initializing VCF reader: " + genotypeDataFile);
		VCFGenotypeData genotypeVCF = new VCFGenotypeData(genotypeDataFile);
		ArrayList<String> genotypeSamples = genotypeVCF.getSamples();
		System.out.println(genotypeSamples.size());
		HashSet<String> expressionSamplesToLoad = new HashSet<String>();
		for (int i = 0; i < genotypeSamples.size(); i++) {
			if (genotypeToExpression == null) {
				expressionSamplesToLoad.add(genotypeSamples.get(i));
			} else {
				String sample = genotypeToExpression.get(genotypeSamples.get(i));
				if (sample != null) {
					expressionSamplesToLoad.add(sample);
				}
			}
		}
		System.out.println(genotypeSamples.size() + " genotype samples, " + expressionSamplesToLoad.size() + " defined by limits");

		// load covariates .. limit to samples present in genotype data
		DoubleMatrixDataset<String, String> cov = null;
		if (covDataFile != null) {
			System.out.println("Loading covariate data: " + covDataFile);
			cov = DoubleMatrixDataset.loadSubsetOfTextDoubleData(covDataFile, "\t", limitCovarariates, expressionSamplesToLoad);
			ArrayList<String> covariateSamples = cov.getColObjects();
			expressionSamplesToLoad = new HashSet<String>();
			expressionSamplesToLoad.addAll(covariateSamples);
			System.out.println(expressionSamplesToLoad.size() + " samples with covariate and genotype data.");
		}

		// load gene expression data .. limit to samples present in both covariate and gene expression data
		System.out.println("Loading expression data: " + expDataFile);
		DoubleMatrixDataset<String, String> exp = DoubleMatrixDataset.loadSubsetOfTextDoubleData(expDataFile, "\t", limitGenes, expressionSamplesToLoad);
		ArrayList<String> expressionSamples = exp.getColObjects();

		if (expressionSamples.isEmpty()) {
			System.err.println("Error: no samples remain after loading expression data. Probably a sample linking issue.");
			System.exit(-1);
		} else {
			System.out.println(expressionSamplesToLoad.size() + " samples with covariate, genotype data, and expression data.");
		}
		expressionSamplesToLoad = new HashSet<String>();
		expressionSamplesToLoad.addAll(expressionSamples);

		// now we need to prune the covariate data with the gene expression samples. Quickest way is to just reload.
		if (covDataFile != null) {
			cov = DoubleMatrixDataset.loadSubsetOfTextDoubleData(covDataFile, "\t", limitCovarariates, expressionSamplesToLoad);
		}

		// put expression data in same order as expression data:
		// -- index genotyped samples that also have covariate/genotype data
		boolean[] includeGenotypeSample = new boolean[genotypeSamples.size()];
		LinkedHashMap<String, Integer> newSampleOrder = new LinkedHashMap<String, Integer>();
		int ctr = 0;
		for (int i = 0; i < genotypeSamples.size(); i++) {
			String genotypeSample = genotypeSamples.get(i);
			String expressionSample = genotypeSample;
			if (genotypeToExpression != null) {
				expressionSample = genotypeToExpression.get(genotypeSample);
			}
			if (expressionSamplesToLoad.contains(expressionSample)) {
				includeGenotypeSample[i] = true;
				newSampleOrder.put(expressionSample, ctr);
				ctr++;
			}
		}
		System.out.println(newSampleOrder.size() + " samples shared between genotypes and expression");

		// reorder covariate and expression data
		if (cov != null) {
			cov.reorderCols(newSampleOrder);
		}
		exp.reorderCols(newSampleOrder);

		// load the gene annotation
		Collection<Gene> allGenes = null;
		if (gtfAnnotationFile.endsWith("gtf")) {
			GTFAnnotation geneAnnotation = new GTFAnnotation(gtfAnnotationFile);
			allGenes = geneAnnotation.getGenes();
		} else {
			TextFile tf = new TextFile(gtfAnnotationFile, TextFile.R);

			String[] elems = tf.readLineElems(TextFile.tab);
			allGenes = new ArrayList<>();
			while (elems != null) {
				String name = elems[1];
				Chromosome chr = Chromosome.parseChr(elems[3]);
				if (chr.equals(Chromosome.NA)) {

				} else {
					Integer start = Integer.parseInt(elems[4]);
					Integer stop = Integer.parseInt(elems[5]);
					Gene g = new Gene(elems[2], chr, Strand.POS, start, stop);
					g.setGeneId(name);
					allGenes.add(g);
//					System.out.println(g.toString());
				}

				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			System.out.println(allGenes.size() + " genes loaded from: " + gtfAnnotationFile);
		}


		ArrayList<Gene> finalGeneSetArr = new ArrayList<>();
		LinkedHashMap<String, Integer> expGeneHash = exp.getHashRows();
		LinkedHashMap<Gene, Integer> geneToExpGene = new LinkedHashMap<Gene, Integer>();

		HashSet<String> geneNamesFound = new HashSet<>();
		for (Gene gene : allGenes) {
			String name = gene.getGeneId();
			if (expGeneHash.containsKey(name)) {
				Integer id = expGeneHash.get(name);

				finalGeneSetArr.add(gene);
				geneToExpGene.put(gene, id);
				geneNamesFound.add(name);
			}
		}

		FeatureComparator comp = new FeatureComparator(false);
		TreeSet<Gene> finalGeneSet = new TreeSet<>(comp);
		finalGeneSet.addAll(finalGeneSetArr);
		comp.setAllowOverlap(true);

		TextFile notfound = new TextFile(outfile + "-expressionNotFound.txt", TextFile.W);
		ArrayList<String> expgenes = exp.getRowObjects();
		for (String gene : expgenes) {
			if (!geneNamesFound.contains(gene)) {
				notfound.writeln(gene);
			}
		}
		notfound.close();

		if (geneToExpGene.isEmpty()) {
			System.err.println("Error: no genes found. Probably a mismatch between annotation and gene expression path");
			System.exit(-1);
		}
		System.out.println(geneToExpGene.size() + " genes with annotation found in gene expression data");
		// rank the gene expression data -- save for later

		// for each chromosome
		// load all genotypes (load only biallelic markers)
		// put the variants in a treemap
		FeatureComparator comp2 = new FeatureComparator(false);
		TreeSet<Feature> variantSet = new TreeSet<Feature>(comp2);
		LinkedHashMap<Feature, VCFVariant> allGenotypes = new LinkedHashMap<>();
		System.out.println("Loading variants");
		HashSet<Gene> genesOverlappingVariants = new HashSet<Gene>();
		int variantctr = 0;
		genotypeVCF.close();
		TextFile tf1 = new TextFile(genotypeDataFile, TextFile.R);
		String vcfln = tf1.readLine();
		while (vcfln != null) {
			if (!vcfln.startsWith("#")) {
				VCFVariant variant = new VCFVariant(vcfln, VCFVariant.PARSE.HEADER);
				if(variant.getChrObj().isAutosome()){
					Double impqual = variant.getImputationQualityScore();
					if (impqual == null || impqual > imputationqualthreshold) {
						variant = new VCFVariant(vcfln, VCFVariant.PARSE.ALL);
						if (variant.getAlleles().length == 2
								&& variant.getMAF() > mafThreshold
								&& variant.getCallrate() > callrateThreshold) {
							int minPosition = variant.getPos() - cisWindow;
							if (minPosition < 0) {
								minPosition = 0;
							}
							int maxPosition = variant.getPos() + cisWindow;
							if (maxPosition > Integer.MAX_VALUE) {
								maxPosition = Integer.MAX_VALUE;
							}

							// load only variants near genes
							Chromosome chr = Chromosome.parseChr(variant.getChr());
							Gene geneStart = new Gene("", chr, Strand.POS, minPosition, minPosition);
							Gene geneStop = new Gene("", chr, Strand.POS, maxPosition, maxPosition);
							SortedSet<Gene> overlappingGenes = finalGeneSet.subSet(geneStart, true, geneStop, true);

							if (!overlappingGenes.isEmpty()) {
								genesOverlappingVariants.addAll(overlappingGenes);
								variantSet.add(variant.asFeature());
								allGenotypes.put(variant.asFeature(), variant);
							}
						}
					}

					variantctr++;
					if (variantctr % 1000 == 0) {
						System.out.print("\r" + variantctr + " variants parsed. " + allGenotypes.size() + " in memory");
					}
				}

			}

			vcfln = tf1.readLine();
		}
		tf1.close();
		comp2.setAllowOverlap(true);
		System.out.println();
		System.out.println("Done loading genotypes");
		genotypeVCF.close();

		if (allGenotypes.isEmpty()) {
			System.err.println("No variants overlap with available genes. Probably a problem with the variant annotation?");
			System.exit(-1);
		}

		// prepare a matrix for the genotypes/covariates
		double[][] x = null;
		if (cov != null) {
			x = new double[newSampleOrder.size()][cov.rows() + 1]; // [sample][genotype + covariates]
			// assume covariates are loaded as [covariate][sample]
			for (int i = 0; i < cov.rows(); i++) {
				for (int j = 0; j < cov.columns(); j++) {
					x[j][i + 1] = cov.getElement(i, j);
				}
			}
		} else {
			x = new double[newSampleOrder.size()][1];
		}

		System.out.println("Size of X: " + x.length + "x" + x[0].length);

		// some libs for the regression...
		cern.jet.random.tdouble.engine.DoubleRandomEngine randomEngine = null;
		cern.jet.random.tdouble.StudentT tDistColt = null;


		System.out.println("Starting calculations...");
		if (tDistColt == null) {
			randomEngine = new cern.jet.random.tdouble.engine.DRand();
			int df = x.length - (x[0].length + 1);
			System.out.println("Creating new tDist with " + df + " degrees of freedom..");
			tDistColt = new cern.jet.random.tdouble.StudentT(df, randomEngine); // df: nrsamples - nr parameters to estimate (== nr covariates + intercept)
		}

		ArrayList<Gene> tmpGeneList = new ArrayList<>();
		for (Gene g : finalGeneSetArr) {
			if (genesOverlappingVariants.contains(g)) {
				tmpGeneList.add(g);
			}
		}
		finalGeneSetArr = tmpGeneList;

		// threadpool
		ExecutorService threadPool = Executors.newFixedThreadPool(options.nrThreads);
		CompletionService<QTLOutput> jobHandler = new ExecutorCompletionService<QTLOutput>(threadPool);

		// iterate the available genes
		System.out.println(finalGeneSetArr.size() + " genes near SNPs for this VCF...");
		int submitted = 0;
		for (Gene g : finalGeneSetArr) {
			QTLTask task = new QTLTask(options, g, geneToExpGene, exp, variantSet, includeGenotypeSample, x, allGenotypes, tDistColt, permutationSeedNumber);
			jobHandler.submit(task);
			submitted++;
		}


		int[] genomeWideNull = new int[options.distsize];
		int[] genomeWideReal = new int[options.distsize];
		int[] topNull = new int[options.distsize];
		int[] topReal = new int[options.distsize];

		// output files
		TextFile outAll = new TextFile(outfile + "-all.txt.gz", TextFile.W);
		TextFile outTop = new TextFile(outfile + "-top.txt.gz", TextFile.W);
		String header = "GeneId" +
				"\tGeneSymbol" +
				"\tChr" +
				"\tStart" +
				"\tStop" +
				"\tStrand" +
				"\tSNP" +
				"\tPos" +
				"\tBeta" +
				"\tSe" +
				"\tZ" +
				"\tP" +
				"\tFractionLowerPvalObservedInPermForSNP" +
				"\tFractionPvalLowerThanLowestPermutedP" +
				"\tWithinGeneFDR";
		outAll.writeln(header);
		outTop.writeln(header);

		int received = 0;
		ProgressBar pb = new ProgressBar(submitted, "Running tasks with " + options.nrThreads + " threads");
		while (received < submitted) {
			try {
				Future<QTLOutput> future = jobHandler.take();
				if (future != null) {
					QTLOutput output = future.get();
					received++;
					if (output != null) {
						String[] textOut = output.getOutput();
						for (int i = 0; i < textOut.length; i++) {
							if (textOut[i] != null) {
								outAll.write(textOut[i]);
							}
						}
						if (output.getTopOut() != null) {
							outTop.write(output.getTopOut());
						}

						updateDist(genomeWideNull, output.getGeneWideNull());
						updateDist(genomeWideReal, output.getGeneWideReal());
						updateDist(topNull, output.getTopNull());
						updateDist(topReal, output.getTopReal());
					}
					pb.set(received);
				}

			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		pb.close();
		outAll.close();
		outTop.close();

		threadPool.shutdown();

		double[] genomewideFDR = calculateFDR(genomeWideReal, genomeWideNull);
		double[] topFDR = calculateFDR(topReal, topNull);

		TextFile distOut = new TextFile(outfile + "-dists.txt.gz", TextFile.W);
		distOut.writeln("Bin\tPval\tNReal\tNPerm\tFDR\tNRealTop\tNPermTop\tFDRTop");

		for (int i = 0; i < genomeWideReal.length; i++) {
			double p = Math.pow(10, -((double) i / 10000));
			if (i == 0) {
				p = 1;
			}
			if (i == genomeWideReal.length - 1) {
				p = 0;
			}

			distOut.writeln(i + "\t" + p
					+ "\t" + genomeWideReal[i]
					+ "\t" + genomeWideNull[i]
					+ "\t" + genomewideFDR[i]
					+ "\t" + topReal[i]
					+ "\t" + topNull[i]
					+ "\t" + topFDR[i]);
		}
		distOut.close();

		int binTop = 0;
		int binGW = 0;
		for (int i = genomeWideReal.length - 1; i > -1; i--) {
			if (genomewideFDR[i] >= 0.05) {
				binGW = i;
				break;
			}
		}
		for (int i = genomeWideReal.length - 1; i > -1; i--) {
			if (topFDR[i] >= 0.05) {
				binTop = i;
				break;
			}
		}

		System.out.println("Thresholds using per gene FDR: " + binToP(binTop));
		System.out.println("Thresholds using GW FDR: " + binToP(binGW));

	}

	private double binToP(int bin) {
		double p = Math.pow(10, -((double) bin / 10000));
		return p;
	}

	private void updateDist(int[] distIn, int[] update) {
		for (int i = 0; i < distIn.length; i++) {
			distIn[i] += update[i];
		}
	}

	private final Pair<Double, Double> NAN_PAIR = new Pair<Double, Double>(Double.NaN, Double.NaN);

	private Pair<Double, Double> convertBetaToP(double beta, double se, StudentT tDistColt) {

		if (Double.isNaN(beta)) {
			return NAN_PAIR;
		}

		double t = beta / se;
		double p = 1;
		double z = 0;
		if (t < 0) {
			p = tDistColt.cdf(t);
			if (p < 2.0E-323) {
				p = 2.0E-323;

			}
			z = cern.jet.stat.Probability.normalInverse(p);
		} else {
			p = tDistColt.cdf(-t);
			if (p < 2.0E-323) {
				p = 2.0E-323;

			}
			z = -cern.jet.stat.Probability.normalInverse(p);
		}
		return new Pair<Double, Double>(p, z);
	}


	class QTLOutput {

		private int[] geneWideReal;
		private int[] geneWideNull;

		private String[] output;

		private String topOut;

		private int[] topReal;
		private int[] topNull;

		public QTLOutput(int[] geneWideReal, int[] geneWideNull, String[] output, String topOut, int[] topReal, int[] topNull) {
			this.geneWideReal = geneWideReal;
			this.geneWideNull = geneWideNull;
			this.output = output;
			this.topOut = topOut;
			this.topReal = topReal;
			this.topNull = topNull;
		}

		public int[] getGeneWideReal() {
			return geneWideReal;
		}

		public int[] getGeneWideNull() {
			return geneWideNull;
		}

		public String[] getOutput() {
			return output;
		}

		public String getTopOut() {
			return topOut;
		}

		public int[] getTopReal() {
			return topReal;
		}

		public int[] getTopNull() {
			return topNull;
		}
	}

	class QTLTask implements Callable<QTLOutput> {


		private final QTLTestOptions options;
		private final Gene g;
		private final HashMap<Gene, Integer> geneToExpGene;
		private final DoubleMatrixDataset<String, String> exp;
		private final TreeSet<Feature> variantSet;
		private final boolean[] includeGenotypeSample;
		private final double[][] xOrig;
		private final LinkedHashMap<Feature, VCFVariant> allGenotypes;
		private final cern.jet.random.tdouble.StudentT tDistColt;

		private final long[] permutationSeedNumber;

		public QTLTask(QTLTestOptions options,
					   Gene g,
					   HashMap<Gene, Integer> geneToExpGene,
					   DoubleMatrixDataset<String, String> exp,
					   TreeSet<Feature> variantSet,
					   boolean[] includeGenotypeSample,
					   double[][] x,
					   LinkedHashMap<Feature, VCFVariant> allGenotypes,
					   cern.jet.random.tdouble.StudentT tDistColt,
					   long[] permutationSeedNumber) {
			this.options = options;
			this.g = g;
			this.geneToExpGene = geneToExpGene;
			this.exp = exp;
			this.variantSet = variantSet;
			this.includeGenotypeSample = includeGenotypeSample;
			this.xOrig = x;
			this.allGenotypes = allGenotypes;
			this.tDistColt = tDistColt;
			this.permutationSeedNumber = permutationSeedNumber;
		}

		@Override
		public QTLOutput call() throws Exception {

			int cisWindow = options.ciswindow;
			int nrPermutationsPerVariant = options.nrpermutationspergene;

			int[] geneWideReal = new int[options.distsize];
			int[] geneWideNull = new int[options.distsize];
			int[] topNull = new int[options.distsize];
			int[] topReal = new int[options.distsize];

			// get expression ID
			Integer expressionId = geneToExpGene.get(g);
			if (expressionId == null) {
				// all genes in the finalGeneSet should have at least one expression gene id
				System.out.println("Something weird is going on: " + g.getGeneId() + " not found in expression data");
				System.exit(-1);
			} else {

				StringBuilder geneHeader = new StringBuilder();
				geneHeader.append(g.getGeneId());
				geneHeader.append("\t").append(g.getName());
				geneHeader.append("\t").append(g.getChromosome().getName());
				geneHeader.append("\t").append(g.getStart());
				geneHeader.append("\t").append(g.getStop());
				geneHeader.append("\t").append(g.getStrand().getName());

				double[] y = exp.getMatrix().viewRow(expressionId).toArray();

				// -- get a list of variants within 1mb
				Strand strand = g.getStrand();
				Feature featureStart = new Feature(g);
				Feature featureStop = new Feature(g);
				if (strand.equals(Strand.NEG)) {
					// let's assume the coordinates of the genes are on the watson strand
					// even though the actual gene is on the crick strand (need to check though)
					featureStart.setStart(g.getStop() - cisWindow);
					featureStart.setStop(g.getStop() - cisWindow);
					featureStop.setStart(g.getStop() + cisWindow);
					featureStop.setStop(g.getStop() + cisWindow);
				} else {
					featureStart.setStart(g.getStart() - cisWindow);
					featureStart.setStop(g.getStart() - cisWindow);
					featureStop.setStart(g.getStart() + cisWindow);
					featureStop.setStop(g.getStart() + cisWindow);
				}
				SortedSet<Feature> overlappingVariants = variantSet.subSet(featureStart, true, featureStop, true);

				if (!overlappingVariants.isEmpty()) {

					String[] output = new String[overlappingVariants.size()];
					String topOut = null;

					// copy x
					double[][] x = new double[xOrig.length][xOrig[0].length];
					for (int i = 0; i < x.length; i++) {
						for (int j = 0; j < x[0].length; j++) {
							x[i][j] = xOrig[i][j];
						}
					}


					Feature[] overlappingVariantsArr = new Feature[overlappingVariants.size()];
					int fctr = 0;
					for (Feature f : overlappingVariants) {
						overlappingVariantsArr[fctr] = f;
						fctr++;
					}

					int[] numberOfTimesPvalueLowerThanLowestPermutedPval = new int[overlappingVariantsArr.length];
					int[] numberOfTimesPvalueLowerThanSNPPermutedPval = new int[overlappingVariantsArr.length];

					AssociationResult[] results = new AssociationResult[overlappingVariantsArr.length];
					// iterate the variants within the region
					VCFVariant topSNP = null;
					OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
					for (int permutation = -1; permutation < nrPermutationsPerVariant; permutation++) {
						double minP = 1;

						if (permutation > -1) {
							// permute Y
							// this is really ugly....
							ArrayList<Double> tmp = new ArrayList<>();
							for (int i = 0; i < y.length; i++) {
								tmp.add(y[i]);
							}
							Collections.shuffle(tmp, new Random(permutationSeedNumber[permutation])); // we can fix the seed per permutation
							y = Primitives.toPrimitiveArr(tmp.toArray(new Double[0]));
						}

						// for each variant
						// TODO: this can be made a lot faster by putting all genotype data into a matrix in permutation round -1
						for (int f = 0; f < overlappingVariantsArr.length; f++) {

							// get the genotypes
							VCFVariant variant = allGenotypes.get(overlappingVariantsArr[f]);

							// update the X matrix
							updateX(x, variant, includeGenotypeSample);

							// -- testNormal variant w/ linear model (use covariates)
							ols.newSampleData(y, x);

							double[] regressionParameters = ols.estimateRegressionParameters();
							double[] regressionStandardErrors = ols.estimateRegressionParametersStandardErrors();

							double betaSNP = regressionParameters[1];
							double seSNP = regressionStandardErrors[1];

							Pair<Double, Double> pair = convertBetaToP(betaSNP, seSNP, tDistColt);
							if (!pair.equals(NAN_PAIR)) {
								double pSNP = pair.getLeft();

								int bin = getBin(pSNP);


								if (permutation == -1) {
									results[f] = new AssociationResult();
									results[f].setBeta(new double[]{betaSNP});
									results[f].setSe(new double[]{seSNP});
									results[f].setPval(pSNP);
									geneWideReal[bin]++;
								} else {
									double realP = results[f].getPval();
									if (realP < pSNP) {
										numberOfTimesPvalueLowerThanSNPPermutedPval[f]++;
									}
									geneWideNull[bin]++;
								}

								// -- keep lowest pval, because meh
								if (pSNP < minP) {
									if (permutation == -1) {
										topSNP = variant;
									}
									minP = pSNP;
								}
							} else {
								if (permutation == -1) {
									results[f] = new AssociationResult();
									results[f].setBeta(new double[]{betaSNP});
									results[f].setSe(new double[]{seSNP});
									results[f].setPval(1);
								}
							}
						}

						for (int f = 0; f < overlappingVariantsArr.length; f++) {
							double p = results[f].getPval();
							if (p < minP) {
								numberOfTimesPvalueLowerThanLowestPermutedPval[f]++;
							}
						}

						int topBin = getBin(minP);
						if (permutation > -1) {
							topNull[topBin]++;
						} else {
							topReal[topBin]++;
						}
					}


					double[] fdr = calculateFDR(geneWideReal, geneWideNull);

					for (int f = 0; f < overlappingVariantsArr.length; f++) {
						// write the result (if it is not null)
						VCFVariant variant = allGenotypes.get(overlappingVariantsArr[f]);
						StringBuilder snpBuilder = new StringBuilder();
						snpBuilder.append(variant.getId());
						snpBuilder.append("\t").append(variant.getPos());
						snpBuilder.append("\t").append(results[f].getBeta()[0]);
						double z = ZScores.pToZ(results[f].getPval());
						if (results[f].getBeta()[0] < 0) {
							z *= -1;
						}
						snpBuilder.append("\t").append(results[f].getSe()[0]);
						snpBuilder.append("\t").append(z);
						snpBuilder.append("\t").append(results[f].getPval());
						double permutationP = (double) numberOfTimesPvalueLowerThanSNPPermutedPval[f] / nrPermutationsPerVariant;
						snpBuilder.append("\t").append(permutationP);

						permutationP = (double) numberOfTimesPvalueLowerThanLowestPermutedPval[f] / nrPermutationsPerVariant;
						snpBuilder.append("\t").append(permutationP);

						int bin = getBin(results[f].getPval());
						snpBuilder.append("\t").append(fdr[bin]);

						StringBuilder outBuilder = new StringBuilder();
						output[f] = outBuilder.append(geneHeader).append("\t").append(snpBuilder).append("\n").toString();

						if (variant.equals(topSNP)) {
							topOut = output[f];
						}
					}
					return new QTLOutput(geneWideReal, geneWideNull, output, topOut, topReal, topNull);
				}
			}
			return null;
		}


		private void updateX(double[][] x, VCFVariant variant, boolean[] includeGenotypeSample) {
			double[][] dosages = variant.getDosage(); // [samples][alleles]
			int ctr = 0;
			for (int i = 0; i < dosages.length; i++) {
				if (includeGenotypeSample[i]) {
					if (dosages[i][0] == -1) {
						System.err.println("ERROR: this method does currently not work with missing genotypes. Please impute the data.");
						System.exit(-1);
					}
					double dosage = dosages[i][0];
					x[ctr][0] = dosage; // x has structure [samples][covariates]
					ctr++;
				}
			}

		}
	}

	private double[] calculateFDR(int[] realDist, int[] nullDist) {
		// determine gene-wide FDR
		int[] cumulativeReal = new int[options.distsize];
		int sumReal = 0;
		int[] cumulativeNull = new int[options.distsize];
		int sumNull = 0;

		// these distributions are inverted..
		for (int i = options.distsize - 1; i > -1; i--) {
			if (i == options.distsize - 1) {
				cumulativeReal[i] = realDist[i];
				cumulativeNull[i] = nullDist[i];
			} else {
				cumulativeReal[i] = cumulativeReal[i + 1] + realDist[i];
				sumReal += realDist[i];
				cumulativeNull[i] = cumulativeNull[i + 1] + nullDist[i];
				sumNull += nullDist[i];
			}

		}

		// check whether the same number of tests have been performed
		if (sumNull / options.nrpermutationspergene != sumReal) {
			System.err.println("Error counting tests: " + sumNull + " null " + sumReal + " real");
			System.exit(-1);
		}

		double[] fdr = new double[options.distsize];
		for (int i = options.distsize - 1; i > -1; i--) {
			if (cumulativeReal[i] > 0 && cumulativeNull[i] > 0) {
				fdr[i] = ((double) cumulativeNull[i] / options.nrpermutationspergene) / cumulativeReal[i];
			} else if (cumulativeReal[i] > 0 && cumulativeNull[i] == 0) {
				fdr[i] = 0;
			}
		}

		return fdr;

	}

	private int getBin(double pSNP) {
		int bin = 0;
		if (pSNP == 1) {
			bin = 0;
		} else if (pSNP == 0) {
			bin = options.distsize - 1;
		} else {
			double log10p = -Math.log10(pSNP);
			if (log10p > 100) {
				log10p = 100;
			}
			bin = (int) Math.ceil(log10p * 10000);
			if (bin >= options.distsize) {
				bin = options.distsize - 1;
			} else if (bin < 0) {
				bin = 0;
			}
		}
		return bin;
	}


}
