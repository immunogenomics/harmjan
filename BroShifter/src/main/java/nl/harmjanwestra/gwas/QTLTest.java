package nl.harmjanwestra.gwas;

import cern.jet.random.tdouble.StudentT;
import nl.harmjanwestra.gwas.CLI.QTLTestOptions;
import nl.harmjanwestra.utilities.association.AssociationResult;
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
		Long permutationSeedNumber = random.nextLong(); // To think about: do we want a fixed order in our permutations? if so, each permutation should get a unique seed, shared between each gene
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
		GTFAnnotation geneAnnotation = new GTFAnnotation(gtfAnnotationFile);
		Collection<Gene> allGenes = geneAnnotation.getGenes();


		ArrayList<Gene> finalGeneSetArr = new ArrayList<>();
		LinkedHashMap<String, Integer> expGeneHash = exp.getHashRows();
		LinkedHashMap<Gene, Integer> geneToExpGene = new LinkedHashMap<Gene, Integer>();

		HashSet<String> geneNamesFound = new HashSet<>();
		for (Gene gene : allGenes) {
			String name = gene.getGeneId();
			if (name.contains("ENSG00000238009")) {
				System.out.println("found it: ENSG00000238009.2");
			}
			if (expGeneHash.containsKey(name)) {
				Integer id = expGeneHash.get(name);

				finalGeneSetArr.add(gene);
				geneToExpGene.put(gene, id);
				geneNamesFound.add(name);
			}
		}

		TreeSet<Gene> finalGeneSet = new TreeSet<>(new FeatureComparator(false));
		finalGeneSet.addAll(finalGeneSetArr);

		TextFile notfound = new TextFile(outfile + "-expressionNotFound.txt", TextFile.W);
		ArrayList<String> expgenes = exp.getRowObjects();
		for (String gene : expgenes) {
			if (!geneNamesFound.contains(gene)) {
				notfound.writeln(gene);
			}
		}
		notfound.close();

		if (geneToExpGene.isEmpty()) {
			System.err.println("Error: no genes found. Probably a mismatch between annotation and gene expression file");
			System.exit(-1);
		}
		System.out.println(geneToExpGene.size() + " genes with annotation found in gene expression data");
		// rank the gene expression data -- save for later

		// for each chromosome
		// load all genotypes (load only biallelic markers)
		// put the variants in a treemap
		TreeSet<Feature> variantSet = new TreeSet<Feature>(new FeatureComparator(false));
		LinkedHashMap<Feature, VCFVariant> allGenotypes = new LinkedHashMap<>();
		System.out.println("Loading variants");
		HashSet<Gene> genesOverlappingVariants = new HashSet<Gene>();
		int variantctr = 0;
		while (genotypeVCF.hasNext()) {
			VCFVariant variant = genotypeVCF.next();
			Double impqual = variant.getInfo().get("AR2");
			if (impqual == null || impqual > imputationqualthreshold) {
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
		OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
		cern.jet.random.tdouble.engine.DoubleRandomEngine randomEngine = null;
		cern.jet.random.tdouble.StudentT tDistColt = null;

		// output files
		TextFile outAll = new TextFile(outfile + "-all.txt", TextFile.W);
		TextFile outTop = new TextFile(outfile + "-top.txt", TextFile.W);
		String header = "Gene\tHugo\tChr\tStart\tStop\tStrand\tSNP\tPos\tBeta\tSe\tZ\tP\tPermutationP";
		outAll.writeln(header);
		outTop.writeln(header);

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

		// iterate the available genes
		ProgressBar pb = new ProgressBar(finalGeneSetArr.size());
		int genectr = 0;
		for (Gene g : finalGeneSetArr) {
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

				System.out.println(genectr + "/" + finalGeneSetArr.size() + "\t" + g.getGeneId() + "\t" + overlappingVariants.size());
				if (!overlappingVariants.isEmpty()) {
					Feature[] overlappingVariantsArr = new Feature[overlappingVariants.size()];
					int fctr = 0;
					for (Feature f : overlappingVariants) {
						overlappingVariantsArr[fctr] = f;
						fctr++;
					}


					int[] nrPValsLowerInPermutation = new int[overlappingVariantsArr.length];
					AssociationResult[] results = new AssociationResult[overlappingVariantsArr.length];

					// iterate the variants within the region
					VCFVariant topSNP = null;
//					double[] topEffectsPerPermutation = new double[nrPermutationsPerVariant];


					for (int permutation = -1; permutation < nrPermutationsPerVariant; permutation++) {
						double minP = 1;

						if (permutation > -1) {
							// permute Y
							// this is really ugly....
							ArrayList<Double> tmp = new ArrayList<>();
							for (int i = 0; i < y.length; i++) {
								tmp.add(y[i]);
							}
							Collections.shuffle(tmp, new Random(permutationSeedNumber)); // we can fix the seed per permutation
							y = Primitives.toPrimitiveArr(tmp.toArray(new Double[0]));
						}

						// for each variant
						// TODO: this can be made a lot faster by putting all genotype data into a matrix in permutation round -1
						for (int f = 0; f < overlappingVariantsArr.length; f++) {

							// get the genotypes
							VCFVariant variant = allGenotypes.get(overlappingVariantsArr[f]);

							// update the X matrix
							updateX(x, variant, includeGenotypeSample);

							// -- test variant w/ linear model (use covariates)
							ols.newSampleData(y, x);

							double[] regressionParameters = ols.estimateRegressionParameters();
							double[] regressionStandardErrors = ols.estimateRegressionParametersStandardErrors();

							double betaSNP = regressionParameters[1];
							double seSNP = regressionStandardErrors[1];

							Pair<Double, Double> pair = convertBetaToP(betaSNP, seSNP, tDistColt);
							if (!pair.equals(NAN_PAIR)) {
								double pSNP = pair.getLeft();

								// create a result object
								if (permutation == -1) {
									results[f] = new AssociationResult();
									results[f].setBeta(new double[]{betaSNP});
									results[f].setSe(new double[]{seSNP});
									results[f].setPval(pSNP);
								}

								// -- keep lowest pval, because meh
								if (pSNP < minP) {
									if (permutation == -1) {
										topSNP = variant;
									}
//								else {
//									topEffectsPerPermutation[permutation] = minP;
//								}
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
								nrPValsLowerInPermutation[f]++;
							}
						}
						System.out.print("\rPermutation: " + permutation);
					}
					System.out.println();

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
						double permutationP = (double) nrPValsLowerInPermutation[f] / nrPermutationsPerVariant;
						snpBuilder.append("\t").append(permutationP);
						outAll.append(geneHeader);
						outAll.append("\t");
						outAll.append(snpBuilder);
						outAll.append("\n");

						if (variant.equals(topSNP)) {
							outTop.append(geneHeader);
							outTop.append("\t");
							outTop.append(snpBuilder);
							outTop.append("\n");
						}
					}

//					// write the top effect if there is any
//					if (topSNP != null) {
//						// permutations are done.
//						// determine likelihood of lowest pval
//						Arrays.sort(topEffectsPerPermutation);
//						int nrLower = 0;
//						for (int i = 0; i < nrPermutationsPerVariant; i++) {
//							if (topResult.getPval() < topEffectsPerPermutation[i]) {
//								nrLower++;
//							} else {
//								break;
//							}
//						}
//						double permutationP = (double) nrLower / nrPermutationsPerVariant;
//
//						// output top effect for this gene.
//						StringBuilder snpBuilder = new StringBuilder();
//						snpBuilder.append(topSNP.getId());
//						snpBuilder.append("\t").append(topSNP.getPos());
//						snpBuilder.append("\t").append(topResult.getBeta());
//						snpBuilder.append("\t").append(topResult.getSe());
//						snpBuilder.append("\t").append(ZScores.pToZ(topResult.getPval()));
//						snpBuilder.append("\t").append(topResult.getPval());
//						snpBuilder.append("\t").append(permutationP);
//
//						outTop.append(geneHeader);
//						outTop.append("\t");
//						outTop.append(snpBuilder);
//						outTop.append("\n");
//					}
				}
			}
			pb.iterate();
			genectr++;
		}
		pb.close();
		outAll.close();
		outTop.close();
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

	private void updateX(double[][] x, VCFVariant variant, boolean[] includeGenotypeSample) {
		double[][] dosages = variant.getImputedDosages(); // [samples][alleles]
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
