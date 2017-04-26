package nl.harmjanwestra.finemappingtools.gwas;

import nl.harmjanwestra.finemappingtools.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.vcf.SampleAnnotation;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.vcf.VCFVariantComparator;
import nl.harmjanwestra.utilities.vcf.VCFVariantLoader;
import nl.harmjanwestra.utilities.vcf.filter.variantfilters.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.ExecutorService;

/**
 * Created by hwestra on 3/22/17.
 */
public class VariantLoader {

	public ArrayList<VCFVariant> load(String logfile,
									  ArrayList<Feature> regions,
									  LRTestOptions options,
									  boolean[] genotypeSamplesWithCovariatesAndDiseaseStatus,
									  SampleAnnotation sampleAnnotation,
									  ExecutorService exService

	) throws IOException {

		VCFVariantLoader loader = new VCFVariantLoader(genotypeSamplesWithCovariatesAndDiseaseStatus, sampleAnnotation, exService);
		VCFVariantFilters filter = new VCFVariantFilters();

		// there is a file that limits the snps to include
		if (options.getSnpLimitFile() != null) {
			filter.addFilter(new VCFVariantSetFilter(options.getSnpLimitFile()));
		}

		filter.addFilter(new VCFVariantCallRateFilter(options.getCallrateThreshold()));
		filter.addFilter(new VCFVariantImpQualFilter(options.getImputationqualitythreshold(), true));
		filter.addFilter(new VCFVariantMAFFilter(options.getMafthresholdD(), VCFVariantMAFFilter.MODE.CONTROLS));
		filter.addFilter(new VCFVariantHWEPFilter(options.getHWEPThreshold(), VCFVariantHWEPFilter.MODE.CONTROLS));
		filter.addFilter(new VCFVariantRegionFilter(regions));

		ArrayList<VCFVariant> variants = loader.run(options.getVcf(), filter);
		Collections.sort(variants, new VCFVariantComparator());


//		if (snpLimit != null) {
//			System.out.println("Limiting to " + snpLimit.size() + " variants.");
//		} else {
//			System.out.println("Not using snplimit.");
//		}

		return variants;

	}

	// original code
//	public ArrayList<VCFVariant> load(String logfile, ArrayList<Feature> regions, LRTestOptions options, ExecutorService exService) throws IOException {
//		TextFile log = new TextFile(logfile, TextFile.W);
//		System.out.println("Log will be written here: " + log.getFileName());
//		log.writeln("Chr\tPos\tId\tMAF\tHWEP\tImpQual\tOverlap\tPassQC");
//
//		int maxQueueSize = options.getNrThreads() * 10;
//		CompletionService<Pair<VCFVariant, String>> jobHandler = new ExecutorCompletionService<Pair<VCFVariant, String>>(exService);
//
//		TextFile vcfIn = new TextFile(options.getVcf(), TextFile.R);
//		String vcfLn = vcfIn.readLine();
//		while (vcfLn != null && vcfLn.startsWith("#")) {
//			vcfLn = vcfIn.readLine();
//		}
//		int linessubmitted = 0;
//		System.out.println("Parsing: " + options.getVcf());
//
//
//		int totalsubmitted = 0;
//		ArrayList<VCFVariant> variants = new ArrayList<>();
//		int read = 0;
//		while (vcfLn != null) {
//			if (linessubmitted == maxQueueSize) {
//				// get some of the output
//				int returned = 0;
//				while (returned < linessubmitted) {
//					try {
//						Future<Pair<VCFVariant, String>> future = jobHandler.take();
//						if (future != null) {
//							Pair<VCFVariant, String> result = future.get();
//
//							if (log != null) {
//								String logStr = result.getRight();
//								log.writeln(logStr);
//							}
//							if (result.getLeft() != null) {
//								VCFVariant variant = result.getLeft();
//								int nrAlleles = variant.getNrAlleles();
//								if (options.splitMultiAllelicIntoMultipleVariants && nrAlleles > 2) {
//									for (int a = 1; a < nrAlleles; a++) {
//										variants.add(variant.variantFromAllele(a));
//									}
//								} else {
//									variants.add(variant);
//								}
//							}
//						}
//						returned++;
//
//					} catch (InterruptedException e) {
//						e.printStackTrace();
//					} catch (ExecutionException e) {
//						e.printStackTrace();
//					}
//
//				}
//				linessubmitted = 0;
//			}
//
//
//			boolean submit = false;
//			if (snpLimit != null) {
//				VCFVariant variantF = new VCFVariant(vcfLn, VCFVariant.PARSE.HEADER);
//				String snpStr = Chromosome.parseChr(variantF.getChr()) + "_" + variantF.getPos() + "_" + variantF.getId();
//				if (snpLimit.contains(snpStr)) {
//					submit = true;
//					System.out.println("Found one");
//				}
//			} else {
//				if (regions != null) {
//					Feature variantF = new VCFVariant(vcfLn, VCFVariant.PARSE.HEADER).asFeature();
//					for (Feature region : regions) {
//						if (variantF.overlaps(region)) {
//							submit = true;
//						}
//					}
//
//				} else {
//					submit = true;
//				}
//			}
//
//
//			if (submit) {
//				LRTestVariantQCTask task = new LRTestVariantQCTask(vcfLn,
//						bedRegions,
//						options,
//						genotypeSamplesWithCovariatesAndDiseaseStatus,
//						sampleAnnotation,
//						null);
//				jobHandler.submit(task);
//				linessubmitted++;
//				totalsubmitted++;
//			}
//			read++;
//			if (read % 1000 == 0) {
//				System.out.print(read + " lines read\t" + totalsubmitted + " submitted. " + variants.size() + " variants in memory " + ManagementFactory.getThreadMXBean().getThreadCount() + " threads\r");
//			}
//			vcfLn = vcfIn.readLine();
//		}
//		System.out.println();
//		System.out.println(totalsubmitted + " total lines read.");
//		vcfIn.close();
//
//		int returned = 0;
//		while (returned < linessubmitted) {
//			try {
//				Future<Pair<VCFVariant, String>> future = jobHandler.take();
//				if (future != null) {
//					Pair<VCFVariant, String> result = future.get();
//
//					if (log != null) {
//						String logStr = result.getRight();
//						log.writeln(logStr);
//					}
//					if (result.getLeft() != null) {
//						VCFVariant variant = result.getLeft();
//						int nrAlleles = variant.getNrAlleles();
//						if (options.splitMultiAllelicIntoMultipleVariants && nrAlleles > 2) {
//							System.out.println("Splitting alleles");
//
//							for (int a = 1; a < nrAlleles; a++) {
//								variants.add(variant.variantFromAllele(a));
//							}
//						} else {
//							variants.add(variant);
//						}
//					}
//				}
//				returned++;
//
//			} catch (InterruptedException e) {
//				e.printStackTrace();
//			} catch (ExecutionException e) {
//				e.printStackTrace();
//			}
//		}
//		System.out.println(variants.size() + "  variants finally included.");
//
//		log.close();
//
//		// why not sort :)
//
//		return variants;
//	}
//
//	public class LRTestVariantQCTask implements Callable<Pair<VCFVariant, String>> {
//
//		private LRTestOptions options;
//		private boolean[] genotypeSamplesWithCovariatesAndDiseaseStatus;
//		private String ln;
//		private ArrayList<Feature> regions;
//		private HashSet<String> limitVariants;
//		private SampleAnnotation sampleAnnotation;
//
//		public LRTestVariantQCTask() {
//		}
//
//		public LRTestVariantQCTask(String ln,
//								   ArrayList<Feature> regions,
//								   LRTestOptions options,
//								   boolean[] genotypeSamplesWithCovariatesAndDiseaseStatus,
//								   SampleAnnotation sampleAnnotation,
//								   HashSet<String> limitVariants) {
//			this.regions = regions;
//			this.options = options;
//			this.ln = ln;
//			this.genotypeSamplesWithCovariatesAndDiseaseStatus = genotypeSamplesWithCovariatesAndDiseaseStatus;
//			this.sampleAnnotation = sampleAnnotation;
//			this.limitVariants = limitVariants;
//		}
//
//		@Override
//		public Pair<VCFVariant, String> call() throws Exception {
//
//			VCFVariant variant = new VCFVariant(ln, VCFVariant.PARSE.HEADER);
//			String variantChr = variant.getChr();
//			String variantId = variant.getId();
//
//			String varStr = variant.toString();
//
//
//			String variantPos = "" + variant.getPos();
//			boolean variantPassesQC = false;
//			double maf = 0;
//			double cr = 0;
//			double hwep = 0;
//			boolean overlap = false;
//			Double impqual = 0d;
//			boolean includedinlist = false;
//			if (limitVariants == null || limitVariants.contains(varStr)) {
//				if (!variantInRegion(variant)) {
//					variant = null;
//					ln = null;
//				} else {
//					variant = new VCFVariant(ln, VCFVariant.PARSE.ALL, genotypeSamplesWithCovariatesAndDiseaseStatus, sampleAnnotation);
//					if (!variant.hasImputationProbabilities()) {
//						impqual = 1d;
//					} else if (variant.getAlleles().length > 2) {
//						VCFImputationQualScoreBeagle vbq = new VCFImputationQualScoreBeagle(variant, true);
//						impqual = vbq.doseR2();
//					} else {
//						VCFImputationQualScoreImpute vsq = new VCFImputationQualScoreImpute();
//						vsq.computeAutosomal(variant);
//						impqual = vsq.getImpinfo();
//					}
//					overlap = true;
//
//					maf = variant.getMAFControls();
//					hwep = variant.getHwepControls();
//					cr = variant.getCallrate();
//
//
//					if (impqual == null || impqual > options.getImputationqualitythreshold()) {
//						// parse the genotype, do some QC checks
//						if (variant.getMAFControls() < options.getMafthresholdD() || variant.getHwepControls() < options.getHWEPThreshold() || variant.getCallrate() < options.getCallrateThreshold()) {
//							variant = null;
//						} else {
//							variantPassesQC = true;
//						}
//
//						ln = null;
//					} else {
//						variant = null;
//						ln = null;
//					}
//				}
//			} else {
//				variant = null;
//				ln = null;
//			}
//
//			String logln = variantChr
//					+ "\t" + variantPos
//					+ "\t" + variantId
//					+ "\t" + maf
//					+ "\t" + hwep
//					+ "\t" + cr
//					+ "\t" + impqual
//					+ "\t" + overlap
//					+ "\t" + variantPassesQC;
//
//			return new Pair<VCFVariant, String>(variant, logln);
//		}
//
//		private boolean variantInRegion(VCFVariant variant) {
//			if (regions == null) {
//				return true;
//			} else {
//				Feature v = variant.asFeature();
//				for (Feature f : regions) {
//					if (f.overlaps(v)) {
//						return true;
//					}
//				}
//			}
//			return false;
//		}
//
//		// recalculate MAF, HWE, etc, using CASE/Control labels
//		public Triple<int[], boolean[], Integer> determineMissingGenotypes(
//				DoubleMatrix2D genotypeAlleles,
//				int nrsamples) {
//
//			int individualCounter = 0;
//			int nrWithMissingGenotypes = 0;
//			boolean[] genotypeMissing = new boolean[nrsamples];
//			for (int i = 0; i < genotypeAlleles.rows(); i++) {
//				int b1 = (int) genotypeAlleles.getQuick(i, 0);
//				if (b1 == -1 || Double.isNaN(b1)) {
//					nrWithMissingGenotypes++;
//					genotypeMissing[individualCounter] = true;
//				}
//				individualCounter++;
//			}
//
//			int[] missingGenotypeIds = new int[nrWithMissingGenotypes];
//			int missingctr = 0;
//			for (int i = 0; i < genotypeAlleles.rows(); i++) {
//				if (genotypeMissing[i]) {
//					missingGenotypeIds[missingctr] = i;
//					missingctr++;
//				}
//			}
//			return new Triple<>(missingGenotypeIds, genotypeMissing, nrWithMissingGenotypes);
//		}
//	}

}
