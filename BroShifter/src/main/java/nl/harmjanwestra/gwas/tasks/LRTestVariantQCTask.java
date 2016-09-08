package nl.harmjanwestra.gwas.tasks;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.vcf.SampleAnnotation;
import nl.harmjanwestra.utilities.vcf.VCFImputationQualScoreBeagle;
import nl.harmjanwestra.utilities.vcf.VCFImputationQualScoreImpute;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.Callable;

/**
 * Created by hwestra on 5/25/16.
 */
public class LRTestVariantQCTask implements Callable<Pair<VCFVariant, String>> {

	private LRTestOptions options;
	private boolean[] genotypeSamplesWithCovariatesAndDiseaseStatus;
	private String ln;
	private ArrayList<Feature> regions;
	private HashSet<String> limitVariants;
	private SampleAnnotation sampleAnnotation;

	public LRTestVariantQCTask() {
	}

	public LRTestVariantQCTask(String ln,
							   ArrayList<Feature> regions,
							   LRTestOptions options,
							   boolean[] genotypeSamplesWithCovariatesAndDiseaseStatus,
							   SampleAnnotation sampleAnnotation,
							   HashSet<String> limitVariants) {
		this.regions = regions;
		this.options = options;
		this.ln = ln;
		this.genotypeSamplesWithCovariatesAndDiseaseStatus = genotypeSamplesWithCovariatesAndDiseaseStatus;
		this.sampleAnnotation = sampleAnnotation;
		this.limitVariants = limitVariants;
	}

	@Override
	public Pair<VCFVariant, String> call() throws Exception {

		VCFVariant variant = new VCFVariant(ln, VCFVariant.PARSE.HEADER);
		String variantChr = variant.getChr();
		String variantId = variant.getId();

		String varStr = variant.toString();


		String variantPos = "" + variant.getPos();
		boolean variantPassesQC = false;
		double maf = 0;
		double hwep = 0;
		boolean overlap = false;
		Double impqual = 0d;
		boolean includedinlist = false;
		if (limitVariants == null || limitVariants.contains(varStr)) {
			if (!variantInRegion(variant)) {
				variant = null;
				ln = null;
			} else {
				variant = new VCFVariant(ln, VCFVariant.PARSE.ALL, genotypeSamplesWithCovariatesAndDiseaseStatus, sampleAnnotation);
				if (!variant.hasImputationDosages()) {
					impqual = 1d;
				} else if (variant.getAlleles().length > 2) {
					VCFImputationQualScoreBeagle vbq = new VCFImputationQualScoreBeagle(variant, true);
					impqual = vbq.doseR2();
				} else {
					VCFImputationQualScoreImpute vsq = new VCFImputationQualScoreImpute();
					vsq.computeAutosomal(variant);
					impqual = vsq.getImpinfo();
				}
				overlap = true;

				if (impqual == null || impqual > options.getImputationqualitythreshold()) {
					// parse the genotype, do some QC checks
					if (variant.getMAFControls() < options.getMafthresholdD() || variant.getHwepControls() < options.getHWEPThreshold()) {
						variant = null;
					} else {
						variantPassesQC = true;
					}
					ln = null;
				} else {
					variant = null;
					ln = null;
				}
			}
		} else {
			variant = null;
			ln = null;
		}

		String logln = variantChr
				+ "\t" + variantPos
				+ "\t" + variantId
				+ "\t" + maf
				+ "\t" + hwep
				+ "\t" + impqual
				+ "\t" + overlap
				+ "\t" + variantPassesQC;

		return new Pair<VCFVariant, String>(variant, logln);
	}

	private boolean variantInRegion(VCFVariant variant) {
		if (regions == null) {
			return true;
		} else {
			Feature v = variant.asFeature();
			for (Feature f : regions) {
				if (f.overlaps(v)) {
					return true;
				}
			}
		}
		return false;
	}

	// recalculate MAF, HWE, etc, using CASE/Control labels
	public Triple<int[], boolean[], Integer> determineMissingGenotypes(
			DoubleMatrix2D genotypeAlleles,
			int nrsamples) {

		int individualCounter = 0;
		int nrWithMissingGenotypes = 0;
		boolean[] genotypeMissing = new boolean[nrsamples];
		for (int i = 0; i < genotypeAlleles.rows(); i++) {
			int b1 = (int) genotypeAlleles.getQuick(i, 0);
			if (b1 == -1 || Double.isNaN(b1)) {
				nrWithMissingGenotypes++;
				genotypeMissing[individualCounter] = true;
			}
			individualCounter++;
		}

		int[] missingGenotypeIds = new int[nrWithMissingGenotypes];
		int missingctr = 0;
		for (int i = 0; i < genotypeAlleles.rows(); i++) {
			if (genotypeMissing[i]) {
				missingGenotypeIds[missingctr] = i;
				missingctr++;
			}
		}
		return new Triple<>(missingGenotypeIds, genotypeMissing, nrWithMissingGenotypes);
	}
}
