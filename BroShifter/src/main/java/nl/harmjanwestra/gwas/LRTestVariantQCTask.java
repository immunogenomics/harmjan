package nl.harmjanwestra.gwas;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.vcf.VCFImputationQualScoreBeagle;
import nl.harmjanwestra.utilities.vcf.VCFImputationQualScoreImpute;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;

import java.util.ArrayList;
import java.util.concurrent.Callable;

/**
 * Created by hwestra on 5/25/16.
 */
public class LRTestVariantQCTask implements Callable<Pair<VCFVariant, String>> {

	private final LRTestOptions options;
	private final boolean[] genotypeSamplesWithCovariatesAndDiseaseStatus;
	private final DoubleMatrix2D finalCovariates;
	private final DiseaseStatus[] finalDiseaseStatus;
	private String ln;
	private ArrayList<Feature> regions;

	public LRTestVariantQCTask(String ln,
							   ArrayList<Feature> regions,
							   LRTestOptions options,
							   boolean[] genotypeSamplesWithCovariatesAndDiseaseStatus,
							   DiseaseStatus[] finalDiseaseStatus,
							   DoubleMatrix2D finalCovariates) {
		this.regions = regions;
		this.options = options;
		this.ln = ln;
		this.genotypeSamplesWithCovariatesAndDiseaseStatus = genotypeSamplesWithCovariatesAndDiseaseStatus;
		this.finalDiseaseStatus = finalDiseaseStatus;
		this.finalCovariates = finalCovariates;
	}

	@Override
	public Pair<VCFVariant, String> call() throws Exception {

		LRTestTask tasktmp = new LRTestTask();
		VCFVariant variant = new VCFVariant(ln, VCFVariant.PARSE.HEADER);
		String variantChr = variant.getChr();
		String variantId = variant.getId();
		String variantPos = "" + variant.getPos();
		boolean variantPassesQC = false;
		double maf = 0;
		double hwep = 0;
		boolean overlap = false;
		Double impqual = 0d;
		if (!variantInRegion(variant)) {
			variant = null;
		} else {
			variant = new VCFVariant(ln, VCFVariant.PARSE.ALL);
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

			if (impqual == null || impqual > options.getImputationqualitythreshold()) {
				overlap = true;

				// parse the genotype, do some QC checks
				Triple<DoubleMatrix2D, boolean[], Triple<Integer, Double, Double>> unfilteredGenotypeData = tasktmp.filterAndRecodeGenotypes(
						genotypeSamplesWithCovariatesAndDiseaseStatus,
						variant.getGenotypeAllelesAsMatrix2D(),
						finalDiseaseStatus,
						variant.getAlleles().length,
						finalCovariates.rows());

				Triple<Integer, Double, Double> qc = unfilteredGenotypeData.getRight();

				maf = qc.getMiddle();
				hwep = qc.getRight();

				if (maf < options.getMafthresholdD() || hwep < options.getHWEPThreshold()) {
					variant = null;
				} else {
					variantPassesQC = true;
				}
			}
		}

		ln = null;
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
}
