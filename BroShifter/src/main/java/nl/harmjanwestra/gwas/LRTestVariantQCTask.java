package nl.harmjanwestra.gwas;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.vcf.VCFImputationQualScoreBeagle;
import nl.harmjanwestra.utilities.vcf.VCFImputationQualScoreImpute;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.math.stats.HWE;

import java.util.ArrayList;
import java.util.concurrent.Callable;

/**
 * Created by hwestra on 5/25/16.
 */
public class LRTestVariantQCTask implements Callable<Pair<VCFVariant, String>> {

	private LRTestOptions options;
	private boolean[] genotypeSamplesWithCovariatesAndDiseaseStatus;
	private DoubleMatrix2D finalCovariates;
	private DiseaseStatus[] finalDiseaseStatus;
	private String ln;
	private ArrayList<Feature> regions;

	public LRTestVariantQCTask() {
	}

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
			ln = null;
		} else {
			variant = new VCFVariant(ln, VCFVariant.PARSE.ALL, genotypeSamplesWithCovariatesAndDiseaseStatus);
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
				Triple<int[], boolean[], Triple<Integer, Double, Double>> unfilteredGenotypeData = filterAndRecodeGenotypes(
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
				ln = null;
			} else {
				variant = null;
				ln = null;
			}
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
	Triple<int[], boolean[], Triple<Integer, Double, Double>> filterAndRecodeGenotypes(
			DoubleMatrix2D genotypeAlleles,
			DiseaseStatus[] diseaseStatus,
			int nrAlleles,
			int nrsamples) {

		int individualCounter = 0;

		// first iterate the genotyped samples to load the genotypes
		// calculate maf and hwep

		// initialize hwep stuff
		int nrCombinations = (nrAlleles * (nrAlleles + 1)) / 2;
		int[] obs = new int[nrCombinations];
		int[] nrHomozygous = new int[nrAlleles];
		int nrCalled = 0;
		int[][] index = new int[nrAlleles][nrAlleles];

		int ctr = 0;
		for (int i = 0; i < nrAlleles; i++) {
			for (int j = i; j < nrAlleles; j++) {
				index[i][j] = ctr;
				index[j][i] = ctr;
				ctr++;
			}
		}

		int nrWithMissingGenotypes = 0;
		boolean[] genotypeMissing = new boolean[nrsamples];
		for (int i = 0; i < genotypeAlleles.rows(); i++) {
			int b1 = (int) genotypeAlleles.getQuick(i, 0);
			if (b1 == -1) {
				nrWithMissingGenotypes++;
				genotypeMissing[individualCounter] = true;
			} else {
				if (diseaseStatus[individualCounter].equals(DiseaseStatus.CONTROL)) { // controls
					int b2 = (int) genotypeAlleles.getQuick(i, 1);
					if (b1 == b2) {
						nrHomozygous[b1]++;
					}

					int id = index[b1][b2];
					obs[id]++;
					nrCalled++;
				}
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


		double[] freqs = new double[nrAlleles];
		for (int i = 0; i < nrAlleles; i++) {
			for (int j = 0; j < nrAlleles; j++) {
				int id = index[i][j];
				if (i == j) {
					freqs[i] += (2 * obs[id]);
				} else {
					freqs[i] += obs[id];
				}
			}
		}
		for (int i = 0; i < nrAlleles; i++) {
			freqs[i] /= (nrCalled * 2);
		}

		double hwep = 0;
		boolean debug = true;
		if (nrAlleles == 2 && !debug) {
			hwep = HWE.calculateExactHWEPValue(obs[1], obs[0], obs[2]);
		} else {
			ctr = 0;
			double chisq = 0;
			for (int i = 0; i < nrAlleles; i++) {
				for (int j = i; j < nrAlleles; j++) {
					double expectedFreq;
					if (i == j) {
						expectedFreq = (freqs[i] * freqs[i]) * nrCalled; // homozygote freq
					} else {
						expectedFreq = (2 * (freqs[i] * freqs[j])) * nrCalled; // heterozygote freq
					}
					double observedFreq = obs[ctr];
					double oe = (observedFreq - expectedFreq);
					if (oe != 0 && expectedFreq != 0) {
						chisq += ((oe * oe) / expectedFreq);
					}
					ctr++;
				}
			}

			int df = (nrCombinations - nrAlleles);
			hwep = umcg.genetica.math.stats.ChiSquare.getP(df, chisq);
		}

		double maf = 1;
		for (int i = 0; i < freqs.length; i++) {
			double d = freqs[i];
			if (d < maf) {
				maf = d;
			}
		}


		return new Triple<>(missingGenotypeIds, genotypeMissing, new Triple<>(nrWithMissingGenotypes, maf, hwep));
	}
}
