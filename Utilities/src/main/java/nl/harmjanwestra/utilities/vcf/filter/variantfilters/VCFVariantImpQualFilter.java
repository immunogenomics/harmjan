package nl.harmjanwestra.utilities.vcf.filter.variantfilters;

import nl.harmjanwestra.utilities.vcf.VCFImputationQualScoreBeagle;
import nl.harmjanwestra.utilities.vcf.VCFImputationQualScoreImpute;
import nl.harmjanwestra.utilities.vcf.VCFVariant;

/**
 * Created by hwestra on 3/21/17.
 */
public class VCFVariantImpQualFilter implements VCFVariantFilter {
	private double threshold = 0.3;
	private boolean recalculate = false;
	
	public VCFVariantImpQualFilter(double threshold, boolean recalculate) {
		this.threshold = threshold;
		this.recalculate = recalculate;
	}
	
	public VCFVariantImpQualFilter(double threshold) {
		this.threshold = threshold;
	}
	
	@Override
	public boolean passesThreshold(VCFVariant variant) {
		
		Double impqual = variant.getImputationQualityScore();
		if (recalculate) {
			if (!variant.hasImputationProbabilities()) {
				impqual = 1d;
				variant.setImputationQualityScore(impqual, true);
			} else if (variant.getAlleles().length > 2) {
				VCFImputationQualScoreBeagle vbq = new VCFImputationQualScoreBeagle(variant, true);
				impqual = vbq.doseR2();
				variant.setImputationQualityScore(impqual, false);
			} else {
				VCFImputationQualScoreImpute vsq = new VCFImputationQualScoreImpute();
				vsq.computeAutosomal(variant);
				impqual = vsq.getImpinfo();
				variant.setImputationQualityScore(impqual, false);
			}
			
		}
		
		if (impqual != null) {
			return impqual > threshold;
		} else {
			return true;
		}
		
	}
	
	@Override
	public String toString() {
		return "VCFVariantImpQualFilter{" +
				"threshold=" + threshold +
				", recalculate=" + recalculate +
				'}';
	}
}
