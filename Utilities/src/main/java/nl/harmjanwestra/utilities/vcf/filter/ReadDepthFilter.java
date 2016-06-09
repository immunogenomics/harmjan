package nl.harmjanwestra.utilities.vcf.filter;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.utilities.vcf.VCFVariant;

/**
 * Created by hwestra on 2/8/16.
 */
public class ReadDepthFilter implements VCFGenotypeFilter {


	private int minimalReadDepth = 10;

	public ReadDepthFilter() {

	}

	public ReadDepthFilter(int minimalReadDepth) {
		this.minimalReadDepth = minimalReadDepth;
	}

	public void filter(VCFVariant variant) {

		int nrSamples = variant.getNrSamples();
		DoubleMatrix2D alleles = variant.getGenotypeAllelesAsMatrix2D();
		short[][] allelicDepth = variant.getAllelicDepth();
		short[] approximateDepth = variant.getApproximateDepth();

		if (approximateDepth != null) {
			if (allelicDepth != null) {
				// account for a bug in newer GATK output
				for (int i = 0; i < nrSamples; i++) {
					short indSum = 0;
					for (int j = 0; j < allelicDepth.length; j++) {
						indSum += allelicDepth[j][i];
					}
					if (approximateDepth[i] != indSum) {
						approximateDepth[i] = indSum;
					}
				}
			}

			for (int i = 0; i < nrSamples; i++) {
				int depth = approximateDepth[i];
				if (depth < minimalReadDepth) {
					alleles.setQuick(i, 0, -1);
					alleles.setQuick(i, 1, -1);
				}
			}
		} else {
			// should everything be -1 now???
		}
	}
}
