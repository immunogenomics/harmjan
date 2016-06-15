package nl.harmjanwestra.utilities.features;

/**
 * Created by hwestra on 11/11/15.
 */
public class SNPFeature extends Feature {
	double p;

	private double imputationQualScore;
	private String[] alleles;
	private String minorAllele;

	public SNPFeature() {

	}

	public SNPFeature(SNPFeature f2) {
		super(f2);
		this.p = f2.getP();
	}

	public SNPFeature(Chromosome chr, int start, int end) {
		super(chr, start, end);
	}

	public double getImputationQualityScore() {
		return imputationQualScore;
	}

	public void setImputationQualityScore(double imputationQualScore) {
		this.imputationQualScore = imputationQualScore;
	}

	public String[] getAlleles() {
		return alleles;
	}

	public void setAlleles(String[] alleles) {
		this.alleles = alleles;
	}

	public String getMinorAllele() {
		return minorAllele;
	}

	public void setMinorAllele(String minorAllele) {
		this.minorAllele = minorAllele;
	}

	public double getP() {
		return p;
	}

	public void setP(double p) {
		this.p = p;
	}

	@Override
	public String toString() {
		return getChromosome().getName() + "_" + getStart() + "_" + name;
	}
}
