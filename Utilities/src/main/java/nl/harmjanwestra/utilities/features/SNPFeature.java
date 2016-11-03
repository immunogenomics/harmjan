package nl.harmjanwestra.utilities.features;

import nl.harmjanwestra.utilities.enums.Chromosome;

/**
 * Created by hwestra on 11/11/15.
 */
public class SNPFeature extends Feature {
	double p;

	private double imputationQualScore;
	private String[] alleles;
	private String minorAllele;
	private double maf;
	private double hwep;
	private double cr;
	private double AFCases;
	private double AFControls;

	public double getMaf() {
		return maf;
	}

	public void setMaf(double maf) {
		this.maf = maf;
	}

	public double getHwep() {
		return hwep;
	}

	public void setHwep(double hwep) {
		this.hwep = hwep;
	}

	public double getCr() {
		return cr;
	}

	public void setCr(double cr) {
		this.cr = cr;
	}

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

	public void setAFCases(double AFCases) {
		this.AFCases = AFCases;
	}

	public void setAFControls(double AFControls) {
		this.AFControls = AFControls;
	}

	public double getAFCases() {
		return AFCases;
	}

	public double getAFControls() {
		return AFControls;
	}

	public static SNPFeature parseSNPFeature(String str) {

		String[] elems = str.split("_");
		if (elems.length == 3) {
			Chromosome chr = Chromosome.parseChr(elems[0]);

			Integer s1 = Integer.parseInt(elems[1]);
			String name = elems[2];

			SNPFeature out = new SNPFeature(chr, s1, s1);
			out.setName(name);

			return out;
		} else {
			return null;
		}
	}
}
