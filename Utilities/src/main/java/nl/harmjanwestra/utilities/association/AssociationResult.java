package nl.harmjanwestra.utilities.association;

import nl.harmjanwestra.utilities.features.Feature;


/**
 * Created by hwestra on 11/3/15.
 */
public class AssociationResult {


	private Feature snp;
	private double or = Double.NaN;
	private double beta = Double.NaN;
	private double se = Double.NaN;
	private double pval = Double.NaN;
	private double maf = Double.NaN;

	private Feature region;
	private double bf = Double.NaN;
	private int n;
	private double devianceNull;
	private double devianceGeno;
	private int df;
	private double orLo;
	private double orHi;
	private double posterior;

	public Feature getSnp() {
		return snp;
	}

	public void setSnp(Feature snp) {
		this.snp = snp;
	}

	public double getOr() {
		return or;
	}

	public void setOr(double or) {
		this.or = or;
	}

	public double getBeta() {
		return beta;
	}

	public void setBeta(double beta) {
		this.beta = beta;
	}

	public double getSe() {
		return se;
	}

	public void setSe(double se) {
		this.se = se;
	}

	public double getPval() {
		return pval;
	}

	public void setPval(double pval) {
		this.pval = pval;
	}

	public double getMaf() {
		return maf;
	}

	public void setMaf(double maf) {
		this.maf = maf;
	}

	public Feature getRegion() {
		return region;
	}

	public void setRegion(Feature region) {
		this.region = region;
	}

	public double getBf() {
		return bf;
	}

	public void setBf(double bf) {
		this.bf = bf;
	}

	public int getN() {
		return n;
	}

	public void setN(int n) {
		this.n = n;
	}

	public double getDevianceNull() {
		return devianceNull;
	}

	public void setDevianceNull(double devianceNull) {
		this.devianceNull = devianceNull;
	}

	public double getDevianceGeno() {
		return devianceGeno;
	}

	public void setDevianceGeno(double devianceGeno) {
		this.devianceGeno = devianceGeno;
	}

	public int getDf() {
		return df;
	}

	public void setDf(int df) {
		this.df = df;
	}

	public double getOrLo() {
		return orLo;
	}

	public void setOrLo(double orLo) {
		this.orLo = orLo;
	}

	public double getOrHi() {
		return orHi;
	}

	public void setOrHi(double orHi) {
		this.orHi = orHi;
	}

	public double getPosterior() {
		return posterior;
	}

	public void setPosterior(double posterior) {
		this.posterior = posterior;
	}

	@Override
	public String toString() {
		// Chr     Pos     Id      CombinedId      N       MAF     DevianceNull    DevianceGeno    Df      Beta(Genotype)  SE(Genotype)    OR      OR-Hi   OR-Lo   Pval    -Log10(pval)
		String str = snp.getChromosome().toString()
				+ "\t" + snp.getStart()
				+ "\t" + snp.getName()
				+ "\t" + snp.toString()
				+ "\t" + n
				+ "\t" + maf
				+ "\t" + devianceNull
				+ "\t" + getDevianceGeno()
				+ "\t" + df
				+ "\t" + beta
				+ "\t" + se
				+ "\t" + or
				+ "\t" + orHi
				+ "\t" + orLo
				+ "\t" + pval
				+ "\t" + getLog10Pval();
		return str;

	}

	public double getLog10Pval() {
		return -Math.log10(pval);
	}
}
