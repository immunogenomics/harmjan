package nl.harmjanwestra.utilities.association;

import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import umcg.genetica.text.Strings;


/**
 * Created by hwestra on 11/3/15.
 */
public class AssociationResult {


	protected SNPFeature snp;

	protected double[] beta = null;
	protected double[] se = null;
	protected double pval = Double.NaN;

	protected Feature region;
	protected double bf = Double.NaN;
	protected int n;
	protected double devianceNull;
	protected double devianceGeno;
	protected int dfnull;
	protected int dfalt;
	protected int df;
	protected double posterior;

	protected double[] OR;
	protected AssociationResult[] subresults;
	private double z;


	public SNPFeature getSnp() {
		return snp;
	}

	public void setSnp(SNPFeature snp) {
		this.snp = snp;
	}

	public double[] getBeta() {
		return beta;
	}

	public void setBeta(double[] beta) {
		this.beta = beta;
	}

	public double[] getSe() {
		return se;
	}

	public void setSe(double[] se) {
		this.se = se;
	}

	public double getPval() {
		return pval;
	}

	public void setPval(double pval) {
		this.pval = pval;
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

	public void setDfnull(int dfnull) {
		this.dfnull = dfnull;
	}

	public void setDfalt(int dfalt) {
		this.dfalt = dfalt;
	}

	public double getPosterior() {
		return posterior;
	}

	public void setPosterior(double posterior) {
		this.posterior = posterior;
	}

	@Override
	public String toString() {
		String str = snp.getChromosome().toString()
				+ "\t" + snp.getStart()
				+ "\t" + snp.getName()
				+ "\t" + snp.toString()
				+ "\t" + Strings.concat(snp.getAlleles(), Strings.comma)
				+ "\t" + snp.getMinorAllele()
				+ "\t" + snp.getImputationQualityScore()
				+ "\t" + n
				+ "\t" + snp.getMaf()
				+ "\t" + Strings.concat(snp.getAFCasesArray(), Strings.comma)
				+ "\t" + Strings.concat(snp.getAFControlsArray(), Strings.comma)
				+ "\t" + snp.getHwep();


		// Chr     Pos     Id      CombinedId      N       MAF     DevianceNull    DevianceGeno    Df      Beta(Genotype)  SE(Genotype)    OR      OR-Hi   OR-Lo   Pval    -Log10(pval)
		str += "\t" + devianceNull
				+ "\t" + devianceGeno
				+ "\t" + dfnull
				+ "\t" + dfalt
				+ "\t" + df;

		if (beta != null) {
			str += "\t" + Strings.concat(beta, Strings.semicolon);
		} else {
			str += "\tnull";
		}
		if (se != null) {
			str += "\t" + Strings.concat(se, Strings.semicolon);
		} else {
			str += "\tnull";
		}

		if (beta != null && se != null) {
			str += "\t" + Strings.concat(getORs(), Strings.semicolon)
					+ "\t" + Strings.concat(getConfHi(), Strings.semicolon)
					+ "\t" + Strings.concat(getConfLo(), Strings.semicolon);
		} else {
			str += "\tnull\tnull\tnull";
		}


		str += "\t" + pval
				+ "\t" + getLog10Pval();
		return str;

	}

	public double getLog10Pval() {
		if (pval == 0d || Double.isNaN(pval) || Double.isInfinite(pval)) {
			return 0;
		}
		return -Math.log10(pval);
	}

	public double[] getConfHi() {
		if (beta != null && se != null) {
			double[] output = new double[beta.length];
			for (int i = 0; i < output.length; i++) {
				output[i] = Math.exp(beta[i] + 1.96 * se[i]);
			}
			return output;
		} else {
			return null;
		}
	}

	public double[] getConfLo() {
		if (beta != null && se != null) {
			double[] output = new double[beta.length];
			for (int i = 0; i < output.length; i++) {
				output[i] = Math.exp(beta[i] - 1.96 * se[i]);
			}
			return output;
		} else {
			return null;
		}
	}

	public double[] getORs() {
		if (beta != null && se != null) {
			double[] output = new double[beta.length];
			for (int i = 0; i < output.length; i++) {
				output[i] = Math.exp(beta[i]);
			}
			return output;
		} else if (OR != null) {
			return OR;
		} else {
			return null;
		}
	}


	public void setOR(double[] OR) {
		this.OR = OR;
	}

	public void setSubresults(AssociationResult[] subresults) {
		this.subresults = subresults;
	}

	public AssociationResult[] getSubresults() {
		return subresults;
	}


	public double getZ() {
		return beta[0] / se[0];
	}
}
