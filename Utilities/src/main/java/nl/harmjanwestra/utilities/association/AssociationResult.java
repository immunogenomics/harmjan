package nl.harmjanwestra.utilities.association;

import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import umcg.genetica.text.Strings;


/**
 * Created by hwestra on 11/3/15.
 */
public class AssociationResult {


	private SNPFeature snp;
	private double[] beta = null;
	private double[] se = null;
	private double pval = Double.NaN;
	private double maf = Double.NaN;
	private Feature region;
	private double bf = Double.NaN;
	private int n;
	private double devianceNull;
	private double devianceGeno;
	private int df;
	private double posterior;

	public Feature getSnp() {
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
				+ "\t" + Strings.concat(snp.getAlleles(), Strings.comma)
				+ "\t" + snp.getMinorAllele()
				+ "\t" + snp.getImputationQualityScore()
				+ "\t" + n
				+ "\t" + maf
				+ "\t" + devianceNull
				+ "\t" + devianceGeno
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
			str += "\tnull";
		}

		str += "\t" + pval
				+ "\t" + getLog10Pval();
		return str;

	}

	public double getLog10Pval() {
		if(pval == 0d || Double.isNaN(pval) || Double.isInfinite(pval)){
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
		} else {
			return null;
		}
	}


}
