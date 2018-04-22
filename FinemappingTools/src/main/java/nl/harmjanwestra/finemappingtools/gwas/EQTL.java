package nl.harmjanwestra.finemappingtools.gwas;

import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;

/**
 * Created by hwestra on 10/28/16.
 */
public class EQTL extends Feature {

	double pval;
	double qval;
	double beta;
	double zscore;
	String genename;
	SNPFeature snp;

	public double getPval() {
		return pval;
	}

	public void setPval(double pval) {
		this.pval = pval;
	}

	public double getQval() {
		return qval;
	}

	public void setQval(double qval) {
		this.qval = qval;
	}

	public double getBeta() {
		return beta;
	}

	public void setBeta(double beta) {
		this.beta = beta;
	}

	public double getZscore() {
		return zscore;
	}

	public void setZscore(double zscore) {
		this.zscore = zscore;
	}

	public String getGenename() {
		return genename;
	}

	public void setGenename(String genename) {
		this.genename = genename;
	}

	public SNPFeature getSnp() {
		return snp;
	}

	public void setSnp(SNPFeature snp) {
		this.snp = snp;
	}

	@Override
	public String toString() {
		return "EQTL{" +
				"genename='" + genename + '\'' +
				", snp=" + snp +
				'}';
	}
}
