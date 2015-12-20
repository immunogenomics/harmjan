package nl.harmjanwestra.utilities.association;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;


/**
 * Created by hwestra on 11/3/15.
 */
public class AssociationResult {
	Feature snp;
	double or = Double.NaN;
	double beta = Double.NaN;
	double se = Double.NaN;
	double pval = Double.NaN;
	double maf = Double.NaN;
	double bf = Double.NaN;
	double abf = Double.NaN;

	public void setBf(double bf) {
		this.bf = bf;
	}

	public void setAbf(double abf) {
		this.abf = abf;
	}

	public double getBf() {
		return bf;
	}

	public void setSnp(Feature snp) {
		this.snp = snp;
	}

	public void setOr(double or) {
		this.or = or;
	}

	public void setBeta(double beta) {
		this.beta = beta;
	}

	public void setSe(double se) {
		this.se = se;
	}

	public void setPval(double pval) {
		this.pval = pval;
	}

	public void setMaf(double maf) {
		this.maf = maf;
	}

	public double getAbf() {
		return abf;
	}

	public AssociationResult() {

	}

	public Feature getSnp() {
		return snp;
	}

	public double getOr() {
		return or;
	}

	public double getBeta() {
		return beta;
	}

	public double getSe() {
		return se;
	}

	public double getPval() {
		return pval;
	}

	public double getMaf() {
		return maf;
	}

	public AssociationResult(Chromosome chr, int pos, double maf, double or, double beta, double se, double pval) {
		snp = new Feature(chr, pos, pos);
		this.maf = maf;

		this.or = or;
		this.beta = beta;
		this.se = se;
		this.pval = pval;
	}

	public void setLocation(Chromosome chr, int pos) {
		snp = new Feature(chr, pos, pos);
	}

	public void setAssociation(double or, double beta, double se, double pval) {
		this.or = or;
		this.beta = beta;
		this.se = se;
		this.pval = pval;
	}

}
