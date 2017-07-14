package nl.harmjanwestra.utilities.association;

import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.ZScores;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;


/**
 * Created by hwestra on 11/3/15.
 */
public class AssociationResult {
	
	
	protected SNPFeature snp;
	
	protected double[][] beta = null; // [disease][allele]
	protected double[][] se = null; // [disease][allele]
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
	
	protected double[][] OR; // [disease][allele]
	protected AssociationResult[] subresults;
	private double z;
	
	
	public SNPFeature getSnp() {
		return snp;
	}
	
	public void setSnp(SNPFeature snp) {
		this.snp = snp;
	}
	
	public double[][] getBeta() {
		return beta;
	}
	
	public void setBeta(double[][] beta) {
		this.beta = beta;
	}
	
	public double[][] getSe() {
		return se;
	}
	
	public void setSe(double[][] se) {
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
				+ "\t" + snp.getCr()
				+ "\t" + snp.getMissingnessP()
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
			
			// diisease1Allele2,disease1Allele3;diisease2Allele2,disease2Allele3
			// format: beta[disease][allele]
			
			str += "\t" + mergeMultiStr(beta);
		} else {
			str += "\tnull";
		}
		if (se != null) {
			str += "\t" + mergeMultiStr(se);
		} else {
			str += "\tnull";
		}
		
		if (beta != null && se != null) {
			
			str += "\t" + mergeMultiStr(getORs())
					+ "\t" + mergeMultiStr(getConf(true))
					+ "\t" + mergeMultiStr(getConf(false));
		} else {
			str += "\tnull\tnull\tnull";
		}
		
		
		str += "\t" + pval
				+ "\t" + getLog10Pval();
		return str;
		
	}
	
	protected String mergeMultiStr(double[][] input) {
		String tmpStr = "";
		for (int i = 0; i < input.length; i++) {
			if (i == 0) {
				tmpStr += Strings.concat(input[i], Strings.comma);
			} else {
				tmpStr += ";" + Strings.concat(input[i], Strings.comma);
			}
		}
		return tmpStr;
	}
	
	
	public double getLog10Pval() {
		if (pval == 0d || Double.isNaN(pval) || Double.isInfinite(pval)) {
			return 0;
		}
		return -Math.log10(pval);
	}
	
	public double[][] getConf(boolean high) {
		if (beta != null && se != null) {
			double[][] output = new double[beta.length][beta[0].length];
			for (int i = 0; i < output.length; i++) {
				for(int j=0;j<output[i].length;j++) {
					if(high) {
						output[i][j] = Math.exp(beta[i][j] + 1.96 * se[i][j]);
					} else {
						output[i][j] = Math.exp(beta[i][j] - 1.96 * se[i][j]);
					}
				}
			}
			return output;
		} else {
			return null;
		}
	}
	
	public double[][] getORs() {
		if (beta != null && se != null) {
			double[][] output = new double[beta.length][beta[0].length];
			for (int i = 0; i < output.length; i++) {
				for (int j = 0; j < output[i].length; j++) {
					output[i][j] = Math.exp(beta[i][j]);
				}
			}
			return output;
		} else if (OR != null) {
			return OR;
		} else {
			return null;
		}
	}
	
	
	public void setOR(double[][] OR) {
		this.OR = OR;
	}
	
	public void setSubresults(AssociationResult[] subresults) {
		this.subresults = subresults;
	}
	
	public AssociationResult[] getSubresults() {
		return subresults;
	}
	
	public double getZFromP() {
		return ZScores.pToZ(pval);
	}
	
	public double getZ() {
		return ZScores.betaToZ(beta[0][0], se[0][0]);
	}
}
