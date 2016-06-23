package nl.harmjanwestra.utilities.association;

import nl.harmjanwestra.utilities.features.SNPFeature;
import umcg.genetica.text.Strings;

/**
 * Created by hwestra on 6/22/16.
 */
public class AssociationResultPairwise extends AssociationResult {


	private SNPFeature snp2;
	private double ldRSquared;
	private Double ldDprime;

	public void setSnp2(SNPFeature snp2) {
		this.snp2 = snp2;
	}

	public double getLdRSquared() {
		return ldRSquared;
	}

	public SNPFeature getSnp2() {
		return snp2;
	}

	public void setLDRSquared(double LD) {
		this.ldRSquared = LD;
	}

	public void setLdDprime(Double ldDprime) {
		this.ldDprime = ldDprime;
	}

	@Override
	public String toString() {

		String str = snp.getChromosome().toString()
				+ "\t" + snp.getStart()
				+ "\t" + snp.getName()
				+ "\t" + snp.toString()
				+ "\t" + snp2.getChromosome().toString()
				+ "\t" + snp2.getStart()
				+ "\t" + snp2.getName()
				+ "\t" + snp2.toString()
				+ "\t" + Strings.concat(snp.getAlleles(), Strings.comma)
				+ "\t" + Strings.concat(snp2.getAlleles(), Strings.comma)
				+ "\t" + snp.getMinorAllele()
				+ "\t" + snp2.getMinorAllele()
				+ "\t" + snp.getImputationQualityScore()
				+ "\t" + snp2.getImputationQualityScore()
				+ "\t" + n
				+ "\t" + snp.getMaf()
				+ "\t" + snp2.getMaf()
				+ "\t" + snp.getHwep()
				+ "\t" + snp2.getHwep()
				+ "\t" + (snp.getStart() - snp2.getStart())
				+ "\t" + ldRSquared
				+ "\t" + ldDprime;
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
			str += "\tnull";
		}

		str += "\t" + pval
				+ "\t" + getLog10Pval();
		return str;
	}


}
