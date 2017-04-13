package nl.harmjanwestra.utilities.association;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 6/22/16.
 */
public class AssociationFilePairwise {


	public String getHeader() {
		String str = "#Chr1" +
				"\tPos1" +
				"\tId1" +
				"\tCombinedId1" +
				"\tChr2" +
				"\tPos2" +
				"\tId2" +
				"\tCombinedId2" +
				"\tAlleles1" +
				"\tAlleles2" +
				"\tMinorAllele1" +
				"\tMinorAllele2" +
				"\tImputationQualScore" +
				"\tImputationQualScore2" +
				"\tN" +
				"\tMAF1" +
				"\tAFCases1" +
				"\tAFControls1" +
				"\tMAF2" +
				"\tAFCases2" +
				"\tAFControls2" +
				"\tHWEP" +
				"\tHWEP2" +
				"\tDistance" +
				"\tLD(RSQ)" +
				"\tLD(D')" +
				"\tDevianceNull" +
				"\tDevianceGeno" +
				"\tDfNull" +
				"\tDfAlt" +
				"\tDiffDf" +
				"\tBeta(Genotype)" +
				"\tSE(Genotype)" +
				"\tOR" +
				"\tOR-Hi" +
				"\tOR-Lo" +
				"\tPval" +
				"\tLog10(p)";
		return str;
	}

	public ArrayList<AssociationResultPairwise> read(String file, Feature region) throws IOException {

		ArrayList<AssociationResultPairwise> output = new ArrayList<>();

		int chr1col = -1;
		int pos1col = -1;
		int Id1col = -1;
		int CombinedId1col = -1;
		int Chr2col = -1;
		int Pos2col = -1;
		int Id2col = -1;
		int CombinedId2col = -1;
		int Alleles1col = -1;
		int Alleles2col = -1;
		int MinorAllele1col = -1;
		int MinorAllele2col = -1;
		int ImputationQualScorecol = -1;
		int ImputationQualScore2col = -1;
		int Ncol = -1;
		int MAF1col = -1;
		int MAF2col = -1;
		int AFCases1Col = -1;
		int AFCases2Col = -1;
		int AFControls1Col = -1;
		int AFControls2Col = -1;
		int HWEPcol = -1;
		int HWEP2col = -1;
		int Distancecol = -1;
		int LDRSQcol = -1;
		int LDDcol = -1;
		int DevianceNullcol = -1;
		int DevianceGenocol = -1;
		int DfNullcol = -1;
		int DfAltcol = -1;
		int DiffDfcol = -1;
		int BetaGenotypecol = -1;
		int SEGenotypecol = -1;
		int ORcol = -1;
		int ORHicol = -1;
		int ORLocol = -1;
		int Pvalcol = -1;
		int Log10pcol = -1;

		System.out.println("Reading assoc path: " + file);
		TextFile tf = new TextFile(file, TextFile.R);
		String ln = tf.readLine();

		while (ln != null) {
			String[] elems = ln.split("\t");
			if (ln.startsWith("#Chr1")) {
				for (int i = 0; i < elems.length; i++) {
					switch (elems[i]) {
						case "#Chr1":
							chr1col = i;
							break;
						case "Pos1":
							pos1col = i;
							break;
						case "Id1":
							Id1col = i;
							break;
						case "CombinedId1":
							CombinedId1col = i;
							break;
						case "Chr2":
							Chr2col = i;
							break;
						case "Pos2":
							Pos2col = i;
							break;
						case "Id2":
							Id2col = i;
							break;
						case "CombinedId2":
							CombinedId2col = i;
							break;
						case "Alleles1":
							Alleles1col = i;
							break;
						case "Alleles2":
							Alleles2col = i;
							break;
						case "MinorAllele1":
							MinorAllele1col = i;
							break;
						case "MinorAllele2":
							MinorAllele2col = i;
							break;
						case "ImputationQualScore":
							ImputationQualScorecol = i;
							break;
						case "ImputationQualScore2":
							ImputationQualScore2col = i;
							break;
						case "N":
							Ncol = i;
							break;
						case "MAF1":
							MAF1col = i;
							break;
						case "MAF2":
							MAF2col = i;
							break;
						case "HWEP":
							HWEPcol = i;
							break;
						case "HWEP2":
							HWEP2col = i;
							break;
						case "Distance":
							Distancecol = i;
							break;
						case "LD(RSQ)":
							LDRSQcol = i;
							break;
						case "LD(D')":
							LDDcol = i;
							break;
						case "DevianceNull":
							DevianceNullcol = i;
							break;
						case "DevianceGeno":
							DevianceGenocol = i;
							break;
						case "DfNull":
							DfNullcol = i;
							break;
						case "DfAlt":
							DfAltcol = i;
							break;
						case "DiffDf":
							DiffDfcol = i;
							break;
						case "Beta(Genotype)":
							BetaGenotypecol = i;
							break;
						case "SE(Genotype)":
							SEGenotypecol = i;
							break;
						case "OR":
							ORcol = i;
							break;
						case "OR-Hi":
							ORHicol = i;
							break;
						case "OR-Lo":
							ORLocol = i;
							break;
						case "Pval":
							Pvalcol = i;
							break;
						case "Log10(p)":
							Log10pcol = i;
							break;


						case "AFCases1":
							AFCases1Col = i;
							break;
						case "AFControls1":
							AFControls1Col = i;
							break;
						case "AFCases2":
							AFCases2Col = i;
							break;
						case "AFControls2":
							AFControls2Col = i;
							break;
					}
				}
			} else {
				if (elems.length > 4) {

					AssociationResultPairwise r = new AssociationResultPairwise();

					SNPFeature snp1 = new SNPFeature();
					SNPFeature snp2 = new SNPFeature();

					if (chr1col != -1) {
						snp1.setChromosome(Chromosome.parseChr(elems[chr1col]));
					}
					if (pos1col != -1) {
						snp1.setStart(Integer.parseInt(elems[pos1col]));
						snp1.setStop(snp1.getStart());
					}
					if (Id1col != -1) {
						snp1.setName(elems[Id1col]);
					}

					if (Chr2col != -1) {
						snp2.setChromosome(Chromosome.parseChr(elems[Chr2col]));
					}
					if (Pos2col != -1) {
						snp2.setStart(Integer.parseInt(elems[Pos2col]));
						snp2.setStop(snp1.getStart());
					}
					if (Id2col != -1) {
						snp2.setName(elems[Id2col]);
					}

					if (Alleles1col != -1) {
						snp1.setAlleles(Strings.comma.split(elems[Alleles1col]));
					}
					if (Alleles2col != -1) {
						snp2.setAlleles(Strings.comma.split(elems[Alleles2col]));
					}

					if (MinorAllele1col != -1) {
						snp1.setMinorAllele(elems[MinorAllele1col]);
					}
					if (MinorAllele2col != -1) {
						snp2.setMinorAllele(elems[MinorAllele1col]);
					}

					if (ImputationQualScorecol != -1) {
						snp1.setImputationQualityScore(Double.parseDouble(elems[ImputationQualScorecol]));
					}
					if (ImputationQualScore2col != -1) {
						snp2.setImputationQualityScore(Double.parseDouble(elems[ImputationQualScore2col]));
					}

					if (Ncol != -1) {
						r.setN(Integer.parseInt(elems[Ncol]));
					}
					if (MAF1col != -1) {
						snp1.setMaf(Double.parseDouble(elems[MAF1col]));
					}
					if (MAF2col != -1) {
						snp2.setMaf(Double.parseDouble(elems[MAF2col]));
					}

					if (AFCases1Col != -1) {
						snp1.setAFCases(Double.parseDouble(elems[AFCases1Col]));
					}
					if (AFControls1Col != -1) {
						snp1.setAFControls(Double.parseDouble(elems[AFControls1Col]));
					}
					if (AFCases2Col != -1) {
						snp2.setAFCases(Double.parseDouble(elems[AFCases2Col]));
					}
					if (AFControls2Col != -1) {
						snp2.setAFControls(Double.parseDouble(elems[AFControls2Col]));
					}

					if (HWEPcol != -1) {
						snp1.setHwep(Double.parseDouble(elems[HWEPcol]));
					}
					if (HWEP2col != -1) {
						snp2.setHwep(Double.parseDouble(elems[HWEP2col]));
					}

					if (LDRSQcol != -1) {
						r.setLDRSquared(Double.parseDouble(elems[LDRSQcol]));
					}
					if (LDDcol != -1) {
						r.setLDRSquared(Double.parseDouble(elems[LDDcol]));
					}
					if (DevianceNullcol != -1) {
						r.setDevianceNull(Double.parseDouble(elems[DevianceNullcol]));
					}
					if (DevianceGenocol != -1) {
						r.setDevianceGeno(Double.parseDouble(elems[DevianceGenocol]));
					}
					if (DfNullcol != -1) {
						r.setDfnull(Integer.parseInt(elems[DfNullcol]));
					}
					if (DfAltcol != -1) {
						r.setDfalt(Integer.parseInt(elems[DfAltcol]));
					}
					if (DiffDfcol != -1) {
						r.setDf(Integer.parseInt(elems[DiffDfcol]));
					}

					if (BetaGenotypecol != -1) {
						String[] subelems = Strings.semicolon.split(elems[BetaGenotypecol]);
						double[] beta = new double[subelems.length];
						for (int e = 0; e < subelems.length; e++) {
							beta[e] = Double.parseDouble(subelems[e]);
						}
						r.setBeta(beta);
					}

					if (SEGenotypecol != -1) {
						String[] subelems = Strings.semicolon.split(elems[SEGenotypecol]);
						double[] se = new double[subelems.length];
						for (int e = 0; e < subelems.length; e++) {
							se[e] = Double.parseDouble(subelems[e]);
						}
						r.setSe(se);
					}

					if (Pvalcol != -1) {
						r.setPval(Double.parseDouble(elems[Pvalcol]));
					}

					r.setSnp(snp1);
					r.setSnp(snp2);
					if (region == null || (region.overlaps(snp1) && region.overlaps(snp2))) {
						output.add(r);
					}
				}
			}
		}

		tf.close();
		return output;
	}

}
