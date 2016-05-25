package nl.harmjanwestra.utilities.association;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by hwestra on 10/27/15.
 */
public class AssociationFile {

	private boolean pairWise;

	public ArrayList<AssociationResult> readVariantPValues(String pvaluefile, Feature region) throws IOException {
		HashSet<String> variantHash = new HashSet<String>();
		TextFile textfile = new TextFile(pvaluefile, TextFile.R);

		ArrayList<AssociationResult> output = new ArrayList<AssociationResult>();

		// skip header
		textfile.readLine();
		String[] elems = textfile.readLineElems(TextFile.tab);
		int pvalctr = 0;
		while (elems != null) {
			// Marker	Chr	Position	PValue	Odds Ratio
			if (pvaluefile.endsWith(".tab")) {
				Chromosome chr = Chromosome.parseChr(elems[1]);
				if (region.getChromosome().equals(chr)) {
					String variant = elems[0] + "_" + elems[1] + "_" + elems[2];
					SNPFeature f2 = new SNPFeature();
					f2.setChromosome(chr);
					f2.setStart(Integer.parseInt(elems[2]));
					f2.setStop(Integer.parseInt(elems[2]));
					if (f2.overlaps(region)) {
						try {

							Double pval = Double.parseDouble(elems[elems.length - 2]); // need to check position...
							double log10 = -Math.log10(pval);
							variantHash.add(variant);

							f2.setName(elems[0]);
							AssociationResult r = new AssociationResult();
							r.setSnp(f2);
							r.setPval(pval);

							output.add(r);
						} catch (NumberFormatException e) {

						}
						pvalctr++;
					}

				}
			} else {
				Chromosome chr = Chromosome.parseChr(elems[1]);
				if (region.getChromosome().equals(chr)) {
					String variant = elems[0] + "_" + elems[1] + "_" + elems[2];
					SNPFeature f2 = new SNPFeature();
					f2.setChromosome(chr);
					Integer position = Integer.parseInt(elems[2]);
					f2.setStart(position);
					f2.setStop(position);
					if (f2.overlaps(region)) {

						Double pval = Double.parseDouble(elems[elems.length - 1]); // need to check position...
						variantHash.add(variant);

						AssociationResult r = new AssociationResult();
						r.setSnp(f2);
						r.setPval(pval);


						output.add(r);
						pvalctr++;
					}

				}
			}
			elems = textfile.readLineElems(TextFile.tab);
		}
		textfile.close();

		System.out.println(pvalctr + " pvalues for " + pvalctr + " positions from file: " + pvaluefile);

		return output;

	}

	private String model = null;

	public String getModel() {
		return model;
	}

	public ArrayList<AssociationResult> read(String file) throws IOException {
		return read(file, null);
	}

	public ArrayList<AssociationResult> read(String file, Feature region) throws IOException {
		TextFile tf = new TextFile(file, TextFile.R);
		String ln = tf.readLine();

		model = null;

		// Chr     Pos     Id      CombinedId      N       MAF     DevianceNull    DevianceGeno    Df      Beta(Genotype)  SE(Genotype)    OR      OR-Hi   OR-Lo   Pval    -Log10(pval)
		int chrcol = -1;
		int poscol = -1;
		int idcol = -1;
		int allelesCol = -1;
		int minorAlleleCol = -1;
		int combinedIdCol = -1;
		int ncol = -1;
		int mafcol = -1;
		int impqualscorecol = -1;
		int deviancenullcol = -1;
		int deviancegenocol = -1;
		int dfcol = -1;
		int betacol = -1;
		int secol = -1;
		int orcol = -1;
		int orhicol = -1;
		int orlocol = -1;
		int pvalcol = -1;
		int log10pvalcol = -1;
		int posteriorcol = -1;
		int bfcol = -1;
		int regioncol = -1;
		int hwepcol = -1;


		ArrayList<AssociationResult> results = new ArrayList<AssociationResult>();
		while (ln != null) {
			if (ln.startsWith("#Chromosome") || ln.startsWith("Chr\tPos")) {
// skip header
//				System.out.println("Found header");

				String[] elems = Strings.tab.split(ln);
				for (int i = 0; i < elems.length; i++) {
					String e = elems[i];
					if (e.equals("#Chromosome") || e.equals("Chr")) {
						chrcol = i;
					} else if (e.equals("Pos")) {
						poscol = i;
					} else if (e.equals("Id")) {
						idcol = i;
					} else if (e.equals("CombinedId")) {
						combinedIdCol = i;
					} else if (e.equals("Alleles")) {
						allelesCol = i;
					} else if (e.equals("MinorAllele")) {
						minorAlleleCol = i;
					} else if (e.equals("N")) {
						ncol = i;
					} else if (e.equals("MAF")) {
						mafcol = i;
					} else if (e.equals("HWEP")) {
						hwepcol = i;
					} else if (e.equals("DevianceNull")) {
						deviancenullcol = i;
					} else if (e.equals("DevianceGeno")) {
						deviancegenocol = i;
					} else if (e.equals("Df")) {
						dfcol = i;
					} else if (e.equals("Beta(Genotype)")) {
						betacol = i;
					} else if (e.equals("ImputationQualScore")) {
						impqualscorecol = i;
					} else if (e.equals("SE(Genotype)")) {
						secol = i;
					} else if (e.equals("OR")) {
						orcol = i;
					} else if (e.equals("OR-Hi")) {
						orhicol = i;
					} else if (e.equals("OR-Lo")) {
						orlocol = i;
					} else if (e.equals("Pval")) {
						pvalcol = i;
					} else if (e.equals("-Log10(pval)")) {
						log10pvalcol = i;
					} else if (e.equals("Region")) {
						regioncol = i;
					} else if (e.equals("BF")) {
						bfcol = i;
					} else if (e.equals("Posterior")) {
						posteriorcol = i;
					}
				}
			} else if (ln.startsWith("#")) {
				// get the model
				model = ln;
			} else {
				String[] elems = Strings.tab.split(ln);
				if (elems.length > 4) {
					// VariantID	N	MAF	DevianceNull	DfNull	DevianceGeno	DfAlt	Beta(Genotype)	SE(Genotype)	OR	OR-Hi	OR-Lo	Pval	-Log10(pval)
					Chromosome chr = Chromosome.NA;
					int pos = -1;
					String id = null;
					int n = 0;
					double maf = 0d;
					double hwep = 0d;
					double deviancenull = 0d;
					double deviancegeno = 0d;
					int df = 0;
					double[] beta = null;
					double[] se = null;

					double pval = 1d;

					double bf = 0d;
					double posterior = 0d;
					Feature assocregion = null;
					double impqualscore = Double.NaN;
					String[] alleles = null;
					String minorAllele = null;

					if (chrcol != -1) {
						chr = Chromosome.parseChr(elems[chrcol]);
					}
					if (poscol != -1) {
						try {
							pos = Integer.parseInt(elems[poscol]);
						} catch (NumberFormatException e) {

						}
					}
					if (idcol != -1) {
						id = elems[idcol];
					}

					if (allelesCol != -1) {
						String alleleStr = elems[allelesCol];
						String[] alleletmp = alleleStr.split(",");
						alleles = new String[alleletmp.length];
						for (int i = 0; i < alleletmp.length; i++) {
							alleles[i] = new String(alleletmp[i]).intern();
						}
					}
					if (minorAlleleCol != -1) {
						minorAllele = new String(elems[minorAlleleCol]).intern();
					}
					if (ncol != -1) {
						try {
							n = Integer.parseInt(elems[ncol]);
						} catch (NumberFormatException e) {

						}
					}
					if (mafcol != -1) {
						try {
							maf = Double.parseDouble(elems[mafcol]);
						} catch (NumberFormatException e) {

						}
					}

					if (hwepcol != -1) {
						try {
							hwep = Double.parseDouble(elems[hwepcol]);
						} catch (NumberFormatException e) {

						}
					}
					if (deviancenullcol != -1) {
						try {
							deviancenull = Double.parseDouble(elems[deviancenullcol]);
						} catch (NumberFormatException e) {

						}
					}
					if (deviancegenocol != -1) {
						try {
							deviancegeno = Double.parseDouble(elems[deviancegenocol]);
						} catch (NumberFormatException e) {

						}
					}
					if (dfcol != -1) {
						try {
							df = Integer.parseInt(elems[dfcol]);
						} catch (NumberFormatException e) {

						}

					}

					if (betacol != -1) {
						String betaStr = elems[betacol];
						String[] betaStrElems = betaStr.split(";");
						beta = new double[betaStrElems.length];
						for (int i = 0; i < betaStrElems.length; i++) {
							try {
								beta[i] = Double.parseDouble(betaStrElems[i]);
							} catch (NumberFormatException e) {

							}
						}
					}

					if (secol != -1) {
						String seStr = elems[secol];
						String[] seStrElems = seStr.split(";");
						se = new double[seStrElems.length];
						for (int i = 0; i < seStrElems.length; i++) {
							try {
								se[i] = Double.parseDouble(seStrElems[i]);
							} catch (NumberFormatException e) {

							}
						}
					}


					if (pvalcol != -1 && pvalcol < elems.length) {
						try {
							pval = Double.parseDouble(elems[pvalcol]);
						} catch (NumberFormatException e) {

						}
					}


					if (impqualscorecol != -1) {
						try {
							impqualscore = Double.parseDouble(elems[impqualscorecol]);
						} catch (NumberFormatException e) {

						}
					}

					if (regioncol != -1) {
						String regionStr = elems[regioncol];
						assocregion = Feature.parseFeature(regionStr);
					}

					if (bfcol != -1) {
						try {
							bf = Double.parseDouble(elems[bfcol]);
						} catch (NumberFormatException e) {

						}
					}
					if (posteriorcol != -1) {
						try {
							posterior = Double.parseDouble(elems[posteriorcol]);
						} catch (NumberFormatException e) {

						}
					}


					SNPFeature snp = new SNPFeature(chr, pos, pos);
					snp.setName(id);
					snp.setImputationQualityScore(impqualscore);
					snp.setAlleles(alleles);
					snp.setMinorAllele(minorAllele);

					if (region == null || region.overlaps(snp)) {
						AssociationResult result = new AssociationResult();
						result.setSnp(snp);
						result.setN(n);
						result.setMaf(maf);
						result.setDevianceNull(deviancenull);
						result.setDevianceGeno(deviancegeno);
						result.setDf(df);
						result.setBeta(beta);
						result.setSe(se);
						result.setPval(pval);
						result.setBf(bf);
						result.setHWEP(hwep);
						result.setPosterior(posterior);
						result.setRegion(assocregion);

						results.add(result);
					}
				}
			}

			ln = tf.readLine();
		}
		tf.close();
		return results;
	}

	public String getHeader() {
		String str = "";
		if (model != null) {
			str = model + "\n";
		}

		if (pairWise) {
			str += "#Chromosome" +
					"\tPos1" +
					"\tId1" +
					"\tPos2" +
					"\tId2" +
					"\tCombinedId1" +
					"\tCombinedId2" +
					"\tAlleles1" +
					"\tAlleles2" +
					"\tMinorAllele1" +
					"\tMinorAllele2" +
					"\tImputationQualScore" +
					"\tImputationQualScore2" +
					"\tN" +
					"\tMAF1" +
					"\tMAF2" +
					"\tHWEP" +
					"\tHWEP2" +
					"\tLD(RSQ)" +
					"\tLD(D')" +
					"\tDevianceNull" +
					"\tDevianceGeno" +
					"\tDfAlt" +
					"\tBeta(Genotype)" +
					"\tSE(Genotype)" +
					"\tOR" +
					"\tOR-Hi" +
					"\tOR-Lo" +
					"\tPval" +
					"\t-Log10(pval)";
		} else {
			str += "#Chromosome" +
					"\tPos" +
					"\tId" +
					"\tCombinedId" +
					"\tAlleles" +
					"\tMinorAllele" +
					"\tImputationQualScore" +
					"\tN" +
					"\tMAF" +
					"\tHWEP" +
					"\tDevianceNull" +
					"\tDevianceGeno" +
					"\tDfAlt" +
					"\tBeta(Genotype)" +
					"\tSE(Genotype)" +
					"\tOR" +
					"\tOR-Hi" +
					"\tOR-Lo" +
					"\tPval" +
					"\t-Log10(pval)";
		}

		return str;
	}

	public void setPairWise(boolean pairWise) {
		this.pairWise = pairWise;
	}
}
