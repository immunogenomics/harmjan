package nl.harmjanwestra.utilities.association;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by hwestra on 10/27/15.
 */
public class AssociationFile {

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
					Feature f2 = new Feature();
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
							r.setPval(log10);
							r.setOr(Double.parseDouble(elems[4]));
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
					Feature f2 = new Feature();
					f2.setChromosome(chr);
					f2.setStart(Integer.parseInt(elems[2]));
					f2.setStop(Integer.parseInt(elems[2]));
					if (f2.overlaps(region)) {
						Integer position = Integer.parseInt(elems[2]);
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

	public ArrayList<AssociationResult> loadConditionalAssocData(String file) throws IOException {
		return loadConditionalAssocData(file, null);
	}

	public ArrayList<AssociationResult> loadConditionalAssocData(String file, Feature region) throws IOException {
		TextFile tf = new TextFile(file, TextFile.R);
		String ln = tf.readLine();

		model = null;

		// Chr     Pos     Id      CombinedId      N       MAF     DevianceNull    DevianceGeno    Df      Beta(Genotype)  SE(Genotype)    OR      OR-Hi   OR-Lo   Pval    -Log10(pval)
		int chrcol = -1;
		int poscol = -1;
		int idcol = -1;
		int combinedIdCol = -1;
		int ncol = -1;
		int mafcol = -1;
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


		ArrayList<AssociationResult> results = new ArrayList<AssociationResult>();
		while (ln != null) {
			if (ln.startsWith("#")) {
				// get the model
				model = ln;
			} else if (ln.startsWith("VariantID")) {
// skip header

				String[] elems = Strings.tab.split(ln);
				for (int i = 0; i < elems.length; i++) {
					String e = elems[i];
					if (e.equals("Chr")) {
						chrcol = i;
					} else if (e.equals("Pos")) {
						poscol = i;
					} else if (e.equals("Id")) {
						idcol = i;
					} else if (e.equals("CombinedId")) {
						combinedIdCol = i;
					} else if (e.equals("N")) {
						ncol = i;
					} else if (e.equals("MAF")) {
						mafcol = i;
					} else if (e.equals("DevianceNull")) {
						deviancenullcol = i;
					} else if (e.equals("DevianceGeno")) {
						deviancegenocol = i;
					} else if (e.equals("Df")) {
						dfcol = i;
					} else if (e.equals("Beta(Genotype)")) {
						betacol = i;
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

			} else {
				String[] elems = Strings.tab.split(ln);
				if (elems.length > 4) {
					// VariantID	N	MAF	DevianceNull	DfNull	DevianceGeno	DfAlt	Beta(Genotype)	SE(Genotype)	OR	OR-Hi	OR-Lo	Pval	-Log10(pval)


					Chromosome chr = Chromosome.NA;
					int pos = -1;
					String id = null;
					String combinedId = null;
					int n = 0;
					double maf = 0d;
					double deviancenull = 0d;
					double deviancegeno = 0d;
					int df = 0;
					double beta = 0d;
					double se = 0d;
					double or = 0d;
					double orhi = 0d;
					double orlo = 0d;
					double pval = 1d;
					double log10pval = 0d;
					double bf = 0d;
					double posterior = 0d;
					Feature assocregion = null;

					if (chrcol != -1) {
						chr = Chromosome.parseChr(elems[chrcol]);
					}
					if (poscol != -1) {
						pos = Integer.parseInt(elems[poscol]);
					}
					if (idcol != -1) {
						id = elems[idcol];
					}
					if (ncol != -1) {
						n = Integer.parseInt(elems[ncol]);
					}
					if (mafcol != -1) {
						maf = Double.parseDouble(elems[mafcol]);
					}
					if (deviancenullcol != -1) {
						deviancenull = Double.parseDouble(elems[deviancenullcol]);
					}
					if (deviancegenocol != -1) {
						deviancegeno = Double.parseDouble(elems[deviancegenocol]);
					}
					if (dfcol != -1) {
						df = Integer.parseInt(elems[dfcol]);
					}

					if (betacol != -1) {
						beta = Double.parseDouble(elems[betacol]);
					}
					if (secol != -1) {
						se = Double.parseDouble(elems[secol]);
					}
					if (orcol != -1) {
						or = Double.parseDouble(elems[orcol]);
					}
					if (orhicol != -1) {
						orhi = Double.parseDouble(elems[orhicol]);
					}
					if (orlocol != -1) {
						orlo = Double.parseDouble(elems[orlocol]);
					}

					if (pvalcol != -1) {
						pval = Double.parseDouble(elems[pvalcol]);
					}

					if (log10pvalcol != -1) {
						log10pval = Double.parseDouble(elems[log10pvalcol]);
					}

					if (regioncol != -1) {
						String regionStr = elems[regioncol];
						assocregion = Feature.parseFeature(regionStr);
					}

					if (bfcol != -1) {
						bf = Double.parseDouble(elems[bfcol]);
					}
					if (posteriorcol != -1) {
						posterior = Double.parseDouble(elems[posteriorcol]);
					}


					Feature snp = new Feature(chr, pos, pos);
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
						result.setOr(or);
						result.setOrHi(orhi);
						result.setOrLo(orlo);
						result.setPval(pval);
						result.setBf(bf);
						result.setPosterior(posterior);
						result.setRegion(assocregion);
					}
				}
			}

			ln = tf.readLine();
		}
		tf.close();
		return results;
	}

	public String getHeader(){
		String str = "VariantID" +
				"\tN" +
				"\tMAF" +
				"\tDevianceNull" +
				"\tDfNull" +
				"\tDevianceGeno" +
				"\tDfAlt" +
				"\tBeta(Genotype)" +
				"\tSE(Genotype)" +
				"\tOR" +
				"\tOR-Hi" +
				"\tOR-Lo" +
				"\tPval" +
				"\t-Log10(pval)";
		return str;
	}
}
