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

	public ArrayList<AssociationResult> loadConditionalAssocData(String file, Feature region) throws IOException {
		TextFile tf = new TextFile(file, TextFile.R);
		String ln = tf.readLine();

		model = null;

		ArrayList<AssociationResult> results = new ArrayList<AssociationResult>();
		while (ln != null) {
			if (ln.startsWith("#")) {
				// get the model
				model = ln;
			} else if (ln.startsWith("VariantID")) {
// skip header
			} else {
				String[] elems = Strings.tab.split(ln);
				if (elems.length > 4) {
					// VariantID	N	MAF	DevianceNull	DfNull	DevianceGeno	DfAlt	Beta(Genotype)	SE(Genotype)	OR	OR-Hi	OR-Lo	Pval	-Log10(pval)

					String snp = elems[0];
					String mafStr = elems[2];
					String betaStr = elems[7];
					String seStr = elems[8];
					String orStr = elems[9];
					String pvalStr = elems[elems.length - 1];

					Double maf = Double.parseDouble(mafStr);

					String[] snpelems = snp.split("-");
					String[] posElems = snpelems[0].split(":");
					Integer pos = Integer.parseInt(posElems[1]);

					Chromosome chr = Chromosome.parseChr(posElems[0]);
					Feature f = new Feature(chr, pos, pos);
					if (region.overlaps(f)) {
						Double pval = Double.parseDouble(pvalStr);

						if (betaStr.equals("null")) {
							AssociationResult r = new AssociationResult(Chromosome.parseChr(posElems[0]), pos, maf, Double.NaN, Double.NaN, Double.NaN, pval);
							if (snpelems.length > 1) {
								r.getSnp().setName(snpelems[1]);
							}
							results.add(r);
						} else {
							Double beta = Double.parseDouble(betaStr);
							if (pval == 0) {
								AssociationResult r = new AssociationResult(Chromosome.parseChr(posElems[0]), pos, maf, Double.NaN, Double.NaN, Double.NaN, pval);
								if (snpelems.length > 1) {
									r.getSnp().setName(snpelems[1]);
								}
								results.add(r);
							} else {
								Double se = Double.parseDouble(seStr);
								Double or = Double.parseDouble(orStr);
								AssociationResult r = new AssociationResult(Chromosome.parseChr(posElems[0]), pos, maf, or, beta, se, pval);
								if (snpelems.length > 1) {
									r.getSnp().setName(snpelems[1]);
								}
								results.add(r);
							}

						}
					}
				}
			}

			ln = tf.readLine();
		}
		tf.close();
		return results;
	}
}
