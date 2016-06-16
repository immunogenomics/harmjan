package nl.harmjanwestra.finemapping;

import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by Harm-Jan on 06/16/16.
 */
public class DeterminePositionWithinExhaustiveResults {

	public void run() throws IOException {

		String regionfile = "";
		String assocfileprefix = "";
		String modelfileprefix = "";
		String outputfile = "";

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(regionfile);


		TextFile outf = new TextFile(outputfile, TextFile.W);
		outf.writeln("Region\tVariant1\tVariant2\tP\tNrLowerPvals\tNrTotalPvals\tTopVariant1\tTopVariant2\tTopP");
		for (int r = 0; r < regions.size(); r++) {
			Feature region = regions.get(r);
			String modelfile = modelfileprefix + region.getChromosome().getNumber() + "-gwas-conditional-models.txt";
			String assocfile = assocfileprefix + region.getChromosome().getNumber() + "-pairwise.txt.gz";
			String[] variants = loadConditionalVariants(region, modelfile);

			if (variants == null) {
				System.err.println("Model file is broken for some reason: " + modelfile);
				System.exit(-1);
			}

			// get all pvalues in region
			TextFile tf2 = new TextFile(assocfile, TextFile.R);
			String ln = tf2.readLine();
			String[] elems = tf2.readLineElems(TextFile.tab);

			Double pval = 0d;
			ArrayList<Double> allPvals = new ArrayList<>();

			String regionStr = region.toString();

			String[] topEffect = null;
			Double maxP = 0d;

			while (elems != null) {

				// #Chromosome     Pos1    Id1     Pos2    Id2
				Feature s1 = new Feature();
				Chromosome chr = Chromosome.parseChr(elems[0]);
				Integer pos1 = Integer.parseInt(elems[1]);
				String snp1Id = elems[2];
				Integer pos2 = Integer.parseInt(elems[3]);
				String snp2Id = elems[4];

				SNPFeature snp1 = new SNPFeature();
				snp1.setChromosome(chr);
				snp1.setStart(pos1);
				snp1.setStop(pos1);
				snp1.setName(snp1Id);

				SNPFeature snp2 = new SNPFeature();
				snp2.setChromosome(chr);
				snp2.setStart(pos2);
				snp2.setStop(pos2);
				snp2.setName(snp2Id);


				if (region.overlaps(snp1) && region.overlaps(snp2)) {

					String snp1str = snp1.toString();
					String snp2str = snp2.toString();
					Double p = Double.parseDouble(elems[elems.length - 1]);
					if ((snp1str.equals(variants[0]) && snp2str.equals(variants[1]))
							|| (snp1str.equals(variants[1]) && snp2str.equals(variants[0]))) {
						pval = p;
					}
					if (p > maxP) {
						topEffect = elems;
						maxP = p;
					}
					allPvals.add(p);
				}
				elems = tf2.readLineElems(TextFile.tab);
			}
			tf2.close();

			int nrLowerPvals = 0;
			for (Double d : allPvals) {
				if (d > pval) {
					nrLowerPvals++;
				}
			}


			Chromosome chr = Chromosome.parseChr(topEffect[0]);
			Integer pos1 = Integer.parseInt(topEffect[1]);
			String snp1Id = topEffect[2];
			Integer pos2 = Integer.parseInt(topEffect[3]);
			String snp2Id = topEffect[4];
			SNPFeature snp1 = new SNPFeature();
			snp1.setChromosome(chr);
			snp1.setStart(pos1);
			snp1.setStop(pos1);
			snp1.setName(snp1Id);

			SNPFeature snp2 = new SNPFeature();
			snp2.setChromosome(chr);
			snp2.setStart(pos2);
			snp2.setStop(pos2);
			snp2.setName(snp2Id);

			String output = regionStr + "\t" + variants[0] + "\t" + variants[1] + "\t" + pval + "\t" + nrLowerPvals + "\t" + allPvals.size() + "\t" + snp1.toString() + "\t" + snp2.toString() + "\t" + maxP;
			System.out.println(region.toString() + "Conditional effect: " + (nrLowerPvals) + " lower pvals out of " + allPvals.size());

			outf.writeln(output);
		}

		outf.close();

	}

	private String[] loadConditionalVariants(Feature region, String modelfile) throws IOException {

		TextFile tf = new TextFile(modelfile, TextFile.R);
		String ln = tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String modelregion = elems[0];
			if (region.toString().equals(modelregion)) {
				String iter = elems[1];
				if (iter.equals("2")) {
					String snps = elems[2];
					return snps.split(";");
				}

			}
			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();

		return null;
	}

}
