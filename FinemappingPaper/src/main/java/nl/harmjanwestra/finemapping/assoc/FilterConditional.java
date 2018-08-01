package nl.harmjanwestra.finemapping.assoc;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 8/23/16.
 */
public class FilterConditional {

	public static void main(String[] args) {
		int nriters = 5;
		double threshold = 7.5E-7;
		FilterConditional q = new FilterConditional();
		try {
			String bonferroni = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/BonferroniThresholds/META.txt";
			String modelfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Conditional/META-assoc0.3-COSMO-gwas-conditional-topvariants.txt";
			String output = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Conditional/output/META.txt";
			q.run(threshold, bonferroni, modelfile, nriters, output);

			bonferroni = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/BonferroniThresholds/RA.txt";
			modelfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Conditional/RA-assoc0.3-COSMO-gwas-conditional-topvariants.txt";
			output = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Conditional/output/RA.txt";
			q.run(threshold, bonferroni, modelfile, nriters, output);

			bonferroni = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/BonferroniThresholds/T1D.txt";
			modelfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Conditional/T1D-assoc0.3-COSMO-gwas-conditional-topvariants.txt";
			output = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Conditional/output/T1D.txt";
			q.run(threshold, bonferroni, modelfile, nriters, output);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(double threshold, String bonferronifile, String modelfile, int nriters, String output) throws IOException {

		threshold = -Math.log10(threshold);

		HashMap<String, Double> bonferronithresholds = loadBonferroni(bonferronifile);
		HashMap<String, ArrayList<String>> variantsPerRegion = new HashMap<String, ArrayList<String>>();
		HashSet<String> allRegions = new HashSet<String>();
		for (int iter = 0; iter < nriters; iter++) {
			TextFile tf = new TextFile(modelfile, TextFile.R);
			tf.readLine();
			String[] elems = tf.readLineElems(TextFile.tab);

			while (elems != null) {
				String region = elems[0];

				Integer i = Integer.parseInt(elems[1]);
				if (i.equals(iter)) {
					Double p = Double.parseDouble(elems[3]);
					String variant = elems[2];
					if (i.equals(0)) {
						if (p > threshold) {
							allRegions.add(region);
							ArrayList<String> variants = new ArrayList<>();
							variants.add(variant);
							variantsPerRegion.put(region, variants);
						}
					} else {
						Double bonferroni = bonferronithresholds.get(region);
						bonferroni = -Math.log10(bonferroni);
						if (p > bonferroni) {
							ArrayList<String> variants = variantsPerRegion.get(region);
							if (variants != null) {
								variants.add(variant);
								variantsPerRegion.put(region, variants);
							}
						}
					}
				}

				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
		}

		System.out.println(variantsPerRegion.size() + " significant loci");
		int ctr = 0;

		TextFile out = new TextFile(output, TextFile.W);
		TextFile bedout = new TextFile(output + ".bed", TextFile.W);
		out.writeln("Region\tn\tVariants");
		for (String region : allRegions) {
			ArrayList<String> variants = variantsPerRegion.get(region);
			if (variants.size() > 1) {
				String[] regionelems = region.split("_");
				String[] poselems = regionelems[1].split("-");
				bedout.writeln(regionelems[0] + "\t" + poselems[0] + "\t" + poselems[1]);
				out.writeln(region + "\t" + variants.size() + "\t" + Strings.concat(variants, Strings.semicolon));
				ctr++;
			}
		}
		out.close();
		bedout.close();
		System.out.println(ctr + " have independent effects.");
	}

	private HashMap<String, Double> loadBonferroni(String bonferronifile) throws IOException {
		HashMap<String, Double> output = new HashMap<>();
		TextFile tf = new TextFile(bonferronifile, TextFile.R);

		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String region = elems[0];
			Double threshold = Double.parseDouble(elems[2]);
			output.put(region, threshold);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		return output;
	}
}
