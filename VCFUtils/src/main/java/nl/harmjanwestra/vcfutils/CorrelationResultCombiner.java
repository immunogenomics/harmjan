package nl.harmjanwestra.vcfutils;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 12/21/15.
 */
public class CorrelationResultCombiner {

	public void run(String indir, String outfilename, String refname, int nriter, int nrBatches) throws IOException {

		// 1kg-chr9-batch2-perc0.05-iter2-correlations.txt.gz
		HashSet<String> uniqueVariants = new HashSet<String>();
		for (int chr = 1; chr < 23; chr++) {
			for (int batch = 0; batch < nrBatches; batch++) {
				for (int iter = 0; iter < nriter; iter++) {
					String file = indir + refname + "-Chr" + chr + "-batch" + batch + "-perc0.05-iter" + iter + "-correlations.txt.gz";
					if (Gpio.exists(file)) {
						ArrayList<String> variants = getListOfVariants(file);
						uniqueVariants.addAll(variants);
					}
				}
			}
		}

		// now all variants have been loaded.. index them
		HashMap<String, Integer> variantToId = new HashMap<String, Integer>();
		ArrayList<String> variantList = new ArrayList<String>();
		int ctr = 0;
		for (String variant : uniqueVariants) {
			variantToId.put(variant, ctr);
			variantList.add(variant);
			ctr++;
		}


		double[] correlations = new double[variantList.size()];
		double[] imputationquals = new double[variantList.size()];
		double[] ctrs = new double[variantList.size()];
		double[] mafs1 = new double[variantList.size()];
		double[] mafs2 = new double[variantList.size()];

		for (int iter = 0; iter < nriter; iter++) {
			// get the list of variants:
			for (int chr = 1; chr < 23; chr++) {
				// for each chromosome
				for (int batch = 0; batch < nrBatches; batch++) {
					String file = indir + refname + "-Chr" + chr + "-batch" + batch + "-perc0.05-iter" + iter + "-correlations.txt.gz";
					if (Gpio.exists(file)) {

						ArrayList<CorrelationResult> results = getData(file);

						for (CorrelationResult result : results) {
							Integer id = variantToId.get(result.variant);
							if (id != null) {
								correlations[id] += result.rsqPearson;
								imputationquals[id] += result.rsqBeagle;
								mafs1[id] += result.maf1;
								mafs2[id] += result.maf2;
								ctrs[id] += 1d;
							}
						}
					}
				}
			}
		}


		// String outfilename = outdir + "/" + dataset + "-" + ref + ".txt.gz";

		TextFile outfile = new TextFile(outfilename, TextFile.W);

		String header = "Variant\tmaf1\tmaf2\tpearsonrsquared\tbeaglersquared\tn";
		outfile.writeln(header);

		for (int i = 0; i < variantList.size(); i++) {
			String variant = variantList.get(i);
			double n = ctrs[i];
			if (n > 0) {
				variant += "\t" + (mafs1[i] / n)
						+ "\t" + (mafs2[i] / n)
						+ "\t" + (correlations[i] / n)
						+ "\t" + (imputationquals[i] / n)
						+ "\t" + n;
			} else {
				variant += "\t" + 0
						+ "\t" + 0
						+ "\t" + 0
						+ "\t" + 0
						+ "\t" + n;
			}
			outfile.writeln(variant);
		}
		outfile.close();


	}

	private ArrayList<String> getListOfVariants(String file) throws IOException {
		ArrayList<String> output = new ArrayList<String>();
		TextFile tf = new TextFile(file, TextFile.R);

		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String variant = elems[0];
			output.add(variant);
			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();

		return output;
	}

	private ArrayList<CorrelationResult> getData(String file) throws IOException {
		ArrayList<CorrelationResult> output = new ArrayList<>();
		TextFile tf = new TextFile(file, TextFile.R);

		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {

//		9-34970885-rs77176551   G       G,A     0.5     A       G,A     0.0019759791583874244   1       2237    0.04873531816737055     0.002375131236874838    0.0     null
			String variant = elems[0];

			double maf1 = Double.parseDouble(elems[3]);
			double maf2 = Double.parseDouble(elems[6]);

			double rsqPearson = Double.parseDouble(elems[10]);
			double rsqBeagle = Double.parseDouble(elems[11]);

			CorrelationResult result = new CorrelationResult();
			result.maf1 = maf1;
			result.maf2 = maf2;
			result.rsqBeagle = rsqBeagle;
			result.rsqPearson = rsqPearson;
			result.variant = variant;

			output.add(result);

			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();

		return output;
	}

	class CorrelationResult {
		double maf1;
		double maf2;
		double rsqPearson;
		double rsqBeagle;
		String variant;
	}

}
