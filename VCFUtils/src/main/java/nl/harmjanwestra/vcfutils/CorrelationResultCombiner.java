package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 12/21/15.
 */
public class CorrelationResultCombiner {

	public static void main(String[] args) {
		CorrelationResultCombiner comb = new CorrelationResultCombiner();
		int nriter = 10;
		String dir = "d:\\tmp\\files\\T1D\\";
		String[] refs = new String[]{"1kg", "1kg-seq-merged", "seq", "swisscheese", "swisscheesebeagle41"};
		for (String ref : refs) {

			String[] files = new String[nriter];
			for (int i = 0; i < files.length; i++) {
				String file = dir + ref + "-iter" + (i + 1) + ".txt";
				files[i] = file;

			}
			String outfile = dir + "" + ref + ".txt";

			try {
				comb.concat(files, outfile);
			} catch (IOException e) {
				e.printStackTrace();
			}

		}

		try {
			String[] files = new String[refs.length];
			dir = "d:\\tmp\\files\\T1D\\";
			for (int i = 0; i < refs.length; i++) {
				files[i] = dir + refs[i] + ".txt";
			}
			String outprefix = "d:\\tmp\\files\\T1D";
			comb.mergeIntoOneBigTable(files, refs, outprefix);

			files = new String[refs.length];
			dir = "d:\\tmp\\files\\RA\\";
			for (int i = 0; i < refs.length; i++) {
				files[i] = dir + refs[i] + ".txt";
			}
			outprefix = "d:\\tmp\\files\\RA";
			comb.mergeIntoOneBigTable(files, refs, outprefix);
		} catch (IOException e) {

		}

	}

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
		int[] samplesizes = new int[variantList.size()];
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
								imputationquals[id] += result.rsqBeagle1;
								samplesizes[id] += result.nrSamples;
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

//		String header = "Variant\tmaf1\tmaf2\tpearsonrsquared\tbeaglersquared\tn";

		// 9-34970885-rs77176551   G       G,A     0.5     A       G,A     0.0019759791583874244   1       2237    0.04873531816737055     0.002375131236874838    0.0     null

		// outfile.writeln(header);

		for (int i = 0; i < variantList.size(); i++) {
			String variant = variantList.get(i);

			double n = ctrs[i];
			if (n > 0) {
				variant += "\t-t\t-\t" + (mafs1[i] / n)
						+ "\t-t\t-\t" + (mafs2[i] / n)
						+ "\t-"
						+ "\t" + samplesizes[i]
						+ "\t-"
						+ "\t" + (correlations[i] / n)
						+ "\t" + (imputationquals[i] / n)
						+ "\t" + n;
			} else {
				variant += "\t-t\t-\t" + 0
						+ "\t-t\t-\t" + 0
						+ "\t-"
						+ "\t" + 0
						+ "\t-"
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

	public ArrayList<CorrelationResult> getData(String file) throws IOException {
		ArrayList<CorrelationResult> output = new ArrayList<>();
		TextFile tf = new TextFile(file, TextFile.R);

		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {


			// var1.getMinorAllele() + "\t" + Strings.concat(var1.getAlleles(), Strings.comma) + "\t" + var1.getMAF() + "\t" + var1.getCallrate()
			// + "\t" + r + "\t" + rsq + "\t" + info1.get("AR2") + "\t" + info2.get("AR2")
//      9-34726524-rs1810820    null    G,A     0.0     0.0     null    G,A     0.0     0.0     1       22624   0.894894178229274       0.8008355902286476      null    0.4125
//		9-34970885-rs77176551   G       G,A     0.5     A       G,A     0.0019759791583874244   1       2237    0.04873531816737055     0.002375131236874838    0.0     null
			if (elems.length > 14) {

				CorrelationResult result = new CorrelationResult();
				result.variant = elems[0];
				try {
					result.maf1 = Double.parseDouble(elems[3]);
				} catch (NumberFormatException e) {

				}
				try {
					result.maf2 = Double.parseDouble(elems[7]);
				} catch (NumberFormatException e) {

				}
				try {
					result.nrSamples = Integer.parseInt(elems[10]);
				} catch (NumberFormatException e) {

				}
				try {
					result.rsqPearson = Double.parseDouble(elems[12]);
				} catch (NumberFormatException e) {

				}

				try {
					result.rsqBeagle1 = Double.parseDouble(elems[13]);
				} catch (NumberFormatException e) {

				}

				try {
					result.rsqBeagle2 = Double.parseDouble(elems[14]);
				} catch (NumberFormatException e) {

				}

				output.add(result);
			}

			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();

		return output;
	}

	public void concat(String[] files, String out) throws IOException {
		TextFile outf = new TextFile(out, TextFile.W);
		for (String f : files) {
			TextFile in = new TextFile(f, TextFile.R);
			String ln = in.readLine();
			while (ln != null) {
				outf.writeln(ln);
				ln = in.readLine();
			}
			in.close();
		}
		outf.close();
	}

	public void mergeIntoOneBigTable(String[] files, String[] names, String outprefix) throws IOException {


		// get list of unique variants
		HashSet<String> uniqueVariants = new HashSet<String>();
		for (String file : files) {
			TextFile tf = new TextFile(file, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				String variant = elems[0];
				uniqueVariants.add(variant);
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
		}
		System.out.println(uniqueVariants.size() + " variants found.");


		// now index
		HashMap<String, Integer> variantToId = new HashMap<String, Integer>();
		ArrayList<String> variants = new ArrayList<>();
		int ctr = 0;
		for (String variant : uniqueVariants) {
			variantToId.put(variant, ctr);
			variants.add(variant);
			ctr++;
		}

		double[][] pearsons = new double[variants.size()][files.length];
		double[][] mafs1 = new double[variants.size()][files.length];
		double[][] mafs2 = new double[variants.size()][files.length];
		double[][] beaglersq1 = new double[variants.size()][files.length];
		double[][] beaglersq2 = new double[variants.size()][files.length];

		for (int f = 0; f < files.length; f++) {
			TextFile tf = new TextFile(files[f], TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				if (elems.length > 14) {
					String variant = elems[0];
					Integer index = variantToId.get(variant);
					if (index != null) {


						double p = Double.NaN;
						double d1 = Double.NaN;
						double d2 = Double.NaN;
						double maf1 = Double.NaN;
						double maf2 = Double.NaN;

						try {
							p = Double.parseDouble(elems[12]);
						} catch (NumberFormatException e) {

						}

						try {
							maf1 = Double.parseDouble(elems[3]);
						} catch (NumberFormatException e) {

						}

						try {
							maf2 = Double.parseDouble(elems[7]);
						} catch (NumberFormatException e) {

						}

						try {
							d1 = Double.parseDouble(elems[13]);
						} catch (NumberFormatException e) {

						}

						try {
							d2 = Double.parseDouble(elems[14]);
						} catch (NumberFormatException e) {

						}

						pearsons[index][f] = p;

						mafs1[index][f] = maf1;
						mafs2[index][f] = maf2;
						beaglersq1[index][f] = d1;
						beaglersq2[index][f] = d2;

					}
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
		}

		String[] variantArr = new String[variants.size()];
		for (int i = 0; i < variants.size(); i++) {
			variantArr[i] = variants.get(i);
		}

		write(pearsons, variantArr, names, outprefix + "-pearsons.txt");
		write(mafs1, variantArr, names, outprefix + "-mafs1.txt");
		write(mafs2, variantArr, names, outprefix + "-mafs2.txt");
		write(beaglersq1, variantArr, names, outprefix + "-beaglersq1.txt");
		write(beaglersq2, variantArr, names, outprefix + "-beaglersq2.txt");


	}

	private void write(double[][] d, String[] rows, String[] cols, String outfilename) throws IOException {
		String header = "-\t" + Strings.concat(cols, Strings.tab);
		TextFile out = new TextFile(outfilename, TextFile.W);
		out.writeln(header);
		for (int r = 0; r < d.length; r++) {
			String ln = rows[r] + "\t" + Strings.concat(d[r], Strings.tab);
			out.writeln(ln);
		}
		out.close();
	}


}
