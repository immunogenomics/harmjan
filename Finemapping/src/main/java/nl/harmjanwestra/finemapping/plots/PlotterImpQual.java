package nl.harmjanwestra.finemapping.plots;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.panels.BoxPlotPanel;
import nl.harmjanwestra.utilities.graphics.panels.HistogramPanel;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

/**
 * Created by hwestra on 5/20/16.
 */
public class PlotterImpQual {

	int width = 640;
	int height = 480;
	boolean onlyIc = false;
	boolean windows = false;

	public static void main(String[] args) {

		PlotterImpQual piq = new PlotterImpQual();
		try {
			piq.run();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}

	}

	public void run() throws IOException, DocumentException {

		String[] files = new String[]{
				"/Data/tmp/2016-05-23/T1D-EUR-stats.vcf.gz",
				"/Data/tmp/2016-05-23/T1D-COSMO-stats.vcf.gz",
				"/Data/tmp/2016-05-23/T1D-HRC-COSMO.vcf.gz",
				"/Data/tmp/2016-05-23/T1D-HRC-COSMO-w100kb.vcf.gz",
				"/Data/tmp/2016-05-23/T1D-HRC-HRC-w100kb.vcf.gz"
		};

		String[] labels = new String[]{"EUR", "COSMO", "HRC-COSMO", "HRC-COSMO-w100kb", "HRC-HRC-w100kb"};
		String variantsOnIC = "/Data/tmp/2016-05-20/T1D-recode-stats.vcf.gz";

		String[] files2 = new String[]{
				"/Data/tmp/2016-05-20/T1D/ImmunoChipGenotyped.txt"
		};
		String bedregions = "/Data/tmp/2016-05-23/AllICLoci.bed";
		String outdir = "/Data/tmp/2016-05-23/T1D-plotsImpQual/";

		if (windows) {

			files = new String[]{
					"D:\\tmp\\2016-05-19\\T1D-EUR-stats.vcf.gz",
					"D:\\tmp\\2016-05-19\\T1D-COSMO-stats.vcf.gz",
					"D:\\tmp\\2016-05-19\\T1D-HRC-COSMO.vcf.gz",
			};
			variantsOnIC = "D:\\tmp\\2016-05-19\\T1D-recode-stats.vcf.gz";
			bedregions = "D:\\tmp\\2016-05-19\\AllICLoci.bed";
			outdir = "D:\\tmp\\2016-05-19\\T1D-plotsImpQual\\";
		}

		if (!Gpio.exists(outdir)) {
			Gpio.createDir(outdir);
		}

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> bedfileRegions = reader.readAsList(bedregions);

		HashSet<String> variantHash = null;
		if (variantsOnIC != null) {
			variantHash = loadVariantHash(variantsOnIC);
		}

		boolean includeindels = true;
		boolean usemafthreshold = false;
		boolean requireabovemaf = false;
		boolean plotOnlyImputed = false;
		double mafthreshold = 0.005;

		String out = "";
		includeindels = true;
		usemafthreshold = true;
		requireabovemaf = true;
		plotOnlyImputed = true;

		out = outdir + "plot1-imputed-impqual-withIndels-mafgt0.005.pdf";
		//plot1(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
		out = outdir + "plot1-imputed-impqual-withoutIndels-mafgt0.005.pdf";
		includeindels = false;
		//plot1(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);

		usemafthreshold = true;
		requireabovemaf = false;
		out = outdir + "plot1-imputed-impqual-withoutIndels-maflt0.005.png";
		plot1(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);


//		out = outdir + "plot1-impqual-unfiltered.png";
//		plot1(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);

//		out = outdir + "plot2-impqual-unfiltered.png";
//		plot2(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);

//		out = outdir +"plot2-immunochip.png";
//		onlyIc = true;
//		plot2(files2, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		onlyIc = false;

//		includeindels = false;
//		usemafthreshold = false;
//		requireabovemaf = false;
//		plotOnlyImputed = false;
//
//		out = outdir + "plot1-impqual-noindels.png";
//		plot1(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-impqual-noindels.pdf";
//		plot2(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//
//		includeindels = false;
//		usemafthreshold = true;
//		requireabovemaf = true;
//		plotOnlyImputed = false;
//		out = outdir + "plot1-impqual-noindels-mafgt" + mafthreshold + ".png";
//		plot1(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot1-impqual-withindels-mafgt" + mafthreshold + ".png";
//		plot1(files, labels, out, true, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-impqual-noindels-mafgt0.01.png";
//		plot2(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = false;
//		usemafthreshold = true;
//		requireabovemaf = false;
//		plotOnlyImputed = false;
//
//		out = outdir + "plot1-impqual-noindels-maflt" + mafthreshold + ".png";
//		plot1(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-impqual-noindels-maflt0.01.png";
//		plot2(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		//
//		includeindels = true;
//		usemafthreshold = true;
//		requireabovemaf = true;
//		plotOnlyImputed = false;
//
//		out = outdir + "plot1-impqual-withindels-mafgt" + mafthreshold + ".png";
//		plot1(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = true;
//		usemafthreshold = true;
//		requireabovemaf = false;
//		plotOnlyImputed = false;
//
//		out = outdir + "plot1-impqual-withindels-maflt" + mafthreshold + ".png";
//		plot1(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//
//		//
//
//		includeindels = true;
//		usemafthreshold = false;
//		requireabovemaf = false;
//		plotOnlyImputed = true;
//
//		out = outdir + "plot1-impqual-imputedonly.png";
//		plot1(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-impqual-imputedonly.png";
//		plot2(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = false;
//		usemafthreshold = false;
//		requireabovemaf = false;
//		plotOnlyImputed = true;
//
//		out = outdir + "plot1-impqual-imputedonly-noindels.png";
//		plot1(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-impqual-imputedonly-noindels.pdf";
//		plot2(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = false;
//		usemafthreshold = true;
//		requireabovemaf = true;
//		plotOnlyImputed = true;
//
//		out = outdir + "plot1-impqual-imputedonly-noindels-mafgt" + mafthreshold + ".png";
//		plot1(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-impqual-imputedonly-noindels-mafgt0.01.png";
//		plot2(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = false;
//		usemafthreshold = true;
//		requireabovemaf = false;
//		plotOnlyImputed = true;
//
//		out = outdir + "plot1-impqual-imputedonly-noindels-maflt" + mafthreshold + ".png";
//		plot1(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-impqual-imputedonly-noindels-maflt0.01.png";
//		plot2(files, labels, out, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);

	}

	// #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
	public boolean isIndel(String[] elems) {

		String alleles = elems[3] + "," + elems[4];
		String[] alleleElems = alleles.split(",");

		for (String s : alleleElems) {
			if (s.length() > 1) {
				return true;
			}
		}
		return false;
	}

	private HashSet<String> loadVariantHash(String variantsOnIC) throws IOException {
		TextFile tf = new TextFile(variantsOnIC, TextFile.R);
		HashSet<String> variantIds = new HashSet<String>();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (!elems[0].startsWith("#")) {
				String id = elems[0] + "_" + elems[1] + "_" + elems[2];
				variantIds.add(id);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		return variantIds;
	}

	private double getMaf(String[] elems) {
		String[] infoElems = elems[7].split(";");
		double maf = 1;
		for (String info : infoElems) {
			if (info.startsWith("AF")) {
				String[] elems2 = info.split("=");
				String[] elems3 = elems2[1].split(",");
				for (String afStr : elems3) {
					double af = Double.parseDouble(afStr);
					if (af > 0.5) {
						af = 1 - af;
					}
					if (af < maf) {
						maf = af;
					}

				}
			}
		}
		return maf;
	}

	private double getInfo(String[] elems) {
		String[] infoElems = elems[7].split(";");
		double infoscore = 0;
		for (String info : infoElems) {
			if (info.startsWith("INFO")) {
				String[] elems2 = info.split("=");
				infoscore = Double.parseDouble(elems2[1]);
			}
		}
		return infoscore;
	}


	public boolean mafbelowfthreshold(String[] elems, double t) {
		double maf = getMaf(elems);
		if (maf < t) {
			return true;
		}
		return false;
	}

	private boolean isVariantOnIC(String[] elems, HashSet<String> variantHash) {
		String variant = elems[0] + "_" + elems[1] + "_" + elems[2];
		return variantHash.contains(variant);
	}

	private boolean isWithinRegion(ArrayList<Feature> list, String[] elems) {
		Feature varfeat = new Feature();

		int pos = Integer.parseInt(elems[1]);
		varfeat.setChromosome(Chromosome.parseChr(elems[0]));
		varfeat.setStart(pos);
		varfeat.setStop(pos);
		for (Feature f : list) {
			if (f.overlaps(varfeat)) {
				return true;
			}
		}
		return false;
	}


	public void plot1(String[] files, String[] labels, String out,
					  boolean includeIndels,
					  boolean usemafthreshold,
					  boolean requireabovemaf,
					  double mafthreshold,
					  boolean plotOnlyImputed,
					  HashSet<String> variantHash,
					  ArrayList<Feature> bedregions
	) throws IOException, DocumentException {
		// plot 1: x-axis nr of variants, y-axis correlation,
		ArrayList<ArrayList<Double>> vals = new ArrayList<ArrayList<Double>>();
		int maxSize = 0;
		for (int i = 0; i < files.length; i++) {
			String file = files[i];
			ArrayList<Double> corvals = new ArrayList<>();
			TextFile tf = new TextFile(file, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				if (elems.length >= 8 && !elems[0].startsWith("#")) {
					boolean isIndel = isIndel(elems);
					boolean belowmaf = mafbelowfthreshold(elems, mafthreshold);

					boolean include = isWithinRegion(bedregions, elems);
					if (!includeIndels && isIndel) {
						include = false;
					}
					if (usemafthreshold && requireabovemaf && belowmaf) {
						include = false;
					} else if (usemafthreshold && !requireabovemaf && !belowmaf) {
						include = false;
					}


					if (variantHash != null && plotOnlyImputed) {
						boolean isOnIc = isVariantOnIC(elems, variantHash);
						if (plotOnlyImputed && isOnIc) {
							include = false;
						} else if (!plotOnlyImputed && !isOnIc) {
							include = false;
						}
					}


					if (include) {
						double info = getInfo(elems);
						corvals.add(info);
					}
				}

				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			Collections.sort(corvals, Collections.reverseOrder());
			System.out.println(corvals.size() + " vals in file " + file);
			if (corvals.size() > maxSize) {
				maxSize = corvals.size();
			}
			vals.add(corvals);
		}
		System.out.println(maxSize + " total vals");

		double[][] x = new double[files.length][maxSize];
		double[][] y = new double[files.length][maxSize];
		for (int ds = 0; ds < files.length; ds++) {

			for (int i = 0; i < maxSize; i++) {
				x[ds][i] = i;
			}
			ArrayList<Double> corvals = vals.get(ds);

			for (int i = 0; i < corvals.size(); i++) {
				y[ds][i] = corvals.get(i);
			}
			for (int i = corvals.size(); i < maxSize; i++) {
				y[ds][i] = Double.NaN;
			}
		}

		vals = null;
		Grid grid = new Grid(width, height, 1, 1, 100, 100);
		ScatterplotPanel panel = new ScatterplotPanel(1, 1);
		panel.setData(x, y);
		Range range = new Range(0, 0, maxSize, 1);
		panel.setDataRange(range);
//		range.roundX();

		panel.setDatasetLabels(labels);
		grid.addPanel(panel);
		grid.draw(out);
	}

	public void plot2(String[] files, String[] datasetLabels, String out,
					  boolean includeIndels,
					  boolean usemafthreshold,
					  boolean requireabovemaf,
					  double mafthreshold,
					  boolean plotOnlyImputed,
					  HashSet<String> variantHash,
					  ArrayList<Feature> bedregions
	) throws IOException, DocumentException {
		// plot 2: x-axis maf, y-axis correlation (boxplot)

		String[] binLabels = new String[]{
				"0 - 0.005",
				"0.005 - 0.01",
				"0.01 - 0.02",
				"0.02 - 0.03",
				"0.03 - 0.04",
				"0.04 - 0.05",
				"0.05 - 0.1",
				"0.1 - 0.2",
				"0.2 - 0.3",
				"0.3 - 0.4",
				"0.4 - 0.5"
		};

		double[] lowerthreshold = new double[]{
				0,
				0.005,
				0.01,
				0.02,
				0.03,
				0.04,
				0.05,
				0.1,
				0.2,
				0.3,
				0.4
		};
		double[] upperthreshold = new double[]{
				0.005,
				0.01,
				0.02,
				0.03,
				0.04,
				0.05,
				0.1,
				0.2,
				0.3,
				0.4,
				0.5
		};

		double[][][] bins = new double[files.length][upperthreshold.length][];

		double[][] binfreqs = new double[files.length][upperthreshold.length];

		for (int i = 0; i < files.length; i++) {
			int loaded = 0;
			String file = files[i];

			TextFile tf = new TextFile(file, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			ArrayList<ArrayList<Double>> filebins = new ArrayList<>();
			for (int bin = 0; bin < upperthreshold.length; bin++) {
				filebins.add(new ArrayList<Double>());
			}
			while (elems != null) {
				if (elems.length >= 8 && !elems[0].startsWith("#")) {
					double val = getInfo(elems);
					double maf = getMaf(elems);
					boolean isIndel = isIndel(elems);
					boolean include = isWithinRegion(bedregions, elems);
					if (!includeIndels && isIndel) {
						include = false;
					}

					if (variantHash != null && plotOnlyImputed) {
						boolean isOnIc = false;
						isOnIc = isVariantOnIC(elems, variantHash);
						if (plotOnlyImputed && isOnIc) {
							include = false;
						} else if (!plotOnlyImputed && !isOnIc) {
							include = false;
						}
					}

					if (onlyIc && variantHash != null) {
						boolean isOnIc = false;
						isOnIc = isVariantOnIC(elems, variantHash);
						include = isOnIc;
					}


					if (include) {
						for (int bin = 0; bin < upperthreshold.length; bin++) {
							if (maf >= lowerthreshold[bin] && maf <= upperthreshold[bin]) {
								filebins.get(bin).add(val);
								loaded++;
							}
						}
					} else {
//					System.out.println("Value excluded?");
					}
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();

			// convert to double[]
			for (int bin = 0; bin < upperthreshold.length; bin++) {
				ArrayList<Double> d = filebins.get(bin);
				binfreqs[i][bin] = d.size();
				double[] arr = Primitives.toPrimitiveArr(d.toArray(new Double[0]));
				bins[i][bin] = arr;
			}
			System.out.println("Loaded: " + loaded + " from file: " + file);
		}


		Grid grid = new Grid(width * 2, height, 1, 2, 100, 100);
		BoxPlotPanel panel = new BoxPlotPanel(1, 1);
		panel.setData(bins);
		panel.setDrawDataPoints(false);
		panel.useTukeysDefault(true);
		panel.setBinLabels(binLabels);
		grid.addPanel(panel);

		HistogramPanel panel2 = new HistogramPanel(1, 1);
		panel2.setData(binfreqs);
		panel2.setDatasetLabels(datasetLabels);
		panel2.setBinLabels(binLabels);

		grid.addPanel(panel2);
		grid.draw(out);


	}


}
