package nl.harmjanwestra.finemapping.plots;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.finemapping.VariantCounter;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.panels.BoxPlotPanel;
import nl.harmjanwestra.utilities.graphics.panels.HistogramPanel;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

/**
 * Created by hwestra on 5/18/16.
 */
public class PlotterAccuracy extends VariantCounter {

	// id	allele1	allele2	maf	callrate	allele21	allele22	maf 	cr	df	nrsamples	r	rsq	beta	se	impqual0	impqual1
	protected int width = 640;
	protected int height = 480;
	boolean onlyIc = false;

	public static void main(String[] args) {
		PlotterAccuracy p = new PlotterAccuracy();
		try {
			p.plotCorr();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
	}


	public void plotCorr() throws IOException, DocumentException {
		String[] files = new String[]{
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-07-10-Accuracy/T1D-EUR.txt",
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-07-10-Accuracy/T1D-COSMO.txt",
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-07-10-Accuracy/T1D-HRC-EAGLE.txt",
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-07-10-Accuracy/T1D-HRC-SHAPEIT.txt",
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-07-10-Accuracy/T1D-HRC-EAGLE-Michigan.txt"
		};

		String[] labels = new String[]{
				"EUR",
				"COSMO",
				"HRC / HRC / EAGLE",
				"HRC / HRC / SHAPEIT",
				"HRC / HRC / EAGLE / MICHIGAN"
		};
		String diseaseprefix = "T1D";

		files = new String[]{
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-07-10-Accuracy/RA-EUR.txt",
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-07-10-Accuracy/RA-COSMO.txt",
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-07-10-Accuracy/RA-HRC-w100kb.txt",
		};
		labels = new String[]{
				"EUR",
				"COSMO",
				"HRC / COSMO / EAGLE"
		};
		diseaseprefix = "RA";

		String seqpanelvcf = "/Data/tmp/2016-05-28/seqpanelfiltered-maf0005-cr0950-rd10-gq30-runNamesFixed-RASampleNamesFixed-badSamplesRemoved-mixupsFixed.vcf.gz";
		String variantsOnIC = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/RAAndT1D-recode-maf0005-ICRegionsW100kb-samplenamefix.vcf.gz-updatedRSId-stats.vcf.gz";

//		String[] files2 = new String[]{
//				"/Data/tmp/2016-05-20/T1D/ImmunoChipGenotyped.txt"
//		};
		String bedregions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
//		bedregions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci.bed";
		bedregions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		String outdir = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-07-10-Accuracy/" + diseaseprefix + "-plots/";
		double mafthreshold = 0.01;

		boolean includeId = true;
		boolean windows = false;
		String ext = "pdf";

		if (windows) {

			files = new String[]{
					"D:\\tmp\\2016-05-19\\T1D\\T1D-EUR-merged.txt",
					"D:\\tmp\\2016-05-19\\T1D\\T1D-COSMO-merged.txt",
					"D:\\tmp\\2016-05-19\\T1D\\T1D-HRC-COSMO-merged.txt",
					"D:\\tmp\\2016-05-19\\T1D\\T1D-HRC-HRC-w100kb-merged.txt",
					"D:\\tmp\\2016-05-19\\T1D\\T1D-HRC-COSMO-w100kb-merged.txt"
			};
			labels = new String[]{"EUR", "COSMO", "HRC-COSMO", "HRC-HRC-w100kb", "HRC-COSMO-w100kb"};
			variantsOnIC = "D:\\tmp\\2016-05-19\\T1D-recode-stats.vcf.gz";
			bedregions = "D:\\tmp\\2016-05-19\\AllICLoci.bed";
			outdir = "D:\\tmp\\2016-05-19\\T1D-plotsAccuracy\\";
//			files2 = new String[]{
//					"D:\\tmp\\2016-05-19\\ImmunoChipGenotyped.txt"
//			};
		}


		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> bedfileRegions = reader.readAsList(bedregions);

		if (!Gpio.exists(outdir)) {
			Gpio.createDir(outdir);
		}


		String out = "";

		HashSet<String> variantHash = null;
		if (variantsOnIC != null) {
			variantHash = loadVariantHash(variantsOnIC, includeId);
		}

		boolean includeindels = true;
		boolean usemafthreshold = false;
		boolean includeICVariants = true;


		// including indels
		int maxNrVariants = 1862;
		includeindels = true;
		usemafthreshold = false;
		includeICVariants = true;
		out = outdir + diseaseprefix + "-plot1-allvariants." + ext;
		plot1(files, labels, out, bedregions, includeindels, usemafthreshold, mafthreshold, includeICVariants, maxNrVariants, includeId, variantsOnIC, seqpanelvcf);

		// maf > 1%
		usemafthreshold = true;
		out = outdir + diseaseprefix + "-plot1-allvariants-mafgt" + mafthreshold + "." + ext;
		plot1(files, labels, out, bedregions, includeindels, usemafthreshold, mafthreshold, includeICVariants, maxNrVariants, includeId, variantsOnIC, seqpanelvcf);

		// maf > 1%, imputed variants
		usemafthreshold = true;
		includeICVariants = false;
		maxNrVariants = 693;
		out = outdir + diseaseprefix + "-plot1-allvariants-imputed-mafgt" + mafthreshold + "." + ext;
		plot1(files, labels, out, bedregions, includeindels, usemafthreshold, mafthreshold, includeICVariants, maxNrVariants, includeId, variantsOnIC, seqpanelvcf);


		// excluding indels
		maxNrVariants = 1717;
		includeindels = false;
		usemafthreshold = false;
		includeICVariants = true;
		out = outdir + diseaseprefix + "-plot1-excludedindels." + ext;
		plot1(files, labels, out, bedregions, includeindels, usemafthreshold, mafthreshold, includeICVariants, maxNrVariants, includeId, variantsOnIC, seqpanelvcf);

		// maf > 1%
		usemafthreshold = true;
		out = outdir + diseaseprefix + "-plot1-excludedindels-mafgt" + mafthreshold + "." + ext;
		plot1(files, labels, out, bedregions, includeindels, usemafthreshold, mafthreshold, includeICVariants, maxNrVariants, includeId, variantsOnIC, seqpanelvcf);

		// maf > 1%, imputed variants
		usemafthreshold = true;
		includeICVariants = false;
		maxNrVariants = 548;
		out = outdir + diseaseprefix + "-plot1-excludedindels-imputed-mafgt" + mafthreshold + "." + ext;
		plot1(files, labels, out, bedregions, includeindels, usemafthreshold, mafthreshold, includeICVariants, maxNrVariants, includeId, variantsOnIC, seqpanelvcf);


	}


	public void plot1(String[] files, String[] labels, String out, String bedregionsfile,
					  boolean includeIndels,
					  boolean usemafthreshold,
					  double mafthreshold,
					  boolean includeICVariants,
					  int maxSize,
					  boolean includeId,
					  String variantsOnIC,
					  String seqpanelvcf


	) throws IOException, DocumentException {
		// plot 1: x-axis nr of variants, y-axis correlation,
		ArrayList<ArrayList<Double>> vals = new ArrayList<ArrayList<Double>>();

		double upperthreshold = 1d;
		HashSet<String> variantsOnICHash = loadVariantHash(variantsOnIC, includeId);
		System.out.println(variantsOnICHash.size() + " total on IC");
		Triple<ArrayList<VCFVariant>, ArrayList<VCFVariant>, ArrayList<VCFVariant>> seqpanelvariants = loadSequencedVariants(
				seqpanelvcf, bedregionsfile, mafthreshold, upperthreshold, variantsOnICHash, includeId, includeIndels
		);

		ArrayList<VCFVariant> seqpanel = seqpanelvariants.getLeft();
		ArrayList<VCFVariant> variantsOnImmunoChip = seqpanelvariants.getMiddle();
		ArrayList<VCFVariant> variantsNotOnImmunoChip = seqpanelvariants.getRight();

		System.out.println(seqpanel.size() + " variants in VCF");
		System.out.println(variantsNotOnImmunoChip.size() + " not on IC");
		System.out.println(variantsOnImmunoChip.size() + " on IC");


		HashSet<String> sequencedVariantsHash = null;
		if (includeICVariants) {
			ArrayList<VCFVariant> allvars = new ArrayList<>();
			allvars.addAll(variantsNotOnImmunoChip);
			allvars.addAll(variantsOnImmunoChip);
			sequencedVariantsHash = hashit(allvars, includeId);
		} else {
			sequencedVariantsHash = hashit(variantsNotOnImmunoChip, includeId);
		}

		String[] newlabels = new String[labels.length];
		for (int f = 0; f < files.length; f++) {
			// get the imputation accuracies for these variants
			TextFile tf2 = new TextFile(files[f], TextFile.R);
			ArrayList<Double> corvals = new ArrayList<>();
			String[] elems = tf2.readLineElems(TextFile.tab);
			while (elems != null) {
				if (!elems[rsqlcol].equals("null")) {
					double val = Double.parseDouble(elems[rsqlcol]);
					double maf = Double.parseDouble(elems[maf2col]);

					String[] varElems = elems[0].split("_");

					boolean sequenced = isVariantOnIC(varElems, sequencedVariantsHash, includeId);
					if (sequenced) {
						if (!usemafthreshold || (usemafthreshold && maf > mafthreshold)) {
							corvals.add(val);
						}
					}
				}

				elems = tf2.readLineElems(TextFile.tab);
			}
			tf2.close();
			Collections.sort(corvals, Collections.reverseOrder());
			System.out.println(corvals.size() + " vals in file " + files[f]);
			newlabels[f] = labels[f] + " - " + corvals.size() + " / " + maxSize;
			vals.add(corvals);
		}

		System.out.println(maxSize + " total vals");

		double[][] x = new double[files.length][maxSize];
		double[][] y = new double[files.length][maxSize];

		for (int ds = 0; ds < files.length; ds++) {

			for (int i = 0; i < maxSize; i++) {
				x[ds][i] = ((double) i / maxSize) * 100;
			}
			ArrayList<Double> corvals = vals.get(ds);

			for (int i = 0; i < corvals.size(); i++) {
				y[ds][i] = corvals.get(i);
			}
			for (int i = corvals.size(); i < maxSize; i++) {
				y[ds][i] = 0;
			}
		}

		vals = null;
		Grid grid = new Grid(width, height, 1, 1, 100, 100);
		ScatterplotPanel panel = new ScatterplotPanel(1, 1);
		panel.setData(y, x);
		Range range = new Range(0, 0, 1, 100);
		range.roundX();
		range.roundY();
		panel.setDataRange(range);
		// panel.setDatasetLabels(newlabels);
		grid.addPanel(panel);
		grid.draw(out);
		System.out.println("Out: " + out);
	}


	public void plot2(String[] files, String[] datasetLabels, String out, int col,
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
				double val = Double.parseDouble(elems[col]);
				double maf = Double.parseDouble(elems[maf2col]);
				boolean isIndel = isIndel2(elems);
				boolean include = isWithinRegion(bedregions, elems);
				if (!includeIndels && isIndel) {
					include = false;
				}

				if (variantHash != null && plotOnlyImputed) {
					boolean isOnIc = false;
					isOnIc = isVariantOnIC(elems, variantHash, true);
					if (plotOnlyImputed && isOnIc) {
						include = false;
					} else if (!plotOnlyImputed && !isOnIc) {
						include = false;
					}
				}

				if (onlyIc && variantHash != null) {
					boolean isOnIc = false;
					isOnIc = isVariantOnIC(elems, variantHash, true);
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

	public void plot3(String[] files, String[] datasetLabels, String out, int col,
					  boolean includeIndels,
					  boolean usemafthreshold,
					  boolean requireabovemaf,
					  double mafthreshold,
					  boolean plotOnlyImputed,
					  HashSet<String> variantHash,
					  ArrayList<Feature> bedregions) throws IOException, DocumentException {
		//			// plot 3: maf vs maf
		Grid grid = new Grid(width, height, 1, files.length, 100, 100);
		for (int i = 0; i < files.length; i++) {
			String file = files[i];

			ArrayList<Double> fy = new ArrayList<>();
			ArrayList<Double> fx = new ArrayList<>();
			int ctr = 0;
			TextFile tf = new TextFile(file, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				double maf1 = Double.parseDouble(elems[maf1col]);
				double maf2 = Double.parseDouble(elems[maf2col]);
				boolean isIndel = isIndel2(elems);
				boolean belowmaf = maf2belowfthreshold(elems, mafthreshold);

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
					boolean isOnIc = isVariantOnIC(elems, variantHash, true);
					if (plotOnlyImputed && isOnIc) {
						include = false;
					} else if (!plotOnlyImputed && !isOnIc) {
						include = false;
					}
				}

				if (include) {
					fx.add(maf1);
					fy.add(maf2);
				}
				ctr++;
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();

			double[] x = Primitives.toPrimitiveArr(fx.toArray(new Double[0]));
			double[] y = Primitives.toPrimitiveArr(fy.toArray(new Double[0]));

			ScatterplotPanel panel2 = new ScatterplotPanel(1, 1);
			panel2.setData(x, y);
			grid.addPanel(panel2);
		}
		grid.draw(out);
	}

	public void plot4(String[] files, String[] datasetLabels, String out, int col,
					  boolean includeIndels,
					  boolean usemafthreshold,
					  boolean requireabovemaf,
					  double mafthreshold,
					  boolean plotOnlyImputed,
					  HashSet<String> variantHash,
					  ArrayList<Feature> bedregions) throws IOException, DocumentException {
		//			// plot 3: impqual vs beta
		Grid grid = new Grid(width, height, 1, files.length, 100, 100);
		for (int i = 0; i < files.length; i++) {
			String file = files[i];

			ArrayList<Double> fy = new ArrayList<>();
			ArrayList<Double> fx = new ArrayList<>();
			int ctr = 0;
			TextFile tf = new TextFile(file, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				double maf1 = Double.parseDouble(elems[impqual2]);
				double maf2 = Double.parseDouble(elems[betacol]);
				boolean isIndel = isIndel2(elems);
				boolean belowmaf = maf2belowfthreshold(elems, mafthreshold);

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
					boolean isOnIc = isVariantOnIC(elems, variantHash, true);
					if (plotOnlyImputed && isOnIc) {
						include = false;
					} else if (!plotOnlyImputed && !isOnIc) {
						include = false;
					}
				}

				if (include) {
					fx.add(maf1);
					fy.add(maf2);
				}
				ctr++;
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();

			double[] x = Primitives.toPrimitiveArr(fx.toArray(new Double[0]));
			double[] y = Primitives.toPrimitiveArr(fy.toArray(new Double[0]));

			ScatterplotPanel panel2 = new ScatterplotPanel(1, 1);
			panel2.setData(x, y);
			grid.addPanel(panel2);
		}
		grid.draw(out);
	}

	public void plot5(String[] files, String[] datasetLabels, String out, int col,
					  boolean includeIndels,
					  boolean usemafthreshold,
					  boolean requireabovemaf,
					  double mafthreshold,
					  boolean plotOnlyImputed,
					  HashSet<String> variantHash,
					  ArrayList<Feature> bedregions) throws IOException, DocumentException {
		//			// plot 3: impqual vs beta
		Grid grid = new Grid(width, height, 1, files.length, 100, 100);
		for (int i = 0; i < files.length; i++) {
			String file = files[i];

			ArrayList<Double> fy = new ArrayList<>();
			ArrayList<Double> fx = new ArrayList<>();
			int ctr = 0;
			TextFile tf = new TextFile(file, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				double maf1 = Double.parseDouble(elems[impqual2]);
				double maf2 = Double.parseDouble(elems[rsqlcol]);
				boolean isIndel = isIndel2(elems);
				boolean belowmaf = maf2belowfthreshold(elems, mafthreshold);

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
					boolean isOnIc = isVariantOnIC(elems, variantHash, true);
					if (plotOnlyImputed && isOnIc) {
						include = false;
					} else if (!plotOnlyImputed && !isOnIc) {
						include = false;
					}
				}

				if (include) {
					fx.add(maf1);
					fy.add(maf2);
				}
				ctr++;
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();

			double[] x = Primitives.toPrimitiveArr(fx.toArray(new Double[0]));
			double[] y = Primitives.toPrimitiveArr(fy.toArray(new Double[0]));

			ScatterplotPanel panel2 = new ScatterplotPanel(1, 1);
			panel2.setData(x, y);
			grid.addPanel(panel2);
		}
		grid.draw(out);
	}

}
