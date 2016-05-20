package nl.harmjanwestra.vcfutils.plots;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.BoxPlotPanel;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

/**
 * Created by hwestra on 5/18/16.
 */
public class Plotter {


	// id	allele1	allele2	maf	callrate	allele21	allele22	maf 	cr	df	nrsamples	r	rsq	beta	se	impqual0	impqual1
	int idcol = 0;
	int minorallele1col = 1;
	int aleleles1col = 2;
	int maf1col = 3;
	int cr1col = 4;
	int minorallele2col = 5;
	int aleleles2col = 6;
	int maf2col = 7;
	int cr2col = 8;
	int dfcol = 9;
	int samplecol = 10;
	int rcol = 11;
	int rsqlcol = 12;
	int betacol = 13;
	int secol = 14;
	int impqual1 = 15;
	int impqual2 = 16;
	int width = 640;
	int height = 480;
	boolean onlyIc = false;

	public static void main(String[] args) {
		Plotter p = new Plotter();
		try {
			p.plotCorr();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
	}

	public boolean isIndel1(String[] elems) {
		String alleles = elems[minorallele1col];
		String[] alleleElems = alleles.split(",");
		boolean b = false;
		for (String s : alleleElems) {
			if (s.length() > 1) {
				b = true;
			}
		}
		return b;
	}

	public boolean isIndel2(String[] elems) {
		String alleles = elems[minorallele2col];
		String[] alleleElems = alleles.split(",");
		boolean b = false;
		for (String s : alleleElems) {
			if (s.length() > 1) {
				b = true;
			}
		}
		return b;
	}

	public boolean maf1belowfthreshold(String[] elems, double t) {
		Double d = Double.parseDouble(elems[maf1col]);
		if (d < t) {
			return true;
		} else {
			return false;
		}
	}

	public boolean maf2belowfthreshold(String[] elems, double t) {
		Double d = Double.parseDouble(elems[maf2col]);
		if (d < t) {
			return true;
		} else {
			return false;
		}
	}

	public void plotCorr() throws IOException, DocumentException {
		String[] files = new String[]{
				"D:\\tmp\\2016-05-19\\T1D\\T1D-EUR-merged.txt",
				"D:\\tmp\\2016-05-19\\T1D\\T1D-COSMO-merged.txt",
				"D:\\tmp\\2016-05-19\\T1D\\T1D-HRC-COSMO-merged.txt"
		};

		String[] files2 = new String[]{
				"D:\\tmp\\2016-05-19\\ImmunoChipGenotyped.txt"
		};
//		files = new String[]{
////				"D:\\tmp\\2016-05-19\\T1D\\T1D-EUR-merged.txt",
//				"D:\\tmp\\2016-05-19\\T1D\\T1D-COSMO-merged.txt",
////				"D:\\tmp\\2016-05-19\\T1D\\T1D-HRC-COSMO-merged.txt"
//		};

		String[] labels = new String[]{"EUR", "COSMO", "HRC-COSMO"};
		String variantsOnIC = "D:\\tmp\\2016-05-19\\T1D\\T1D-recode-stats.vcf.gz";

		String out = "";

		HashSet<String> variantHash = null;
		if (variantsOnIC != null) {
			variantHash = loadVariantHash(variantsOnIC);
		}

		boolean includeindels = true;
		boolean usemafthreshold = false;
		boolean requireabovemaf = false;
		boolean plotOnlyImputed = false;
		double mafthreshold = 0.01;

		includeindels = true;
		usemafthreshold = false;
		requireabovemaf = false;
		plotOnlyImputed = false;
		mafthreshold = 0.01;
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot1-rsquared-unfiltered.png";
		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);

		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot2-rsquared-unfiltered.png";
		plot2(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot3-rsquared-unfiltered.png";
		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);

		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot4-betaVsImpQual-unfiltered.png";
		plot4(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot4-betaVsImpQual-imputedonly.png";
		plot4(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, true, variantHash);

		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot4-CorrVsImpQual-unfiltered.png";
		plot5(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot4-CorrVsImpQual-imputedonly.png";
		plot5(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, true, variantHash);

		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot2-immunochip.png";
		onlyIc = true;
		plot2(files2, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		onlyIc = false;


		includeindels = false;
		usemafthreshold = false;
		requireabovemaf = false;
		plotOnlyImputed = false;
		mafthreshold = 0.01;
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot1-rsquared-noindels.png";
		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot2-rsquared-noindels.png";
		plot2(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot3-rsquared-noindels.png";
		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);

		includeindels = false;
		usemafthreshold = true;
		requireabovemaf = true;
		plotOnlyImputed = false;
		mafthreshold = 0.01;
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot1-rsquared-noindels-mafgt0.01.png";
		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot2-rsquared-noindels-mafgt0.01.png";
		plot2(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot3-rsquared-noindels-mafgt0.01.png";
		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);

		includeindels = false;
		usemafthreshold = true;
		requireabovemaf = false;
		plotOnlyImputed = false;
		mafthreshold = 0.01;
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot1-rsquared-noindels-maflt0.01.png";
		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot2-rsquared-noindels-maflt0.01.png";
		plot2(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot3-rsquared-noindels-maflt0.01.png";
		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);

		includeindels = true;
		usemafthreshold = false;
		requireabovemaf = false;
		plotOnlyImputed = true;
		mafthreshold = 0.01;
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot1-rsquared-imputedonly.png";
		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot2-rsquared-imputedonly.png";
		plot2(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot3-rsquared-imputedonly.png";
		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);

		includeindels = false;
		usemafthreshold = false;
		requireabovemaf = false;
		plotOnlyImputed = true;
		mafthreshold = 0.01;
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot1-rsquared-imputedonly-noindels.png";
		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot2-rsquared-imputedonly-noindels.png";
		plot2(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot3-rsquared-imputedonly-noindels.png";
		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);

		includeindels = false;
		usemafthreshold = true;
		requireabovemaf = true;
		plotOnlyImputed = true;
		mafthreshold = 0.01;
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot1-rsquared-imputedonly-noindels-mafgt0.01.png";
		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot2-rsquared-imputedonly-noindels-mafgt0.01.png";
		plot2(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot3-rsquared-imputedonly-noindels-mafgt0.01.png";
		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);

		includeindels = false;
		usemafthreshold = true;
		requireabovemaf = false;
		plotOnlyImputed = true;
		mafthreshold = 0.01;
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot1-rsquared-imputedonly-noindels-maflt0.01.png";
		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot2-rsquared-imputedonly-noindels-maflt0.01.png";
		plot2(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot3-rsquared-imputedonly-noindels-maflt0.01.png";
		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);

		/// impquals
		includeindels = true;
		usemafthreshold = false;
		requireabovemaf = false;
		plotOnlyImputed = false;
		mafthreshold = 0.01;
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot1-impqual-unfiltered.png";
		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot2-impqual-unfiltered.png";
		plot2(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot3-impqual-unfiltered.png";
		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);

		includeindels = false;
		usemafthreshold = false;
		requireabovemaf = false;
		plotOnlyImputed = false;
		mafthreshold = 0.01;
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot1-impqual-noindels.png";
		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot2-impqual-noindels.png";
		plot2(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot3-impqual-noindels.png";
		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);

		includeindels = false;
		usemafthreshold = true;
		requireabovemaf = true;
		plotOnlyImputed = false;
		mafthreshold = 0.01;
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot1-impqual-noindels-mafgt0.01.png";
		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot2-impqual-noindels-mafgt0.01.png";
		plot2(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot3-impqual-noindels-mafgt0.01.png";
		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);

		includeindels = false;
		usemafthreshold = true;
		requireabovemaf = false;
		plotOnlyImputed = false;
		mafthreshold = 0.01;
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot1-impqual-noindels-maflt0.01.png";
		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot2-impqual-noindels-maflt0.01.png";
		plot2(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot3-impqual-noindels-maflt0.01.png";
		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);

		includeindels = true;
		usemafthreshold = false;
		requireabovemaf = false;
		plotOnlyImputed = true;
		mafthreshold = 0.01;
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot1-impqual-imputedonly.png";
		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot2-impqual-imputedonly.png";
		plot2(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot3-impqual-imputedonly.png";
		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);

		includeindels = false;
		usemafthreshold = false;
		requireabovemaf = false;
		plotOnlyImputed = true;
		mafthreshold = 0.01;
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot1-impqual-imputedonly-noindels.png";
		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot2-impqual-imputedonly-noindels.png";
		plot2(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot3-impqual-imputedonly-noindels.png";
		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);

		includeindels = false;
		usemafthreshold = true;
		requireabovemaf = true;
		plotOnlyImputed = true;
		mafthreshold = 0.01;
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot1-impqual-imputedonly-noindels-mafgt0.01.png";
		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot2-impqual-imputedonly-noindels-mafgt0.01.png";
		plot2(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot3-impqual-imputedonly-noindels-mafgt0.01.png";
		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);

		includeindels = false;
		usemafthreshold = true;
		requireabovemaf = false;
		plotOnlyImputed = true;
		mafthreshold = 0.01;
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot1-impqual-imputedonly-noindels-maflt0.01.png";
		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot2-impqual-imputedonly-noindels-maflt0.01.png";
		plot2(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);
		out = "D:\\tmp\\2016-05-19\\plotsacc\\plot3-impqual-imputedonly-noindels-maflt0.01.png";
		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash);

//

//

//
//			// plot 5: imputation qual vs correlation
//			grid = new Grid(1000, 1000, 1, files.length, 0, 0);
//			for (int i = 0; i < files.length; i++) {
//				String file = files[i];
//				int ctr = 0;
//				TextFile tf = new TextFile(file, TextFile.R);
//				String[] elems = tf.readLineElems(TextFile.tab);
//				ArrayList<Double> fy = new ArrayList<>();
//				ArrayList<Double> fx = new ArrayList<>();
//				while (elems != null) {
//					double impqual = Double.parseDouble(elems[1]);
//					double correlation = Double.parseDouble(elems[1]);
//					fx.add(impqual);
//					fy.add(correlation);
//					ctr++;
//					elems = tf.readLineElems(TextFile.tab);
//				}
//				tf.close();
//
//				ScatterplotPanel panel3 = new ScatterplotPanel(1, 1);
//				grid.addPanel(panel3);
//			}
//			grid.draw(out);


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

	public void plot1(String[] files, String[] labels, String out, int col,
	                  boolean includeIndels,
	                  boolean usemafthreshold,
	                  boolean requireabovemaf,
	                  double mafthreshold,
	                  boolean plotOnlyImputed,
	                  HashSet<String> variantHash
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
				boolean isIndel = isIndel2(elems);
				boolean belowmaf = maf2belowfthreshold(elems, mafthreshold);

				boolean include = true;
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
					double val = Double.parseDouble(elems[col]);
					corvals.add(val);
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
		panel.setDatasetLabels(labels);
		grid.addPanel(panel);
		grid.draw(out);
	}

	private boolean isVariantOnIC(String[] elems, HashSet<String> variantHash) {
		String variant = elems[idcol];
		return variantHash.contains(variant);

	}

	public void plot2(String[] files, String[] datasetLabels, String out, int col,
	                  boolean includeIndels,
	                  boolean usemafthreshold,
	                  boolean requireabovemaf,
	                  double mafthreshold,
	                  boolean plotOnlyImputed,
	                  HashSet<String> variantHash
	) throws IOException, DocumentException {
		// plot 2: x-axis maf, y-axis correlation (boxplot)

		String[] labels = new String[]{
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
				boolean include = true;
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
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();

			// convert to double[]
			for (int bin = 0; bin < upperthreshold.length; bin++) {
				ArrayList<Double> d = filebins.get(bin);
				double[] arr = Primitives.toPrimitiveArr(d.toArray(new Double[0]));
				bins[i][bin] = arr;
			}
			System.out.println("Loaded: " + loaded + " from file: " + file);
		}

		Grid grid = new Grid(width, height, 1, 1, 100, 100);
		BoxPlotPanel panel = new BoxPlotPanel(1, 1);
		panel.setData(bins);
		panel.setDrawDataPoints(false);
		panel.useTukeysDefault(true);
		panel.setLabels(labels);
		grid.addPanel(panel);
		grid.draw(out);
	}

	public void plot3(String[] files, String[] datasetLabels, String out, int col,
	                  boolean includeIndels,
	                  boolean usemafthreshold,
	                  boolean requireabovemaf,
	                  double mafthreshold,
	                  boolean plotOnlyImputed,
	                  HashSet<String> variantHash) throws IOException, DocumentException {
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

				boolean include = true;
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
	                  HashSet<String> variantHash) throws IOException, DocumentException {
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

				boolean include = true;
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
	                  HashSet<String> variantHash) throws IOException, DocumentException {
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

				boolean include = true;
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
