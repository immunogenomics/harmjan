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
 * Created by hwestra on 5/18/16.
 */
public class PlotterAccuracy {


	// id	allele1	allele2	maf	callrate	allele21	allele22	maf 	cr	df	nrsamples	r	rsq	beta	se	impqual0	impqual1
	int idcol = 0;
	int minorallele1col = 1;
	int aleleles1col = 2;
	int maf1col = 3;
	int cr1col = 4;
	int minorallele2col = 5;
	int aleleles2col = 6;
	public int maf2col = 7;
	int cr2col = 8;
	int dfcol = 9;
	int samplecol = 10;
	int rcol = 11;
	public int rsqlcol = 12;
	int betacol = 13;
	int secol = 14;
	int impqual1 = 15;
	int impqual2 = 16;
	int width = 640;
	int height = 480;
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
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/Accuracy/T1D-EUR.txt",
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/Accuracy/T1D-COSMO.txt",
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/Accuracy/T1D-HRC-COSMO.txt",
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/Accuracy/T1D-HRC-w100kb.txt"
		};
		String diseaseprefix = "T1D";

//		files = new String[]{
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/Accuracy/RA-EUR.txt",
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/Accuracy/RA-COSMO.txt",
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/Accuracy/RA-HRC-COSMO.txt",
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/Accuracy/RA-HRC-w100kb.txt"
//		};
//		diseaseprefix = "RA";

		String[] labels = new String[]{"EUR", "COSMO", "HRC-COSMO", "HRC-COSMO-w100kb"};
		String variantsOnIC = "/Data/tmp/2016-05-20/T1D-recode-stats.vcf.gz";

		String[] files2 = new String[]{
				"/Data/tmp/2016-05-20/T1D/ImmunoChipGenotyped.txt"
		};
		String bedregions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		bedregions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci.bed";
		String outdir = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/AccuracyPlots/";
		double mafthreshold = 0.0;

		boolean windows = false;

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
			files2 = new String[]{
					"D:\\tmp\\2016-05-19\\ImmunoChipGenotyped.txt"
			};
		}

		String ext = "pdf";


		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> bedfileRegions = reader.readAsList(bedregions);

		if (!Gpio.exists(outdir)) {
			Gpio.createDir(outdir);
		}


		String out = "";

		HashSet<String> variantHash = null;
		if (variantsOnIC != null) {
			variantHash = loadVariantHash(variantsOnIC);
		}

		boolean includeindels = true;
		boolean usemafthreshold = false;
		boolean requireabovemaf = false;
		boolean plotOnlyImputed = false;


		includeindels = true;
		usemafthreshold = true;
		requireabovemaf = true;
		plotOnlyImputed = true;
		out = outdir + diseaseprefix + "-plot1-imputedonly-rsquared-withindels-mafgt" + mafthreshold + "." + ext;
		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);

		includeindels = false;
		usemafthreshold = true;
		requireabovemaf = true;
		plotOnlyImputed = true;
		out = outdir + diseaseprefix + "-plot1-imputedonly-rsquared-withoutindels-mafgt" + mafthreshold + "." + ext;
		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);

		includeindels = true;
		usemafthreshold = true;
		requireabovemaf = true;
		plotOnlyImputed = true;
		out = outdir + diseaseprefix + "-plot1-imputedonly-impqual-withindels-mafgt" + mafthreshold + "." + ext;
		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);

		includeindels = false;
		usemafthreshold = true;
		requireabovemaf = true;
		plotOnlyImputed = true;
		out = outdir + diseaseprefix + "-plot1-imputedonly-impqual-withoutindels-mafgt" + mafthreshold + "." + ext;
		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);


		// include unimputed variants as well
		includeindels = true;
		usemafthreshold = true;
		requireabovemaf = true;
		plotOnlyImputed = false;
		out = outdir + diseaseprefix + "-plot1-rsquared-withindels-mafgt" + mafthreshold + "." + ext;
		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);

		includeindels = false;
		usemafthreshold = true;
		requireabovemaf = true;
		plotOnlyImputed = false;
		out = outdir + diseaseprefix + "-plot1-rsquared-withoutindels-mafgt" + mafthreshold + "." + ext;
		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);

		includeindels = true;
		usemafthreshold = true;
		requireabovemaf = true;
		plotOnlyImputed = false;
		out = outdir + diseaseprefix + "-plot1-impqual-withindels-mafgt" + mafthreshold + "." + ext;
		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);

		includeindels = false;
		usemafthreshold = true;
		requireabovemaf = true;
		plotOnlyImputed = false;
		out = outdir + diseaseprefix + "-plot1-impqual-withoutindels-mafgt" + mafthreshold + "." + ext;
		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);


//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
// out = outdir + "plot1-rsquared-unfiltered." + ext;
//		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		out = outdir + "plot2-rsquared-unfiltered." + ext;
//		plot2(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot3-rsquared-unfiltered." + ext;
//		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		out = outdir + "plot4-betaVsImpQual-unfiltered." + ext;
//		plot4(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot4-betaVsImpQual-imputedonly." + ext;
//		plot4(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, true, variantHash, bedfileRegions);
//
//		out = outdir + "plot4-CorrVsImpQual-unfiltered." + ext;
//		plot5(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot4-CorrVsImpQual-imputedonly." + ext;
//		plot5(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, true, variantHash, bedfileRegions);
//
//		out = outdir + "plot2-immunochip." + ext;
//		onlyIc = true;
//		plot2(files2, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		onlyIc = false;
//
//
//		includeindels = false;
//		usemafthreshold = false;
//		requireabovemaf = false;
//		plotOnlyImputed = false;
//
//		out = outdir + "plot1-rsquared-noindels." + ext;
//		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-rsquared-noindels." + ext;
//		plot2(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot3-rsquared-noindels." + ext;
//		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = false;
//		usemafthreshold = true;
//		requireabovemaf = true;
//		plotOnlyImputed = false;
//
//		out = outdir + "plot1-rsquared-noindels-mafgt" + mafthreshold + "." + ext;
//		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-rsquared-noindels-mafgt" + mafthreshold + "." + ext;
//		plot2(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot3-rsquared-noindels-mafgt" + mafthreshold + "." + ext;
//		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = false;
//		usemafthreshold = true;
//		requireabovemaf = false;
//		plotOnlyImputed = false;
//
//		out = outdir + "plot1-rsquared-noindels-maflt" + mafthreshold + "." + ext;
//		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-rsquared-noindels-maflt" + mafthreshold + "." + ext;
//		plot2(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot3-rsquared-noindels-maflt" + mafthreshold + "." + ext;
//		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = true;
//		usemafthreshold = false;
//		requireabovemaf = false;
//		plotOnlyImputed = true;
//
//		out = outdir + "plot1-rsquared-imputedonly." + ext;
//		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-rsquared-imputedonly." + ext;
//		plot2(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot3-rsquared-imputedonly." + ext;
//		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = false;
//		usemafthreshold = false;
//		requireabovemaf = false;
//		plotOnlyImputed = true;
//
//		out = outdir + "plot1-rsquared-imputedonly-noindels." + ext;
//		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-rsquared-imputedonly-noindels." + ext;
//		plot2(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot3-rsquared-imputedonly-noindels." + ext;
//		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = false;
//		usemafthreshold = true;
//		requireabovemaf = true;
//		plotOnlyImputed = true;
//
//		out = outdir + "plot1-rsquared-imputedonly-noindels-mafgt" + mafthreshold + "." + ext;
//		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-rsquared-imputedonly-noindels-mafgt" + mafthreshold + "." + ext;
//		plot2(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot3-rsquared-imputedonly-noindels-mafgt" + mafthreshold + "." + ext;
//		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = false;
//		usemafthreshold = true;
//		requireabovemaf = false;
//		plotOnlyImputed = true;
//
//		out = outdir + "plot1-rsquared-imputedonly-noindels-maflt" + mafthreshold + "." + ext;
//		plot1(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-rsquared-imputedonly-noindels-maflt" + mafthreshold + "." + ext;
//		plot2(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot3-rsquared-imputedonly-noindels-maflt" + mafthreshold + "." + ext;
//		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		/// impquals
//		includeindels = true;
//		usemafthreshold = false;
//		requireabovemaf = false;
//		plotOnlyImputed = false;
//
//		out = outdir + "plot1-impqual-unfiltered."+ext;
//		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-impqual-unfiltered."+ext;
//		plot2(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot3-impqual-unfiltered."+ext;
//		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = false;
//		usemafthreshold = false;
//		requireabovemaf = false;
//		plotOnlyImputed = false;
//
//		out = outdir + "plot1-impqual-noindels."+ext;
//		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-impqual-noindels."+ext;
//		plot2(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot3-impqual-noindels."+ext;
//		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = false;
//		usemafthreshold = true;
//		requireabovemaf = true;
//		plotOnlyImputed = false;
//
//		out = outdir + "plot1-impqual-noindels-mafgt" + mafthreshold + "."+ext;
//		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-impqual-noindels-mafgt" + mafthreshold + "."+ext;
//		plot2(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot3-impqual-noindels-mafgt" + mafthreshold + "."+ext;
//		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = false;
//		usemafthreshold = true;
//		requireabovemaf = false;
//		plotOnlyImputed = false;
//
//		out = outdir + "plot1-impqual-noindels-maflt" + mafthreshold+ "."+ext;
//		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-impqual-noindels-maflt" + mafthreshold + "."+ext;
//		plot2(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot3-impqual-noindels-maflt" + mafthreshold + "."+ext;
//		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = true;
//		usemafthreshold = false;
//		requireabovemaf = false;
//		plotOnlyImputed = true;
//
//		out = outdir + "plot1-impqual-imputedonly."+ext;
//		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-impqual-imputedonly."+ext;
//		plot2(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot3-impqual-imputedonly."+ext;
//		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = false;
//		usemafthreshold = false;
//		requireabovemaf = false;
//		plotOnlyImputed = true;
//
//		out = outdir + "plot1-impqual-imputedonly-noindels."+ext;
//		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-impqual-imputedonly-noindels."+ext;
//		plot2(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot3-impqual-imputedonly-noindels."+ext;
//		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = false;
//		usemafthreshold = true;
//		requireabovemaf = true;
//		plotOnlyImputed = true;
//
//		out = outdir + "plot1-impqual-imputedonly-noindels-mafgt" + mafthreshold + "."+ext;
//		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-impqual-imputedonly-noindels-mafgt" + mafthreshold + "."+ext;
//		plot2(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot3-impqual-imputedonly-noindels-mafgt" + mafthreshold + "."+ext;
//		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//
//		includeindels = false;
//		usemafthreshold = true;
//		requireabovemaf = false;
//		plotOnlyImputed = true;
//
//		out = outdir + "plot1-impqual-imputedonly-noindels-maflt" + mafthreshold + "."+ext;
//		plot1(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot2-impqual-imputedonly-noindels-maflt" + mafthreshold + "."+ext;
//		plot2(files, labels, out, impqual2, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);
//		out = outdir + "plot3-impqual-imputedonly-noindels-maflt" + mafthreshold + "."+ext;
//		plot3(files, labels, out, rsqlcol, includeindels, usemafthreshold, requireabovemaf, mafthreshold, plotOnlyImputed, variantHash, bedfileRegions);

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
		panel.setData(y, x);
		Range range = new Range(0, 0, 1, maxSize);
		range.roundX();
		panel.setDataRange(range);
		panel.setDatasetLabels(labels);
		grid.addPanel(panel);
		grid.draw(out);
	}

	private boolean isVariantOnIC(String[] elems, HashSet<String> variantHash) {
		String variant = elems[idcol];
		return variantHash.contains(variant);
	}

	private boolean isWithinRegion(ArrayList<Feature> list, String[] elems) {
		Feature varfeat = new Feature();
		String[] idelems = elems[idcol].split("_");

		int pos = Integer.parseInt(idelems[1]);
		varfeat.setChromosome(Chromosome.parseChr(idelems[0]));
		varfeat.setStart(pos);
		varfeat.setStop(pos);
		for (Feature f : list) {
			if (f.overlaps(varfeat)) {
				return true;
			}
		}
		return false;
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
