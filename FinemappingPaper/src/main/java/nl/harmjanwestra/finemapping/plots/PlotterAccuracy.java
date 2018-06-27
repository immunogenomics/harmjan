package nl.harmjanwestra.finemapping.plots;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.finemapping.genotypesandimputation.VariantCounter;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.panels.BoxPlotPanel;
import nl.harmjanwestra.utilities.graphics.panels.HistogramPanel;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Triple;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import nl.harmjanwestra.utilities.legacy.genetica.util.Primitives;

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
	final double mafthreshold = 0.01;
	final boolean onlyIc = false;
	final boolean windows = false;
	final boolean includeId = false;
	double qualthreshold = 0.5;
	String ext = "png";
	boolean plotticks = false;
	boolean plotlegend = false;
	
	public static void main(String[] args) {
		PlotterAccuracy p = new PlotterAccuracy();
		String[] files = null;
		String[] labels = null;
		String diseaseprefix = "";
		
		String seqpanelvcf = "";
		String variantsOnIC = "";
		String bedregions = "";
		String samplelist = null;
		String outdir = null;
		
		variantsOnIC = "d:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2016-06-21-ImputationQuality\\RAAndT1D-recode-maf0005-ICRegionsW100kb-samplenamefix.vcf.gz-updatedRSId-stats.vcf.gz";
		bedregions = "d:/Sync/SyncThing/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		
		try {
			outdir = "d:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\outputRA\\";
			Gpio.createDir(outdir);
			files = new String[]{
					"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\RA-cosmo-HC.txt",
					"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\RA-PBWT-HC.txt",
					"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\RA-eur-HC.txt",
					"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\RA-HRC-HC.txt",
				
			};
			labels = new String[]{
					"COSMO",
					"COSMO-PBWT",
					"EUR",
					"HRC",
			};
			diseaseprefix = "RA-HC";
			seqpanelvcf = "d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\panels\\hapcaller-maf0005-cr0950-rd10-gq30.vcf.gz-samplenamefix-mixupfix-nonmatchingremoved.vcf.gz";
//			p.plotCorr(files, labels, diseaseprefix, seqpanelvcf, variantsOnIC, bedregions, samplelist, outdir);

//			files = new String[]{
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\RA-cosmo-UG.txt",
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\RA-PBWT-UG.txt",
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\RA-eur-UG.txt",
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\RA-HRC-UG.txt",
//
//			};
//			labels = new String[]{
//					"COSMO",
//					"COSMO-PBWT",
//					"EUR",
//					"HRC",
//
//			};
//			diseaseprefix = "RA-UG";
//			seqpanelvcf = "c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\panels\\unifiedgenotyper-maf0005-cr0950-rd10-gq30.vcf.gz-samplenamefix-mixupfix-nonmatchingremoved.vcf.gz";
//			p.plotCorr(files, labels, diseaseprefix, seqpanelvcf, variantsOnIC, bedregions, samplelist, outdir);
//
//			files = new String[]{
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\RA-cosmo-ST.txt",
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\RA-PBWT-ST.txt",
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\RA-eur-ST.txt",
//					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\RA-HRC-ST.txt",
//
//			};
//			labels = new String[]{
//					"COSMO",
//					"COSMO-PBWT",
//					"EUR",
//					"HRC",
//
//			};
//			diseaseprefix = "RA-ST";
//			seqpanelvcf = "c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\panels\\samtools-maf0005-cr0950-rd10-gq30.vcf.gz-samplenamefix-mixupfix-nonmatchingremoved.vcf.gz";
//			p.plotCorr(files, labels, diseaseprefix, seqpanelvcf, variantsOnIC, bedregions, samplelist, outdir);
////


//////////////////////////////
// T1D

			outdir = "d:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\outputT1Dtest\\";
			Gpio.createDir(outdir);
			files = new String[]{
					"d:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-cosmo-HC.txt",
//					"d:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-PBWT-HC.txt",
//					"d:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-eur-HC.txt",
//					"d:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-HRC-HC.txt",
//					"d:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-Michigan-HRC-HC.txt",
//					"d:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-SHAPEIT-HRC-HC.txt",
			};
			labels = new String[]{
					"COSMO",
//					"COSMO-PBWT",
//					"EUR",
//					"HRC-Sanger",
//					"HRC-Michigan",
//					"HRC-SHAPEIT",

			};
			diseaseprefix = "T1D-HC";
			seqpanelvcf = "d:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\panels\\hapcaller-maf0005-cr0950-rd10-gq30.vcf.gz-samplenamefix-mixupfix-nonmatchingremoved.vcf.gz";
			p.plotCorr(files, labels, diseaseprefix, seqpanelvcf, variantsOnIC, bedregions, samplelist, outdir);

//			files = new String[]{
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-cosmo-UG.txt",
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-PBWT-UG.txt",
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-eur-UG.txt",
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-HRC-UG.txt",
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-Michigan-HRC-UG.txt",
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-SHAPEIT-HRC-UG.txt",
//
//			};
//			labels = new String[]{
//					"COSMO",
//					"COSMO-PBWT",
//					"EUR",
//					"HRC-Sanger",
//					"HRC-Michigan",
//					"HRC-SHAPEIT",
//
//			};
//			diseaseprefix = "T1D-UG";
//			seqpanelvcf = "c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\panels\\unifiedgenotyper-maf0005-cr0950-rd10-gq30.vcf.gz-samplenamefix-mixupfix-nonmatchingremoved.vcf.gz";
//			p.plotCorr(files, labels, diseaseprefix, seqpanelvcf, variantsOnIC, bedregions, samplelist, outdir);
//
//			files = new String[]{
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-cosmo-ST.txt",
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-PBWT-ST.txt",
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-eur-ST.txt",
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-HRC-ST.txt",
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-Michigan-HRC-ST.txt",
//					"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-SHAPEIT-HRC-ST.txt",
//
//			};
//			labels = new String[]{
//					"COSMO",
//					"COSMO-PBWT",
//					"EUR",
//					"HRC-Sanger",
//					"HRC-Michigan",
//					"HRC-SHAPEIT",
//
//			};
//			diseaseprefix = "T1D-ST";
//			seqpanelvcf = "c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\panels\\samtools-maf0005-cr0950-rd10-gq30.vcf.gz-samplenamefix-mixupfix-nonmatchingremoved.vcf.gz";
//			p.plotCorr(files, labels, diseaseprefix, seqpanelvcf, variantsOnIC, bedregions, samplelist, outdir);
//
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
	}
	
	
	public void plotCorr(String[] files, String[] labels, String diseaseprefix, String seqpanelvcf, String variantsOnIC, String bedregions, String samplelist, String outdir) throws IOException, DocumentException {
//		String[] files = new String[]{
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-08-31-Accuracy/T1D-EUR.txt",
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-08-31-Accuracy/T1D-COSMO.txt",
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-08-31-Accuracy/T1D-HRC-EAGLE.txt",
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-08-31-Accuracy/T1D-HRC-SHAPEIT.txt",
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-08-31-Accuracy/T1D-HRC-EAGLE-Michigan.txt"
//		};
//
//		String[] labels = new String[]{
//				"EUR",
//				"COSMO",
//				"HRC-EAGLE",
//				"HRC-SHAPEIT",
//				"HRC-EAGLE-MICHIGAN"
//		};
//		String samplelist = null;
//		String diseaseprefix = "T1D";
//
//		files = new String[]{
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-08-31-Accuracy/RA-EUR.txt",
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-08-31-Accuracy/RA-COSMO.txt",
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-08-31-Accuracy/RA-HRC-w100kb.txt",
//		};
//		labels = new String[]{
//				"EUR",
//				"COSMO",
//				"HRC-EAGLE"
//		};
//		samplelist = null;
//		diseaseprefix = "RA";
//
//
//		String seqpanelvcf = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/seqpanel/seqpanelfiltered-maf0005-cr0950-rd10-gq30-runNamesFixed-RASampleNamesFixed-badSamplesRemoved-mixupsFixed.vcf.gz-updatedRSId.vcf.gz";
//		String variantsOnIC = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/RAAndT1D-recode-maf0005-ICRegionsW100kb-samplenamefix.vcf.gz-updatedRSId-stats.vcf.gz";
//
////		String[] files2 = new String[]{
////				"/Data/tmp/2016-05-20/T1D/ImmunoChipGenotyped.txt"
////		};
//		String bedregions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
////		bedregions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci.bed";
//		bedregions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
//		String outdir = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-08-31-Accuracy/" + diseaseprefix + "-plots/";
//


//		if (windows) {
//
//			files = new String[]{
//					"D:\\tmp\\2016-05-19\\T1D\\T1D-EUR-merged.txt",
//					"D:\\tmp\\2016-05-19\\T1D\\T1D-COSMO-merged.txt",
//					"D:\\tmp\\2016-05-19\\T1D\\T1D-HRC-COSMO-merged.txt",
//					"D:\\tmp\\2016-05-19\\T1D\\T1D-HRC-HRC-w100kb-merged.txt",
//					"D:\\tmp\\2016-05-19\\T1D\\T1D-HRC-COSMO-w100kb-merged.txt"
//			};
//			labels = new String[]{"EUR", "COSMO", "HRC-COSMO", "HRC-HRC-w100kb", "HRC-COSMO-w100kb"};
//			variantsOnIC = "D:\\tmp\\2016-05-19\\T1D-recode-stats.vcf.gz";
//			bedregions = "D:\\tmp\\2016-05-19\\AllICLoci.bed";
//			outdir = "D:\\tmp\\2016-05-19\\T1D-plotsAccuracy\\";
////			files2 = new String[]{
////					"D:\\tmp\\2016-05-19\\ImmunoChipGenotyped.txt"
////			};
//		}
		
		
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
		
		
		// count variants in reference
		
		
		// determine number of imputed variants
		
		
		// including indels
//		int maxNrVariants = 3000;
		includeindels = true;
		usemafthreshold = false;
		includeICVariants = true;
		out = outdir + diseaseprefix + "-plot1-allvariants." + ext;
		
		plot1(files, labels, out, bedregions, includeindels, usemafthreshold, 0, includeICVariants, includeId, variantsOnIC, seqpanelvcf, samplelist);
		
		// maf > 1%
		usemafthreshold = true;
		out = outdir + diseaseprefix + "-plot1-allvariants-mafgt" + mafthreshold + "." + ext;
		plot1(files, labels, out, bedregions, includeindels, usemafthreshold, mafthreshold, includeICVariants, includeId, variantsOnIC, seqpanelvcf, samplelist);
		
		
		// maf > 1%, imputed variants
		usemafthreshold = true;
		includeICVariants = false;
//		maxNrVariants = 692; // if including IDs
		out = outdir + diseaseprefix + "-plot1-allvariants-imputed-mafgt" + mafthreshold + "." + ext;
		plot1(files, labels, out, bedregions, includeindels, usemafthreshold, mafthreshold, includeICVariants, includeId, variantsOnIC, seqpanelvcf, samplelist);
		
	}
	
	
	public void plot1(String[] files, String[] labels, String out, String bedregionsfile,
					  boolean includeIndels,
					  boolean usemafthreshold,
					  double mafthreshold,
					  boolean includeICVariants,
					  boolean includeId,
					  String variantsOnIC,
					  String seqpanelvcf,
					  String samplelist
	
	) throws IOException, DocumentException {
		// plot 1: x-axis nr of variants, y-axis correlation,
		ArrayList<ArrayList<Double>> vals = new ArrayList<ArrayList<Double>>();
		
		System.out.println("Using MAF threshold: " + usemafthreshold);
		
		double upperthreshold = 1d;
		HashSet<String> variantsOnICHash = loadVariantHash(variantsOnIC, includeId);
		System.out.println(variantsOnICHash.size() + " total on IC");
		Triple<ArrayList<VCFVariant>, ArrayList<VCFVariant>, ArrayList<VCFVariant>> seqpanelvariants = loadSequencedVariants(
				seqpanelvcf, bedregionsfile, 0.01, upperthreshold, variantsOnICHash, includeId, includeIndels, samplelist
		);
		
		ArrayList<VCFVariant> seqpanel = seqpanelvariants.getLeft();
		
		System.out.println("Setting max size: " + seqpanelvcf);
		ArrayList<VCFVariant> variantsOnImmunoChip = seqpanelvariants.getMiddle();
		ArrayList<VCFVariant> variantsNotOnImmunoChip = seqpanelvariants.getRight();
		
		System.out.println(seqpanel.size() + " variants in VCF");
		System.out.println(variantsNotOnImmunoChip.size() + " not on IC");
		System.out.println(variantsOnImmunoChip.size() + " on IC");
		
		
		HashSet<String> sequencedVariantsHash = null;
		ArrayList<VCFVariant> allvars = new ArrayList<>();
		
		if (includeICVariants) {
			allvars.addAll(variantsNotOnImmunoChip);
			allvars.addAll(variantsOnImmunoChip);
			sequencedVariantsHash = hashit(allvars, includeId);
		} else {
			sequencedVariantsHash = hashit(variantsNotOnImmunoChip, includeId);
		}
		
		int maxSize = sequencedVariantsHash.size();
		
		String[] newlabels = new String[labels.length];
		for (int f = 0; f < files.length; f++) {
			// get the imputation accuracies for these variants
			TextFile tf2 = new TextFile(files[f], TextFile.R);
//			System.out.println("parsing: " + files[f]);
			ArrayList<Double> corvals = new ArrayList<>();
			String[] elems = tf2.readLineElems(TextFile.tab);
			int nrAboveQualThreshold = 0;
			
			ArrayList<String> missing = new ArrayList<String>();
			
			
			ArrayList<String> variantsSeen = new ArrayList<>();
			ArrayList<String[]> variantsBelowThresholds = new ArrayList<>();
			ArrayList<String> variantsAboveThresholds = new ArrayList<>();
			while (elems != null) {
				if (elems.length > 12) {
					if (!elems[rsqlcol].equals("null")) {
						double val = Double.parseDouble(elems[rsqlcol]);
						double maf = Double.parseDouble(elems[maf2col]);
						
						String[] varElems = elems[0].split("_");
						
						boolean sequenced = isVariantInHash(varElems, sequencedVariantsHash, includeId);
						
						boolean abovethreshold = true;
						if (sequenced) {
							
							if (!usemafthreshold || (usemafthreshold && maf > mafthreshold)) {
								corvals.add(val);
								if (val > qualthreshold) {
									nrAboveQualThreshold++;
									
								} else {
									abovethreshold = false;
								}
							} else {
								abovethreshold = false;
							}
							
							if (!abovethreshold) {
								variantsBelowThresholds.add(elems);
							} else {
//							String variant = Chromosome.parseChr(elems[0]).toString() + "_" + elems[1] + "_" + elems[2];
								variantsAboveThresholds.add(elems[0]);
							}
							
							
						}
						
						variantsSeen.add(elems[0]);
						
					}
				}
				elems = tf2.readLineElems(TextFile.tab);
			}
			tf2.close();
			Collections.sort(corvals, Collections.reverseOrder());
			System.out.println(nrAboveQualThreshold + " > " + qualthreshold + "\t" + corvals.size() + " vals in path " + files[f]);
			newlabels[f] = labels[f] + " - " + corvals.size() + " / " + maxSize;
			vals.add(corvals);
			
			
			TextFile outf = new TextFile(out + labels[f] + "-belowthresholds.txt", TextFile.W);
			System.out.println("below thresh:\t" + outf.getFullPath());
			for (String[] elem : variantsBelowThresholds) {
				outf.writeln(Strings.concat(elem, Strings.tab));
			}
			outf.close();
			
			HashSet<String> variantsAboveThresholdHash = new HashSet<String>();
			variantsAboveThresholdHash.addAll(variantsAboveThresholds);
			System.out.println(variantsAboveThresholds.size() + " variants above threshold");
			TextFile outf2 = new TextFile(out + labels[f] + "-missing.txt", TextFile.W);
			System.out.println("missing var:\t" + outf2.getFullPath());
			int missingnr = 0;
			for (VCFVariant var : allvars) {
				String variant = var.getChrObj().getNumber() + "_" + var.getPos() + "_" + var.getId();
				if (!variantsAboveThresholdHash.contains(variant)) {
					outf2.writeln(variant + "\t" + var.getMinorAllele() + "\t" + Strings.concat(var.getAlleles(), Strings.comma) + "\t" + var.getMAF());
					missingnr++;
				}
			}
			System.out.println(missingnr + " missing");
			outf2.close();
			
			
			TextFile outf3 = new TextFile(out + labels[f] + "-neverimputed.txt", TextFile.W);
			System.out.println("missing var:\t" + outf3.getFullPath());
			HashSet<String> variantsImputedHash = new HashSet<String>();
			variantsImputedHash.addAll(variantsSeen);
			missingnr = 0;
			for (VCFVariant var : allvars) {
				String variant = var.getChrObj().getNumber() + "_" + var.getPos() + "_" + var.getId();
				if (!variantsImputedHash.contains(variant)) {
					outf3.writeln(variant + "\t" + var.getMinorAllele() + "\t" + Strings.concat(var.getAlleles(), Strings.comma) + "\t" + var.getMAF());
					missingnr++;
				}
			}
			System.out.println(missingnr + " missing");
			outf3.close();
			
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
		panel.setPlotElems(plotticks, plotlegend);
		Range range = new Range(0, 0, 1, 100);
		range.roundX();
		range.roundY();
		panel.setDataRange(range);
		panel.setDatasetLabels(newlabels);
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
					isOnIc = isVariantInHash(elems, variantHash, true);
					if (plotOnlyImputed && isOnIc) {
						include = false;
					} else if (!plotOnlyImputed && !isOnIc) {
						include = false;
					}
				}
				
				if (onlyIc && variantHash != null) {
					boolean isOnIc = false;
					isOnIc = isVariantInHash(elems, variantHash, true);
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
			System.out.println("Loaded: " + loaded + " from path: " + file);
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
					boolean isOnIc = isVariantInHash(elems, variantHash, true);
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
					boolean isOnIc = isVariantInHash(elems, variantHash, true);
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
					boolean isOnIc = isVariantInHash(elems, variantHash, true);
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
