package nl.harmjanwestra.ngs;

import com.lowagie.text.DocumentException;
import nl.harmjanwestra.ngs.GenotypeFormats.PedAndMapFunctions;
import nl.harmjanwestra.ngs.GenotypeFormats.VCF.VCFGenotypeData;
import nl.harmjanwestra.ngs.GenotypeFormats.VCFFunctions;
import nl.harmjanwestra.ngs.GenotypeFormats.VCFVariant;
import nl.harmjanwestra.ngs.graphics.VariantPlot;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Created by hwestra on 4/21/15.
 */
public class PlayGround {

	public static void main(String[] args) {

		PlayGround g = new PlayGround();
//		g.prepareImputationPanels();
		// g.prepareGWASDatasets();

//		Integer chrint = Integer.parseInt(args[0]);
//		String vcfsort = args[1];
//		String refvcf = args[2];
//		String testvcf = args[3];
//		String outdir = args[4];
//		boolean phased = Boolean.parseBoolean(args[5]);
		try {
//			// phased
//			if (phased) {
//				g.mergeAndIntersect(chrint, vcfsort, refvcf, testvcf, outdir, "|");
//			} else {
//				//unphased
//				g.mergeAndIntersect(chrint, vcfsort, refvcf, testvcf, outdir, "/");
//
//			}
//			g.immunochipMatchJobs();
//			g.immunoChipImputationjobs();

//			g.matchfiles("/Data/ImmunoChip/T1D/binary/UK.cov", "/Data/ImmunoChip/T1D/binary/UK.fam");
//			TextFile tf1 = new TextFile("/Data/ImmunoChip/T1D/binary/chr10test/plinkassoc.ped", TextFile.R);
//			String[] elems = tf1.readLineElems(Strings.whitespace);
//			while (elems != null) {
//
//				if (!elems[5].equals("1") && !elems[5].equals("2")) {
//					System.out.println(Strings.concat(elems, Strings.tab, 0, 6));
//				}
//				elems = tf1.readLineElems(Strings.whitespace);
//			}
//			tf1.close();
//			System.exit(-1);

//			g.immunoChipCovariates();
//			Integer nrthreads = Integer.parseInt(args[0]);
//			g.testVCF(nrthreads);


//			g.testVCFLocal(1);
//			g.immunoChipRSquaredPlots();
//


//			String cytoband = "/Data/Annotation/UCSC/cytoBand.txt";
//			ChromosomePlot plot = new ChromosomePlot("/Data/tmp/chrplot.pdf", 2400, 1400);
//			plot.setMargin(200);
//			String[] gfffiles = new String[]{"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/CandidateRegions/Hs_GRCh38-RA-assoc_genesGFF",
//					"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/CandidateRegions/Hs_GRCh38-T1D-assoc_genesGFF"};
//			String sequencedRegionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/2015-03-30-allRegionsWith50PercentOfSamplesAbove20x.bed";
//			plot.plot(cytoband, gfffiles, true, sequencedRegionFile);

//			VCFFunctions v = new VCFFunctions();
//			v.summarizeVCF("/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-03-30-CalledGenotypes/merged.vcf", "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-03-30-CalledGenotypes/merged.summary.txt");


//			g.splitVCFs();
//			g.imputeImmunoChipBatches();
//			g.testVCFLocal(4);
//			g.testVCFLocal(8);
//			String prefix = args[0];
//			String outfilename = args[1];
//			Integer nrBatches = Integer.parseInt(args[2]);
//			Integer chr = Integer.parseInt(args[3]);
//			g.mergeImputedBatches(prefix, outfilename, nrBatches, chr);

//			VCFFunctions v = new VCFFunctions();
//			System.out.println(v.getVCFSamples("/Data/ImmunoChip/T1D/merged/merged-Chr1.vcf.gz").size());

			if (args.length < 4) {
				System.out.println("usage: d r theads pseudo");
				System.out.println("d0 = T1D");
				System.out.println("d1 = RA");
				System.out.println("before kg seq kgseq seqvaronly");
			} else {
				g.testVCF(Integer.parseInt(args[0]), Integer.parseInt(args[1]), Integer.parseInt(args[2]), Boolean.parseBoolean(args[3]));
				//"before", "kg", "seq", "kgseq", "seqvaronly"
			}

//			VCFFunctions v = new VCFFunctions();
//			v.isolateGT(args[0], args[1]);

//			g.testVCF(64, false);
//			g.testVCF(64, true);
//			g.testVCF(64);
//			g.immunoChipAssociationPlots();
//			Integer chring = Integer.parseInt(args[0]);
//			String vcfsort = "/Data/Tools/vcftools/bin/vcf-sort";
//			String ref = args[1];
//			String test = args[2];
//			String out = args[3];
//			g.mergeAndIntersect(false, chring, vcfsort, ref, test, out, "/");

//			g.immunoChipPhaseJobsSequencedVariantsOnly();
//

//			g.panelSortAndConcatJobs();
//
//			g.immunoChipImputationjobs();
//			g.mergeImmunoChipVCF();
//
//			g.prepareImputationPanels();
//			g.panelImputationJobs();
//			g.rsquaredfilter(args[0], args[1]);
//			g.immunochipMatchJobs();
//			g.panelMergeJobs();
//			g.datasetmerger();
//			g.panelSortAndConcatJobs();
////			g.prepareGWASDatasets();

//			String dataset1 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-26-CD28-2/4-ICAndSeqVariantMerged/merged-Chr2.vcf.gz";
//			String dataset2 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19variantsinsequencedregions/merge";
//			String out = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-26-CD28-2/9-ImputationPanelsComparedToUnimputed/";
//			g.compareGenotypesBetweenDatasets(dataset1, dataset2, out);

//			g.filterMarkersFile("/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-26-CD28-2/6-PanelsMatched/1kg-matched-sorted.markers",
//					"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-26-CD28-2/6-PanelsMatched/seq-matched.markers",
//					"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-26-CD28-2/6-PanelsMatched/1kg-matched-sorted-unique.markers");

//			g.countVariantsWithR2("/Data/tmp/tmp/test-phased-refimputed-Chr2.vcf.gz",
//					null, "/Data/tmp/tmp/chr2-impsummary.txt", 50);
//			g.countVariantsWithR2("/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-26-CD28-WithNewSettings-EUROnly/7-PanelsImputed/seqpanel-1kgimputed-Chr2.vcf.gz",
//					null,"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-26-CD28-WithNewSettings-EUROnly/7-PanelsImputed/seqpanel-1kgimputed-Chr2-r2info.txt",10);
//

//			g.countVariantsWithR2("/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-28-Chr2-WithNewSettings-EurOnly/7-PanelsImputed/seqpanel-1kgimputed-Chr2.vcf.gz",
//					null, "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-28-Chr2-WithNewSettings-EurOnly/7-PanelsImputed/seqpanel-1kgimputed-Chr2-r2values.txt", 50);


//			String sequencedRegionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/2015-03-30-allRegionsWith50PercentOfSamplesAbove20x.bed";
//			String before = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-28-Chr2-WithNewSettings-EurOnly/6-PanelsMatched/seq-matched-sorted-Chr2.vcf.gz";
//			String after = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-28-Chr2-WithNewSettings-EurOnly/7-PanelsImputed/seqpanel-1kgimputed-Chr2.vcf.gz";
//			String regionfile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-28-Chr2-WithNewSettings-EurOnly/chr2.bed";
//			String outputdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-28-Chr2-WithNewSettings-EurOnly/7-PanelsImputed/imputationr2plots/";

//			String sequencedRegionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/2015-03-30-allRegionsWith50PercentOfSamplesAbove20x.bed";
//			String before = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-28-Chr2-WithNewSettings-EurOnly/6-PanelsMatched/1kg-matched-Chr2.vcf.gz";
//			String after = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-28-Chr2-WithNewSettings-EurOnly/7-PanelsImputed/1kg-seqpanelimputed-Chr2.vcf.gz";
//			String regionfile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-28-Chr2-WithNewSettings-EurOnly/chr2.bed";
//			String outputdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-28-Chr2-WithNewSettings-EurOnly/7-PanelsImputed/imputationr2plots-1kg/";
//			g.regionR2Plots(regionfile, outputdir, before, after, sequencedRegionFile, 0);

//			g.compareGenotypesBetweenDatasets("/Data/Ref/BeagleRef/1kg.phase3.v5-filtered.merged.vcf.gz",
//					"/Data/Ref/BeagleRef/chr22.1kg.phase3.v5-filtered.vcf.gz",
//					"/Data/Ref/BeagleRef/comp/1kg-mergedvs-chr22");

//			g.filterEuropean("/Data/Ref/BeagleRef/1kg.phase3.v5-filtered.merged.vcf.gz","/Data/Ref/BeagleRef/1kg.phase3.v5-filtered.eur.merged.vcf.gz","/Data/Ref/BeagleRef/europeanpopulations.txt");

//			String filter = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/T1D-MatchFiles/kg/";
//			String rakg = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/Assoc/T1D/kg-pseudo/";
//			String rakgout = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/Assoc-OnlyImputed/T1D/kg-pseudo/";
//			Gpio.createDir(rakgout);
//			g.filterAssoc(rakg, filter, rakgout);
//

//			filter = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-30-AllChr-EurOnly/1-SeqVariantFilter/filtered-mendelianerrorsremoved.vcf";
//			String raseq = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/Assoc/T1D/seqvaronly-pseudo/";
//			String raseqout = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/Assoc-OnlyImputed/T1D/seqvaronly-pseudo/";
//
//			Gpio.createDir(raseqout);
//			g.filterAssoc(raseq, filter, raseqout);
//
//			g.immunoChipRSquaredPlots();
//
//			g.immunoChipAssociationPlots();


//			String assocdir1 = "/Sync/Dropbox/PlotsWoICStudy/AssociationResults/T1D/seqvaronly/";
//			String assocdir2 = "/Sync/Dropbox/PlotsWoICStudy/AssociationResults/T1D/kg-pseudo/";
//			g.matchAssoc(assocdir1, assocdir2);
//
//			String rsqdir1 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/T1D/seqvaronly/";
//			String rsqdir2 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/T1D/kg/";
//			String out = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/T1D/rsqcompare.txt";
//			g.matchRSq(rsqdir1, rsqdir2, out);


//			for (int chr = 1; chr < 23; chr++) {
//				String[] files = new String[]{"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/Assoc/RA/before/Chr" + chr + "-assoc.txt",
//						"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/Assoc/RA/kg/Chr" + chr + "-assoc.txt",
//						"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/Assoc/RA/seqvaronly/Chr" + chr + "-assoc.txt"};
//				String[] names = new String[]{"before", "kg", "seqvaronly"};
//				String[] impquals = new String[]{null, "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/RA/kg/merged-variantinfo-Chr" + chr + ".txt",
//						"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/RA/seqvaronly/merged-variantinfo-Chr" + chr + ".txt"};
//				String out = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationAndQualsMerged/RA-merged-Chr" + chr + ".txt";
//
//				g.mergeAssociationResults(names, files, impquals, out);
//			}


//			String vcf = "/Data/tmp/kg/merged-Chr1.vcf.gz";
//			String outdir = "/Data/tmp/kg/";
//			String diseasestatus = "/Data/ImmunoChip/RA/covar/covarmerged.txtmergeddisease.txt";
//			String covar = "/Data/ImmunoChip/RA/covar/covarmerged.txtmergedCovariates.txt";
//			String fam = "/Data/ImmunoChip/RA/binary/RA.fam";
//
//			HashSet<String> snplimit = null; //new HashSet<String>();
////			snplimit.add("rs1028182");
//
//			LogitTestR test = new LogitTestR(vcf, outdir, diseasestatus, covar, null, snplimit, null, fam, false, 0.0, 0.0, 1);
//			test.call();

		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
////			g.compare();

//		VCFFunctions v = new VCFFunctions();
//		String vcfIn = "/Data/Ref/BeagleRef/1kg.phase3.v5-filtered.merged.vcf.gz";
//		String vcfOut = "";
//		try {
//			v.rewriteVariantsAtSamePositionAsMultiAllelic(vcfIn, vcfOut);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}


		// panelSortAndConcat(int chrint, String imputedprefix, String outprefix, String unimputedprefix, String vcfsort, String vcfConcat)
//		try {
//			Integer chrint = Integer.parseInt(args[0]);
//			String imputedprefix = args[1];
//			String outprefix = args[2];
//			String unimputedprefix = args[3];
//			String vcfsort = args[4];
//			String vcfConcat = args[5];
//
//			g.panelSortAndConcat(chrint, imputedprefix, outprefix, unimputedprefix, vcfsort, vcfConcat);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
	}

	private void matchFreq(String freq, String assoc, String log) throws IOException {
		TextFile tf = new TextFile(freq, TextFile.R);

		String ln = tf.readLine();
		HashMap<String, Double> plink = new HashMap<String, Double>();
		while (ln != null) {
			//    1   rs1028182    T    C       0.2116    54430

			while (ln.contains("  ")) {
				ln = ln.replaceAll("  ", " ");
			}
			String[] elems = Strings.whitespace.split(ln);

			String snp = elems[1];
			Double d = Double.parseDouble(elems[4]);
			plink.put(snp, d);

			ln = tf.readLine();
		}
		tf.close();

		TextFile tf2 = new TextFile(assoc, TextFile.R);

		String[] elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[0];
			Double d = Double.parseDouble(elems[3]);
			Double otherD = plink.get(snp);
			if (otherD != null) {
				System.out.println(snp + "\t" + d + "\t" + otherD + "\t" + (d - otherD));
			}
			elems = tf2.readLineElems(TextFile.tab);
		}

		tf2.close();


		// variant genotypes aliased: 0.0 below threshold rs181610198      1       2153588 0.5
		TextFile tf3 = new TextFile(log, TextFile.R);
		elems = tf3.readLineElems(Strings.whitespace);
		while (elems != null) {


			elems = tf3.readLineElems(Strings.whitespace);
		}
//
		tf3.close();
	}


	private HashMap<Feature, Double> loadAssoc(String dir) throws IOException {
		HashMap<Feature, Double> assocValues = new HashMap<Feature, Double>();
		for (int chr = 1; chr < 23; chr++) {
			String f = dir + "Chr" + chr + "-assoc.txt";
			if (Gpio.exists(f)) {

				TextFile tf = new TextFile(f, TextFile.R);
				tf.readLine();
				String[] elems = tf.readLineElems(TextFile.tab);
				while (elems != null) {
					Chromosome chrom = Chromosome.parseChr(elems[1]);
					Feature feat = new Feature();
					feat.setChromosome(chrom);
					feat.setStart(Integer.parseInt(elems[2]));
					feat.setStop(feat.getStart());
					feat.setName(elems[0]);
					assocValues.put(feat, Double.parseDouble(elems[elems.length - 1]));
					elems = tf.readLineElems(TextFile.tab);
				}
				tf.close();
			}
		}
		return assocValues;
	}

	private HashMap<Feature, Double> loadRsq(String dir) throws IOException {
		HashMap<Feature, Double> assocValues = new HashMap<Feature, Double>();
		for (int chr = 1; chr < 23; chr++) {
			String f = dir + "merged-variantinfo-Chr" + chr + ".txt";
			if (Gpio.exists(f)) {

				TextFile tf = new TextFile(f, TextFile.R);
				tf.readLine();
				tf.readLine();
				String[] elems = tf.readLineElems(Strings.whitespace);
				while (elems != null) {
					Chromosome chrom = Chromosome.parseChr(elems[0]);
					Feature feat = new Feature();
					feat.setChromosome(chrom);
					feat.setStart(Integer.parseInt(elems[1]));
					feat.setStop(feat.getStart());
					feat.setName(elems[2]);


					String[] info = elems[elems.length - 1].split("=");
					assocValues.put(feat, Double.parseDouble(info[1]));
					elems = tf.readLineElems(Strings.whitespace);
				}
				tf.close();
			}
		}
		return assocValues;
	}

	private void matchAssoc(String assocdir1, String assocdir2) throws IOException {

		HashMap<Feature, Double> vals1 = loadAssoc(assocdir1);
		HashMap<Feature, Double> vals2 = loadAssoc(assocdir2);

		ArrayList<Double> x = new ArrayList<Double>();
		System.out.println(x.size());
		ArrayList<Double> y = new ArrayList<Double>();
		System.out.println(y.size());
		Set<Feature> k = vals1.keySet();
		for (Feature key : k) {
			Double d = vals2.get(key);
			if (d != null) {
				x.add(vals1.get(key));
				y.add(d);
			}
		}

		System.out.println(x.size() + " associations shared");

		double[] xarr = Primitives.toPrimitiveArr(x.toArray(new Double[0]));
		double[] yarr = Primitives.toPrimitiveArr(y.toArray(new Double[0]));

		double corr = JSci.maths.ArrayMath.correlation(xarr, yarr);
		System.out.println(corr);

	}

	private void matchRSq(String assocdir1, String assocdir2, String out) throws IOException {

		HashMap<Feature, Double> vals1 = loadRsq(assocdir1);
		HashMap<Feature, Double> vals2 = loadRsq(assocdir2);

		ArrayList<Double> x = new ArrayList<Double>();
		System.out.println(x.size());
		ArrayList<Double> y = new ArrayList<Double>();
		System.out.println(y.size());
		Set<Feature> k = vals1.keySet();
		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln("\t" + assocdir1 + "\t" + assocdir2);
		for (Feature key : k) {
			Double d = vals2.get(key);
			if (d != null) {
				x.add(vals1.get(key));
				y.add(d);
				outf.writeln(key.getChromosome() + "_" + key.getStart() + "_" + key.getName() + "\t" + vals1.get(key) + "\t" + d);
			}
		}
		outf.close();
		System.out.println(x.size() + " associations shared");

		double[] xarr = Primitives.toPrimitiveArr(x.toArray(new Double[0]));
		double[] yarr = Primitives.toPrimitiveArr(y.toArray(new Double[0]));


		double corr = JSci.maths.ArrayMath.correlation(xarr, yarr);
		System.out.println(corr);

	}

	private void filterAssoc(String assoc, String filter, String assocOut) throws IOException {
		for (int chr = 1; chr < 23; chr++) {

			if (chr != 3 && chr != 8) {
				System.out.println(chr);
				HashSet<Feature> filterset = new HashSet<Feature>();
				String fin = filter + "matched-mergelog-Chr" + chr + ".txt";
				if (!Gpio.exists(fin)) {
					fin = filter + "mergelog-Chr" + chr + ".txt";
				}
				TextFile tf = new TextFile(fin, TextFile.R);
				tf.readLine();
				String[] elems = tf.readLineElems(TextFile.tab);

				while (elems != null) {
					Chromosome chrom = Chromosome.parseChr(elems[0]);
					Integer pos = Integer.parseInt(elems[1]);
					Feature f = new Feature();
					f.setChromosome(chrom);
					f.setStart(pos);
					f.setStop(pos);
					filterset.add(f);
					elems = tf.readLineElems(TextFile.tab);
				}

				tf.close();

				TextFile associn = new TextFile(assoc + "Chr" + chr + "-assoc.txt", TextFile.R);
				TextFile assocout = new TextFile(assocOut + "Chr" + chr + "-assoc.txt", TextFile.W);
				assocout.writeln(associn.readLine());
				String[] aselems = associn.readLineElems(TextFile.tab);
				while (aselems != null) {
					Feature f = new Feature();
					Chromosome chrom = Chromosome.parseChr(aselems[1]);
					Integer pos = Integer.parseInt(aselems[2]);

					f.setChromosome(chrom);
					f.setStart(pos);
					f.setStop(pos);

					if (!filterset.contains(f)) {
						assocout.writeln(Strings.concat(aselems, Strings.tab));
					}
					aselems = associn.readLineElems(TextFile.tab);
				}

				assocout.close();
				associn.close();
			}
		}
	}

	private void mergeImputedBatches(String prefix, String outfilename, int nrbatches, int chr) throws IOException {
		VCFFunctions f = new VCFFunctions();
		f.mergeImputationBatches(prefix, outfilename, nrbatches, chr);
	}

	public void imputeImmunoChipBatches() throws IOException {

		VCFFunctions f = new VCFFunctions();

		String beagle = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/beagle.r1399.jar";
		String vcfsort = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/vcftools/vcftools_0.1.12b/bin/vcf-sort";
		String jobdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-30-AllChr-EurOnly/15-ImmunoChipSplitJobs/jobs/";
		Gpio.createDir(jobdir);
		System.out.println(jobdir);
		for (int reference = 0; reference < 4; reference++) {
			String refStr = "kg";

			if (reference == 1) {
				refStr = "kgseq";
			} else if (reference == 2) {
				refStr = "seq";
			} else if (reference == 3) {
				refStr = "seqvaronly";
			}
			for (int dataset = 0; dataset < 3; dataset++) {

				String datasetStr = "RA";
				String famfile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/RA.fam";
				int nrbatches = 27;
				if (dataset == 1) {
					datasetStr = "T1D";
					famfile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/T1D.fam";
					nrbatches = 26;
				}
				for (int chromosome = 1; chromosome < 23; chromosome++) {

					String referenceDataset = "";

					if (reference == 0) {
						referenceDataset = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/9-ImputationPanels/ref-matched-sorted-Chr" + chromosome + ".vcf.gz";
					} else if (reference == 1) {
						referenceDataset = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/9-ImputationPanels/ref-test-merged-sorted-Chr" + chromosome + ".vcf.gz";
					} else if (reference == 2) {
						referenceDataset = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/9-ImputationPanels/test-matched-sorted-Chr" + chromosome + ".vcf.gz";
					} else if (reference == 3) {
						referenceDataset = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/4-ICAndSeqMerged/merged-sorted-phased-Chr" + chromosome + ".vcf.gz";
					}


// test-match-sorted-Chr1-batch-0.vcf.gz
					for (int batch = 0; batch < nrbatches; batch++) {
						String input = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/" + datasetStr + "/matched/" + refStr + "/split/Chr" + chromosome + "/test-match-sorted-Chr" + chromosome + "-batch-" + batch + ".vcf.gz";


						String imputeout = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/15-EverythingImputedAgain/" + datasetStr + "/" + refStr + "/Chr" + chromosome + "/";
						String imputedFileOut = imputeout + "imputed-batch-" + batch;
						String jobout = jobdir + datasetStr + "-" + refStr + "-batch-" + batch + "-Chr" + chromosome + ".sh";
						int threads = 16;
						String xmx = "-Xmx20g";

						String command3 = "mkdir -p " + imputeout + "\n"
								+ "nice -n20 java " + xmx + " -jar " + beagle + " \\\n"
								+ "\tref=" + referenceDataset + " \\\n"
								+ "\tped=" + famfile + " \\\n"
								+ "\tgt=" + input + " \\\n"
								+ "\tnthreads=" + threads + " \\\n"
								+ "\tout=" + imputedFileOut + "-Chr" + chromosome + " \\\n"
								+ "\t&> " + imputedFileOut + "-Chr" + chromosome + ".log.txt\n"
								+ "zcat " + imputedFileOut + "-Chr" + chromosome + ".vcf.gz | \\\n"
								+ "\t" + vcfsort + " > "
								+ "\t\t" + imputedFileOut + "-sorted-Chr" + chromosome + ".vcf\n"
								+ "rm " + imputedFileOut + "-sorted-Chr" + chromosome + ".vcf.gz\n"
								+ "gzip " + imputedFileOut + "-sorted-Chr" + chromosome + ".vcf";
						makeJob(command3, jobout + "-impute.sh");

						if (chromosome == 1 && batch == 0) {
							jobout = jobdir + "testscript-" + datasetStr + "-" + refStr + ".sh";
							String command4 = "mkdir -p " + imputeout + "\n"
									+ "nice -n20 java " + xmx + " -jar " + beagle + " \\\n"
									+ "\tref=" + referenceDataset + " \\\n"
									+ "\tped=" + famfile + " \\\n"
									+ "\tgt=" + input + " \\\n"
									+ "\tnthreads=" + threads + " \\\n"
									+ "\tout=" + imputedFileOut + "-Chr" + chromosome + " \n";
							makeJob(command4, jobout + "-impute.sh");
						}
					}


				}
			}
		}

	}

	public void splitVCFs() throws IOException {

		VCFFunctions f = new VCFFunctions();

		for (int reference = 3; reference < 4; reference++) {
			String refStr = "kg";
			if (reference == 1) {
				refStr = "kgseq";
			} else if (reference == 2) {
				refStr = "seq";
			} else if (reference == 3) {
				refStr = "seqvaronly";
			}
			for (int dataset = 0; dataset < 2; dataset++) {

				String datasetStr = "RA";
				String famfile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/RA.fam";
				if (dataset == 1) {
					datasetStr = "T1D";
					famfile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/T1D.fam";
				}

				for (int chromosome = 1; chromosome < 23; chromosome++) {

					String input = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/" + datasetStr + "/matched/" + refStr + "/test-matched-sorted-Chr" + chromosome + ".vcf.gz";
					if (reference == 3) {
						input = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/" + datasetStr + "/matched/" + refStr + "/matched-test-matched-sorted-Chr" + chromosome + ".vcf.gz";
					}
					if (Gpio.exists(input)) {
						String outputprefix = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/" + datasetStr + "/matched/" + refStr + "/split/Chr" + chromosome + "/";
						Gpio.createDir(outputprefix);
						outputprefix += "test-match-sorted-Chr" + chromosome;
						f.splitVCFOverRandomBatches(input, famfile, outputprefix, 1000);
					} else {
						System.out.println("cannot find: " + input);
					}
				}

			}
		}
	}


	private void immunoChipPhaseJobsSequencedVariantsOnly() throws IOException {

		String beagle = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/beagle.r1399.jar";
		String vcfsort = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/vcftools/vcftools_0.1.12b/bin/vcf-sort";
		String jobdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-30-AllChr-EurOnly/14-ImmunoChipSequencedVariantsOnly/jobs/";
		System.out.println(jobdir);

		for (int chromosome = 1; chromosome < 23; chromosome++) {
			for (int dataset = 0; dataset < 3; dataset++) {

				String xmx = "-Xmx20g";
				if (chromosome == 2) {
					xmx = "-Xmx40g";
				}
				String reference = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/4-ICAndSeqMerged/merged-sorted-phased-Chr" + chromosome + ".vcf.gz";
				String immunochipinput = "";
				String famfile = "";
				String datasetName = "";

				String jobout = "";
				if (dataset == 0) {
					datasetName = "T1D";
					jobout = jobdir + "T1D";
					famfile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/T1D.fam";
				} else {
					jobout = jobdir + "RA";
					datasetName = "RA";
					famfile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/RA.fam";
				}

				immunochipinput = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/14-ImmunoChipSequencedVariantsOnly/Matched/" + datasetName + "/matched-test-matched-sorted-Chr" + chromosome + ".vcf.gz";
				jobout += "-" + chromosome;

				int threads = 64;
				String phasedout = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/14-ImmunoChipSequencedVariantsOnly/" + datasetName + "/phased";
				String outdir = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/seqVariantsOnly/";
				String imputeout = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/14-ImmunoChipSequencedVariantsOnly/Imputed/" + datasetName + "/";

				String command2 = "mkdir -p " + phasedout + " \n"

						+ "nice -n20 java " + xmx + " -jar " + beagle + " \\\n"
						+ "\tped=" + famfile + " \\\n"
						+ "\tgt=" + immunochipinput + " \\\n"
						+ "\tnthreads=" + threads + " \\\n"
						+ "\tout=" + phasedout + "-Chr" + chromosome + " \\\n"
						+ "\t&> " + phasedout + ".log.txt\n"

						+ "zcat " + phasedout + "-Chr" + chromosome + ".vcf.gz | \\\n"
						+ "\t" + vcfsort + " > \\\n"
						+ "\t\t" + phasedout + "-sorted-Chr" + chromosome + ".vcf\n"
						+ "rm " + phasedout + "-sorted-Chr" + chromosome + ".vcf.gz\n"
						+ "gzip " + phasedout + "-sorted-Chr" + chromosome + ".vcf";
				makeJob(command2, jobout + "-phase.sh");

				String command3 = "mkdir -p " + imputeout + "\n"
						+ "nice -n20 java " + xmx + " -jar " + beagle + " \\\n"
						+ "\tref=" + reference + " \\\n"
						+ "\tped=" + famfile + " \\\n"
						+ "\tgt=" + phasedout + "-sorted-Chr" + chromosome + ".vcf.gz \\\n"
						+ "\tnthreads=" + threads + " \\\n"
						+ "\tout=" + imputeout + "-Chr" + chromosome + " \\\n"
						+ "\t&> " + imputeout + "-Chr" + chromosome + ".log.txt\n"
						+ "zcat " + imputeout + "-Chr" + chromosome + ".vcf.gz | \\\n"
						+ "\t" + vcfsort + " > "
						+ "\t\t" + imputeout + "-sorted-Chr" + chromosome + ".vcf\n"
						+ "rm " + imputeout + "-sorted-Chr" + chromosome + ".vcf.gz\n"
						+ "gzip " + imputeout + "-sorted-Chr" + chromosome + ".vcf";
				makeJob(command3, jobout + "-impute.sh");


			}
		}
	}

	private void immunoChipAssociationPlots() throws IOException {
		try {
			String sequencedRegionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/2015-03-30-allRegionsWith50PercentOfSamplesAbove20x.bed";
			String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";
//			String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/interesting.bed";
			String gtf = "/Data/Annotation/UCSC/genes.gtf";
			/*

			 */

			//
			double mafthreshold = 0.005;
			boolean logtransform = false;
			boolean targetregionsonly = false;
			boolean pseudo = false;
			if (targetregionsonly) {
				regionFile = sequencedRegionFile;
			}

			String[] referenceSets = new String[]{"before", "kg", "seqvaronly"};
			for (int dataset = 0; dataset < 2; dataset++) {
				String datasetStr = "RA";
				if (dataset == 1) {
					datasetStr = "T1D";
				}
				for (int chr = 1; chr < 23; chr++) {
					Chromosome chrin = Chromosome.parseChr("" + chr);
					String[] datasetFiles = new String[referenceSets.length];
					String output = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/Assoc/PlotsNovelVariantsPerReference/" + datasetStr + "/";
					Gpio.createDir(output);

					for (int reference = 0; reference < referenceSets.length; reference++) {

						String refStr = referenceSets[reference];
						if (pseudo && datasetStr.equals("T1D")) {
							refStr += "-pseudo";
						}

						if (dataset == 0 && refStr.equals("ICStudy")) {
							datasetFiles[datasetFiles.length - 1] = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/RA-Eyre/hg19_gwas_ic_ra_eyre_4_18_0.tab";

						} else if (refStr.equals("ICStudy")) {
							datasetFiles[datasetFiles.length - 1] = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/T1D-Onengut/hg19_gwas_ic_t1d_onengut_cc_4_18_1.tab";
						} else {
							String afterImputation = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/Assoc-OnlyImputed/" + datasetStr + "/" + refStr + "/" + chrin.toString() + "-assoc.txt";
							datasetFiles[reference] = afterImputation;
						}
					}


					TextFile tf2 = new TextFile(regionFile, TextFile.R);
					String[] elems = tf2.readLineElems(TextFile.tab);

					while (elems != null) {
						Feature region = new Feature();
						if (elems.length > 2) {
							if (targetregionsonly) {
								region.setChromosome(Chromosome.parseChr(elems[0]));
								region.setStart(Integer.parseInt(elems[1]) - 10000);
								region.setStop(Integer.parseInt(elems[2]) + 10000);
							} else {
								region.setChromosome(Chromosome.parseChr(elems[0]));
								region.setStart(Integer.parseInt(elems[1]));
								region.setStop(Integer.parseInt(elems[2]));
							}
							if (region.getChromosome().equals(chrin)) {

								String plotout = output + region.toString() + ".pdf";
								VariantPlot plot = new VariantPlot(plotout, 1400, 1500);
								plot.setMargin(200);

								plot.plotAssocPvalue(gtf,
										datasetFiles,
										referenceSets,
										sequencedRegionFile,
										region,
										mafthreshold,
										logtransform);
							}
						}
						elems = tf2.readLineElems(TextFile.tab);
					}

					tf2.close();
				}
			}


			// String[] variantFiles, String[] variantFileNames, String sequencedRegionFile, String regionFile


//			plot.close();

		} catch (DocumentException e) {
			e.printStackTrace();
		}

	}

	public void mergeImmunoChipVCF() throws IOException {

		String data1 = "/Data/ImmunoChip/T1D/AllRegions/eur/genotypes-filtered-sorted-";
		String data2 = "/Data/ImmunoChip/T1D/AllRegions/uk/genotypes-filtered-sorted-";

		String outdir = "/Data/ImmunoChip/T1D/AllRegions/merged/genotypes-filtered-sorted-";

		VCFFunctions v = new VCFFunctions();
		for (int i = 1; i < 23; i++) {

			String ref = data1 + "Chr" + i + ".vcf.gz";
			String test = data2 + "Chr" + i + ".vcf.gz";
			if (Gpio.exists(ref)) {
				String refout = outdir + "-ref.vcf.gz";
				String testout = outdir + "-ref.vcf.gz";
				String mergedout = outdir + "Chr" + i + ".vcf.gz";
				String log = outdir + "Chr" + i + ".log.gz";
				v.mergeAndIntersectVCFVariants(ref, test, refout, testout, mergedout, "/", log, false);
			}

		}

	}

	public void testVCFLocal(int nrThreads) throws IOException {
		String covariatefile = "";
		String outputdir = "/Data/tmp/pseudocontroltest/";
		String diseasestatus = "";

		HashSet<String> covariatestoinclude = null;
		String samplestoexclude = null;
		boolean filterimputationquality = false;
		double imputationqualitythreshold = 0.8;
		int mafthreshold = 5;

		System.out.println("Opening threadpool for " + nrThreads + " threads.");
		ExecutorService threadPool = Executors.newFixedThreadPool(nrThreads);
		CompletionService<Boolean> pool = new ExecutorCompletionService<Boolean>(threadPool);

		HashSet<String> snpLimit = null;

		int submit = 0;
		int start = 1;
		int stop = 23;
		String famfile = null;
		for (int d = 0; d < 1; d++) {

//			if (d == 0) {
//				start = 1;
//				stop = 23;
//			} else {
//				start = 1;
//				stop = 23;
//			}

			for (int chr = start; chr < stop; chr++) {
				String vcf = "";
				String out = "";
				for (int reference = 2; reference < 3; reference++) {
					String ref = "beforeImputation";


					if (d == 0) {
						covariatefile = "/Data/ImmunoChip/T1D/binary/covarmerged.txtmergedCovariates.txt";
						//vcf = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationOutput/T1D/merged-filtered-Chr" + chr + ".vcf.gz";

						vcf = "/Data/ImmunoChip/T1D/merged/merged-Chr" + chr + ".vcf.gz";
						diseasestatus = "/Data/ImmunoChip/T1D/binary/covarmerged.txtmergeddisease.txt";
						out = outputdir + "T1D/" + ref + "/";
						famfile = "/Data/ImmunoChip/T1D/binary/covarmerged.txtmergedfam.fam";
					} else {
						covariatefile = "/Data/ImmunoChip/RA/covar/covarmerged.txtmergedCovariates.txt";
						vcf = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationOutput/RA/merged-filtered-Chr" + chr + ".vcf.gz";
						diseasestatus = "/Data/ImmunoChip/RA/covar/covarmerged.txtmergeddisease.txt";
						out = outputdir + "RA/" + ref + "/";
						famfile = "/Data/ImmunoChip/RA/covar/covarmerged.txtmergedfam.fam";
					}
					Gpio.createDir(out);


					String pseudoOut = "";
					if (d == 0) {
						pseudoOut = out + "/pseudo/";
						Gpio.createDir(pseudoOut);
						pseudoOut += "Chr" + chr + "-";
					}
					out += "Chr" + chr + "-";

					if (chr != 8 && chr != 3) {
						if (Gpio.exists(vcf)) {
							if (d == 0) {
								LogitTestR testObj = new LogitTestR(vcf,
										out,
										diseasestatus,
										covariatefile,
										covariatestoinclude,
										snpLimit,
										samplestoexclude,
										null,
										filterimputationquality,
										imputationqualitythreshold,
										mafthreshold, submit);

								pool.submit(testObj);
								submit++;
								testObj = new LogitTestR(vcf,
										pseudoOut,
										diseasestatus,
										covariatefile,
										covariatestoinclude,
										snpLimit,
										samplestoexclude,
										famfile,
										filterimputationquality,
										imputationqualitythreshold,
										mafthreshold, submit);

								pool.submit(testObj);
								submit++;
							} else {
								LogitTestR testObj = new LogitTestR(vcf,
										out,
										diseasestatus,
										covariatefile,
										covariatestoinclude,
										snpLimit,
										samplestoexclude,
										famfile,
										filterimputationquality,
										imputationqualitythreshold,
										mafthreshold, submit);

								pool.submit(testObj);
								submit++;
							}

						} else {
							System.err.println("file does not exist: " + vcf);
						}
					}
				}

			}

		}

		int returned = 0;
		while (returned < submit) {

			try {


				Boolean result = pool.take().get();

				if (result) {
					returned++;
				} else {

					System.exit(-1);
				}


				System.out.println(returned + " / " + submit + " returned");
			} catch (Exception e) {
				e.printStackTrace();
			}

		}

		threadPool.shutdown();

	}

	public void testVCF(int d, int reference, int nrThreads, boolean pseudo) throws IOException {

//		String vcf = "/Data/tmp/test-matched-Chr10.vcf.gz";
//		String covariatefile = "/Data/ImmunoChip/RA/covar/covarmerged.txtmergedCovariates.txt";
//		String outputdir = "/Data/tmp/";
//		String diseasestatus = "/Data/ImmunoChip/RA/covar/covarmerged.txtmergeddisease.txt";
//
//		HashSet<String> covariatestoinclude = null;
//		String samplestoexclude = null;
//		boolean filterimputationquality = false;
//		double imputationqualitythreshold = 0.8;
//		double mafthreshold = 0.005;
//
//		test.logitTest(vcf, outputdir, diseasestatus, covariatefile, covariatestoinclude, samplestoexclude, filterimputationquality, imputationqualitythreshold, mafthreshold);
//"/Data/ImmunoChip/T1D/binary/plinkvcf/plink.vcf";

		String covariatefile = "";
		String outputdir = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/15-EverythingImputedAgain/2015-06-22-Assoc/";
		String diseasestatus = "";

		HashSet<String> covariatestoinclude = null;
		String samplestoexclude = null;
		boolean filterimputationquality = true;
		double imputationqualitythreshold = 0.8;
		int mafthreshold = 5;

		System.out.println("Opening threadpool for " + nrThreads + " threads.");
		ExecutorService threadPool = Executors.newFixedThreadPool(nrThreads);
		CompletionService<Boolean> pool = new ExecutorCompletionService<Boolean>(threadPool);


		int submit = 0;
		int start = 1;
		int stop = 23;
		String famfile = null;
		String[] refs = new String[]{"before", "kg", "seq", "kgseq", "seqvaronly"};

		int endD = 2;
		if (pseudo) {
			endD = 1;
		}
//		for (int d = 0; d < endD; d++) {

			if (d == 0) {
				start = 1;
				stop = 23;
			} else {
				start = 1;
				stop = 23;
			}

			for (int chr = start; chr < stop; chr++) {
				String vcf = "";
				String out = "";
				//for (int reference = 0; reference < refs.length; reference++) {
				String ref = refs[reference];

				if (ref.equals("before")) {
					filterimputationquality = false;
				} else {
					filterimputationquality = true;
				}

				if (d == 0) {
					covariatefile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/T1D-covarmerged.txtmergedCovariates.txt";
					vcf = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/15-EverythingImputedAgain/T1D/" + ref + "/merged-filtered.Chr" + chr + ".vcf.gz";
					if (ref.equals("before")) {
						vcf = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip2/T1D/merged/merged-Chr" + chr + ".vcf.gz";
					}
					diseasestatus = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/T1D-covarmerged.txtmergeddisease.txt";
					if (pseudo) {
						famfile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/T1D.fam";
						out = outputdir + "T1D/" + ref + "-pseudo/";
					} else {
						out = outputdir + "T1D/" + ref + "/";
					}

				} else {
					covariatefile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/RA-covarmerged.txtmergedCovariates.txt";
					vcf = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/15-EverythingImputedAgain/RA/" + ref + "/merged-filtered.Chr" + chr + ".vcf.gz";
					diseasestatus = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/RA-covarmerged.txtmergeddisease.txt";
					if (ref.equals("before")) {
						vcf = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip2/RA/merged/merged-Chr" + chr + ".vcf.gz";
					}
					out = outputdir + "RA/" + ref + "/";
				}
				Gpio.createDir(out);
				out += "Chr" + chr + "-";

				HashSet<String> snpLimit = null;

				if (chr != 8 && chr != 3) {
					if (Gpio.exists(vcf)) {
						LogitTestR testObj = new LogitTestR(vcf,
								out,
								diseasestatus,
								covariatefile,
								covariatestoinclude,
								snpLimit,
								samplestoexclude,
								famfile,
								filterimputationquality,
								imputationqualitythreshold,
								mafthreshold, submit);
						pool.submit(testObj);
						submit++;
					}
				}
//			}

		}

//		}

		int returned = 0;
		while (returned < submit) {

			try {
				Boolean result = pool.take().get();
				if (result) {
					returned++;
				} else {

					System.exit(-1);
				}

				System.out.println(returned + " / " + submit + " returned");
			} catch (Exception e) {
				e.printStackTrace();
			}

		}

		threadPool.shutdown();


	}

	public void matchfiles(String cov, String fam) throws IOException {
		TextFile f1 = new TextFile(cov, TextFile.R);
		HashSet<String> c = new HashSet<String>();
		String[] elems = f1.readLineElems(Strings.whitespace);
		while (elems != null) {
			String sample = elems[0] + "_" + elems[1];
			c.add(sample);
			elems = f1.readLineElems(Strings.whitespace);
		}
		f1.close();

		System.out.println(c.size() + " covariates");
		TextFile f2 = new TextFile(fam, TextFile.R);

		elems = f2.readLineElems(Strings.whitespace);
		int match = 0;
		HashSet<String> c2 = new HashSet<String>();
		while (elems != null) {
			String sample = elems[0] + "_" + elems[1];
			if (c.contains(sample)) {
				match++;
			}
			c2.add(sample);
			elems = f2.readLineElems(Strings.whitespace);
		}
		System.out.println(c2.size() + " fam");
		System.out.println(match);
		f2.close();
	}

	public void immunoChipCovariates() throws IOException {
		PlayGround g = new PlayGround();
		String[] famfiles = new String[]{
				"/Data/ImmunoChip/RA/binary/iChip_RACI_PhaseII_ES_QCgp.fam",
				"/Data/ImmunoChip/RA/binary/iChip_RACI_PhaseII_NL_QCgp.fam",
				"/Data/ImmunoChip/RA/binary/iChip_RACI_PhaseII_SE-E_QCgp.fam",
				"/Data/ImmunoChip/RA/binary/iChip_RACI_PhaseII_SE-U_QCgp.fam",
				"/Data/ImmunoChip/RA/binary/iChip_RACI_PhaseII_UK_QCgp.fam",
				"/Data/ImmunoChip/RA/binary/iChip_RACI_PhaseII_US_QCgp.fam"};

		String[] covariatefiles = new String[]{"/Data/ImmunoChip/RA/covar/ES-pca.covar",
				"/Data/ImmunoChip/RA/covar/NL-pca.covar",
				"/Data/ImmunoChip/RA/covar/SE-E-pca.covar",
				"/Data/ImmunoChip/RA/covar/SE-U-pca.covar",
				"/Data/ImmunoChip/RA/covar/UK-pca.covar",
				"/Data/ImmunoChip/RA/covar/US-pca.covar"
		};

		String[] cohortnames = new String[]{"ES", "NL", "SE-E", "SE-U", "UK", "US"};
		String output = "/Data/ImmunoChip/RA/covar/covarmerged.txt";
		g.createCovariateMatrixAndMergeFamFiles(famfiles, covariatefiles, cohortnames, output);


		famfiles = new String[]{"/Data/ImmunoChip/T1D/binary/UK.fam",
				"/Data/ImmunoChip/T1D/binary/eur.fam",
		};
		covariatefiles = new String[]{"/Data/ImmunoChip/T1D/binary/UK.cov",
				"/Data/ImmunoChip/T1D/binary/eur.cov"};


		cohortnames = new String[]{"EUR", "T1D"};
		output = "/Data/ImmunoChip/T1D/binary/covarmerged.txt";
		g.createCovariateMatrixAndMergeFamFiles(famfiles, covariatefiles, cohortnames, output);
	}

	private void immunoChipImputationjobs() throws IOException {
		String jobdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-30-AllChr-EurOnly/12-ImmunoChipImputeJobs/";
		Gpio.createDir(jobdir);
		System.out.println("Outdir: " + jobdir);
		String vcfsort = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/vcftools/vcftools_0.1.12b/bin/vcf-sort";


		String beagle = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/beagle.r1399.jar";
		String merger = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip2/merger.jar";

		String kgPrefix = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/9-ImputationPanels/ref-matched-sorted";
		String kgseqPrefix = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/9-ImputationPanels/ref-test-merged-sorted";
		String seqPrefix = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/9-ImputationPanels/test-matched-sorted";

		String outdir = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/";


//		String immunochipmergeddir = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip2/";

		String immunochipdir = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/";

		for (int dataset = 0; dataset < 2; dataset++) {
			for (int reference = 0; reference < 3; reference++) {
				for (int chromosome = 1; chromosome < 23; chromosome++) {

					String jobout = "";
					String datasetName = "";
					String refprefix = "";
					String imputeout = "";
					String famfile = "";

					if (dataset == 0) {
						datasetName = "T1D";
						jobout = jobdir + "T1D";
						famfile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/T1D.fam";
					} else {
						jobout = jobdir + "RA";
						datasetName = "RA";
						famfile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/RA.fam";
					}

					String immunochipinput = ""; //immunochipmergeddir + datasetName + "/merged/merged-Chr" + chromosome + ".vcf.gz";

					if (reference == 0) {
						refprefix = kgPrefix;
						immunochipinput = immunochipdir + datasetName + "/matched/kg/test-matched-sorted-Chr" + chromosome + ".vcf.gz";
						imputeout = outdir + datasetName + "/imputed/kg/test-matched-imputed";
						jobout += "-kg-" + chromosome + ".sh";
					} else if (reference == 1) {
						refprefix = seqPrefix;
						immunochipinput = immunochipdir + datasetName + "/matched/seq/test-matched-sorted-Chr" + chromosome + ".vcf.gz";
						imputeout = outdir + datasetName + "/imputed/seq/test-matched-imputed";
						jobout += "-seq-" + chromosome + ".sh";
					} else {
						refprefix = kgseqPrefix;
						immunochipinput = immunochipdir + datasetName + "/matched/kgseq/test-matched-sorted-Chr" + chromosome + ".vcf.gz";
						imputeout = outdir + datasetName + "/imputed/kgseq/test-matched-imputed";
						jobout += "-kgseq-" + chromosome + ".sh";
					}

					String xmx = "-Xmx40g";
					int threads = 16;
					if (chromosome == 2) {
						xmx = "-Xmx35g";
					}

					String phasedout = immunochipinput + "-phased";

					String command2 = "mkdir -p " + outdir + datasetName + "/imputed/seq/\n" +
							"mkdir -p " + outdir + datasetName + "/imputed/kgseq/\n" +
							"mkdir -p " + outdir + datasetName + "/imputed/kg/\n"

							+ "nice -n20 java " + xmx + " -jar " + beagle + " "
							+ "ped=" + famfile + " "
							+ "gt=" + immunochipinput + " "
							+ "nthreads=" + threads + " "
							+ "out=" + phasedout + "-Chr" + chromosome + " "
							+ "&> " + phasedout + ".log.txt\n"

							+ "zcat " + phasedout + "-Chr" + chromosome + ".vcf.gz | "
							+ " " + vcfsort + " > "
							+ " " + phasedout + "-sorted-Chr" + chromosome + ".vcf\n"
							+ "rm " + phasedout + "-sorted-Chr" + chromosome + ".vcf.gz\n"
							+ "gzip " + phasedout + "-sorted-Chr" + chromosome + ".vcf";
					makeJob(command2, jobout + "-phase.sh");

					String command3 = "nice -n20 java " + xmx + " -jar " + beagle + " "
							+ "ref=" + refprefix + "-Chr" + chromosome + ".vcf.gz "
							+ "usephase=true burnin-its=0 phase-its=0 "
							+ "gt=" + phasedout + "-sorted-Chr" + chromosome + ".vcf.gz "
							+ "nthreads=" + threads + " "
							+ "out=" + imputeout + "-Chr" + chromosome + " "
							+ "&> " + imputeout + "-Chr" + chromosome + ".log.txt\n"
							+ "zcat " + imputeout + "-Chr" + chromosome + ".vcf.gz | "
							+ " " + vcfsort + " > "
							+ " " + imputeout + "-sorted-Chr" + chromosome + ".vcf\n"
							+ "rm " + imputeout + "-sorted-Chr" + chromosome + ".vcf.gz\n"
							+ "gzip " + imputeout + "-sorted-Chr" + chromosome + ".vcf";
					makeJob(command3, jobout);

				}
			}
		}

	}

	public void filterEuropean(String vcfIn, String vcfOut, String sampleList) throws IOException {
		VCFFunctions v = new VCFFunctions();
		v.filterSamples(vcfIn, vcfOut, sampleList);
	}


	public void countVariantsWithR2(String file1, String limitToTheseVars, String outf, int distsize) throws IOException {
		int[] ar2dist = new int[distsize + 1];
		int[] dr2dist = new int[distsize + 1];
		int[] af2dist = new int[distsize + 1];


		HashSet<String> vars = null;


		if (limitToTheseVars != null) {
			TextFile t = new TextFile(limitToTheseVars, TextFile.R);
			vars = new HashSet<String>();

			String[] elems = t.readLineElems(TextFile.tab);
			while (elems != null) {

				vars.add(elems[0] + "_" + elems[1]);

				elems = t.readLineElems(TextFile.tab);
			}

			t.close();
		}


		if (file1.endsWith(".r2")) {

			TextFile tf = new TextFile(file1, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {

				String rs = elems[0];
				Double r2 = Double.parseDouble(elems[1]);
				if (Double.isNaN(r2)) {
					dr2dist[0]++;
				} else {
					int dr2int = (int) Math.ceil(r2 * distsize);
					dr2dist[dr2int]++;
				}

				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();

		} else {
			VCFGenotypeData data = new VCFGenotypeData(file1);
			TextFile out = new TextFile(outf, TextFile.W);
			while (data.hasNext()) {

				VCFVariant var = data.next();

				String id = var.getId();
				String pos = var.getChr() + "_" + var.getPos();
				if (vars == null || !vars.contains(pos)) {
					HashMap<String, Double> info = var.getInfo();
					Double ar2 = info.get("AR2");
					Double dr2 = info.get("DR2");
					Double af = info.get("AF");


					if (ar2 > 1 || dr2 > 1) {
						System.out.println(id + "\t" + ar2 + "\t" + dr2);
					} else {
						int dr2int = (int) Math.ceil(dr2 * distsize);
						int ar2int = (int) Math.ceil(ar2 * distsize);

						if (af != null) {
							if (af > 0.005) {
								int af2int = (int) Math.ceil(af * distsize);
								out.writeln(af + "\t" + ar2 + "\t" + dr2);

								ar2dist[ar2int]++;
								dr2dist[dr2int]++;
								af2dist[af2int]++;
							}

						} else {
							// System.err.println("AF2 == null for " + pos + "-" + id);
						}


					}
				}
			}
			data.close();
			out.close();
		}


		System.out.println();
		double d = 1d / distsize;
		for (int i = 0; i < ar2dist.length; i++) {
			System.out.println((d * i) + "\t" + ar2dist[i] + "\t" + dr2dist[i] + "\t" + af2dist[i]);
		}
	}

	private void filterMarkersFile(String filterthis, String withthis, String outputhere) throws IOException {
		TextFile tf = new TextFile(withthis, TextFile.R);
		HashSet<String> set = (HashSet<String>) tf.readAsSet(0, Strings.pipe);

		tf.close();

		TextFile tfout = new TextFile(outputhere, TextFile.W);
		TextFile tfin = new TextFile(filterthis, TextFile.R);
		String ln = tfin.readLine();
		while (ln != null) {
			if (!set.contains(ln)) {
				tfout.writeln();
			}
			ln = tfin.readLine();
		}

		tfin.close();
		tfout.close();

	}

	public void compareGenotypesBetweenDatasets(String set1, String set2, String output) throws IOException {

		//
		GenotypeTools t = new GenotypeTools();
		t.compareVCFGenotypesToPedAndMap(set1, set2, output, true);
	}

	public void immunochipMatchJobs() throws IOException {
		String jobdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-30-AllChr-EurOnly/11-ImmunoChipMatchJobs/";
		String vcfSort = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/vcftools/vcftools_0.1.12b/bin/vcf-sort";

		String merger = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip2/merger.jar";

		String kgPrefix = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/9-ImputationPanels/ref-matched-sorted";
		String seqPrefix = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/9-ImputationPanels/ref-test-merged-sorted";
		String kgseqPrefix = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/9-ImputationPanels/test-matched-sorted";

		String outdir = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/";


		String immunochipmergeddir = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip2/";

		for (int dataset = 0; dataset < 2; dataset++) {
			for (int reference = 0; reference < 3; reference++) {
				for (int chromosome = 1; chromosome < 23; chromosome++) {

					String jobout = "";
					String datasetName = "";
					String refprefix = "";
					String matchout = "";

					if (dataset == 0) {
						datasetName = "T1D";
						jobout = jobdir + "T1D";
					} else {
						jobout = jobdir + "RA";
						datasetName = "RA";
					}

					String immunochipinput = immunochipmergeddir + datasetName + "/merged/merged-Chr" + chromosome + ".vcf.gz";

					if (reference == 0) {
						refprefix = kgPrefix;
						matchout = outdir + datasetName + "/matched/kg/";
						jobout += "-kg-" + chromosome + ".sh";
					} else if (reference == 1) {
						refprefix = seqPrefix;
						matchout = outdir + datasetName + "/matched/seq/";
						jobout += "-seq-" + chromosome + ".sh";
					} else {
						refprefix = kgseqPrefix;
						matchout = outdir + datasetName + "/matched/kgseq/";
						jobout += "-kgseq-" + chromosome + ".sh";
					}

					String xmx = "-Xmx10g";
					if (chromosome == 2) {
						xmx = "-Xmx35g";
					}

					String command1 = "mkdir -p " + matchout + "\n" +
							"nice -n20 java " + xmx + " -jar " + merger + " " + chromosome + " \\\n"
							+ vcfSort + " \\\n"
							+ refprefix + "-Chr" + chromosome + ".vcf.gz \\\n"
							+ immunochipinput + " \\\n"
							+ matchout + " false \\\n&> " + matchout + "log-Chr" + chromosome + ".txt\n";

					makeJob(command1, jobout);


				}
			}
		}


	}

	public void panelMergeJobs() throws IOException {
		String jobdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-30-AllChr-EurOnly/9-PanelMergeJobs/";
		Gpio.createDir(jobdir);

		for (int i = 1; i < 23; i++) {

			String xmx = "-Xmx10g";
			if (i == 2) {
				xmx = "-Xmx35g";
			}

			String command1 = "nice -n20 java " + xmx + " -jar /medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/mergerphased.jar " + i + " "
					+ "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/vcftools/vcftools_0.1.12b/bin/vcf-sort "
					+ "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/8-PanelsImputedMerged/ref-merged-unimputedConcatenated-sorted-Chr" + i + ".vcf.gz "
					+ "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/8-PanelsImputedMerged/test-merged-unimputedConcatenated-sorted-Chr" + i + ".vcf.gz "
					+ "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/9-ImputationPanels/ true &> /medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/9-ImputationPanels/log-Chr" + i + ".txt";

			String bashCommandName1 = jobdir + "merge-Chr" + i + ".sh";
			makeJob(command1, bashCommandName1);

		}
	}

	public void panelImputationJobs() throws IOException {

		String jobdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-30-AllChr-EurOnly/7-PanelImputeJobs/";

		String refdirOnServer = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/6-PanelsMatched/";
		String famfileOnServer = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/mergefam.txt";
		String outdirOnServer = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/7-PanelsImputed/";
		String vcfsort = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/vcftools/vcftools_0.1.12b/bin/vcf-sort";

		for (int i = 1; i < 23; i++) {

			String xmx = "-Xmx10g";
			if (i == 2) {
				xmx = "-Xmx35g";
			}

			String command1 = "nice -n20 java " + xmx + " -jar /medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/beagle.r1399.jar "
					+ "ref=" + refdirOnServer + "1kg-matched-Chr" + i + ".vcf.gz "
					+ "ped=" + famfileOnServer + " "
					+ "gt=" + refdirOnServer + "test-matched-sorted-Chr" + i + ".vcf.gz "
					+ "nthreads=32 "
					+ "out=" + outdirOnServer + "test-unphased-refimputed-Chr" + i + " "
					+ "&> " + outdirOnServer + "test-unphased-refimputed-Chr" + i + ".imputelog.txt\n"
					+ "zcat " + outdirOnServer + "test-unphased-refimputed-Chr" + i + ".vcf.gz | "
					+ " " + vcfsort + " > " + outdirOnServer + "test-unphased-refimputed-sorted-Chr" + i + ".vcf\n"
					+ "rm " + outdirOnServer + "test-unphased-refimputed-sorted-Chr" + i + ".vcf.gz\n"
					+ "gzip " + outdirOnServer + "test-unphased-refimputed-sorted-Chr" + i + ".vcf";

			String bashCommandName1 = jobdir + "test-unphased-Chr" + i + ".sh";
			// makeJob(command1, bashCommandName1);

			String command2 = "nice -n20 java " + xmx + " -jar /medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/beagle.r1399.jar "
					+ "ref=" + refdirOnServer + "1kg-matched-Chr" + i + ".vcf.gz "
					+ "usephase=true burnin-its=0 phase-its=0 "
					+ "gt=" + refdirOnServer + "seqpanel-matched-sorted-phased-sorted-Chr" + i + ".vcf.gz "
					+ "nthreads=32 "
					+ "out=" + outdirOnServer + "test-phased-refimputed-Chr" + i + " "
					+ "&> " + outdirOnServer + "test-phased-refimputed-Chr" + i + ".imputelog.txt\n"
					+ "zcat " + outdirOnServer + "test-phased-refimputed-Chr" + i + ".vcf.gz | "
					+ " " + vcfsort + " > " +
					" " + outdirOnServer + "test-phased-refimputed-sorted-Chr" + i + ".vcf\n"
					+ "rm " + outdirOnServer + "test-phased-refimputed-sorted-Chr" + i + ".vcf.gz\n"
					+ "gzip " + outdirOnServer + "test-phased-refimputed-sorted-Chr" + i + ".vcf";

			String bashCommandName2 = jobdir + "test-phased-Chr" + i + ".sh";
			makeJob(command2, bashCommandName2);
			String command3 = "nice -n20 java " + xmx + " -jar /medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/beagle.r1399.jar "
					+ "ref=" + refdirOnServer + "seqpanel-matched-sorted-phased-sorted-Chr" + i + ".vcf.gz "
					+ "usephase=true burnin-its=0 phase-its=0 "
					+ "gt=" + refdirOnServer + "1kg-matched-Chr" + i + ".vcf.gz "
					+ "nthreads=32 "
					+ "out=" + outdirOnServer + "ref-testimputed-Chr" + i + " "
					+ "&> " + outdirOnServer + "ref-testimputed-Chr" + i + ".imputelog.txt\n"
					+ "zcat " + outdirOnServer + "ref-testimputed-Chr" + i + ".vcf.gz | "
					+ " " + vcfsort + " > " +
					" " + outdirOnServer + "ref-testimputed-sorted-Chr" + i + ".vcf\n"
					+ "rm " + outdirOnServer + "ref-testimputed-sorted-Chr" + i + ".vcf.gz\n"
					+ "gzip " + outdirOnServer + "ref-testimputed-sorted-Chr" + i + ".vcf";

			String bashCommandName3 = jobdir + "ref-Chr" + i + ".sh";
			makeJob(command3, bashCommandName3);
		}

		/*


						pb = new ProcessBuilder("java",



						// impute reference panel into 1kg
						pb = new ProcessBuilder("java",
								"-Xmx6g",
								"-jar", beagle,
								"ref=" + matchedPanelsOut + "seqpanel-matched-sorted-phased-sorted-" + chr.getName() + ".vcf.gz",
								"ped=" + mergeFam,
								"gt=" + matchedPanelsOut + "1kg-matched-sorted-" + chr.getName() + ".vcf.gz",
								"nthreads=32",
								"out=" + imputedPanelsOut + "1kg-seqpanelimputed-" + chr.getName()
						);
						t.run(pb);
						t.sortVCF(vcfsort, imputedPanelsOut + "1kg-seqpanelimputed-" + chr.getName() + ".vcf.gz",
								imputedPanelsOut + "1kg-seqpanelimputed-sorted-" + chr.getName() + ".vcf.gz", imputedPanelsOut + "sort.sh");
		 */

	}

	public void panelSortAndConcatJobs() throws IOException {

		// imp:
		// unimp:

		// imp:
		// unimp: -Chr


		String jobdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-30-AllChr-EurOnly/8-PanelMergeJobs/";

		String vcfsort = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/vcftools/vcftools_0.1.12b/bin/vcf-sort";
		String vcfconcat = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/vcftools/vcftools_0.1.12b/bin/vcf-concat";

		String inputdir = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/7-PanelsImputed/";
		String unimputeddir = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/6-PanelsMatched/";
		String outputdir = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/8-PanelsImputedMerged/";

		for (int chr = 1; chr < 23; chr++) {

			String inputprefix = inputdir + "test-phased-refimputed-sorted-rsquarefiltered";
			String unimputedprefix = unimputeddir + "seqpanel-matched-sorted-phased-sorted";
			String outputprefix = outputdir + "test-merged";


			String command = "nice -n20 java -Xmx5g -jar /medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/panelSortAndConcat.jar " + chr + " "
					+ inputprefix + " "
					+ outputprefix + " "
					+ unimputedprefix + " "
					+ vcfsort + " "
					+ vcfconcat;
			String bashout = jobdir + "test-Chr" + chr + ".sh";

			makeJob(command, bashout);

			inputprefix = inputdir + "ref-testimputed-sorted-rsquarefiltered";
			String kgunimputeddir = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/ref/5-1KG-SequencedRegions-HighFreq/";
			unimputedprefix = kgunimputeddir + "filtered";
			outputprefix = outputdir + "ref-merged";


			command = "nice -n20 java -Xmx5g -jar /medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/panelSortAndConcat.jar " + chr + " "
					+ inputprefix + " "
					+ outputprefix + " "
					+ unimputedprefix + " "
					+ vcfsort + " "
					+ vcfconcat;
			bashout = jobdir + "ref-Chr" + chr + ".sh";
			makeJob(command, bashout);

		}


	}

	public void immunochipSortAndConcatJobs() throws IOException {

		// imp:
		// unimp:

		// imp:
		// unimp: -Chr


		String jobdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-30-AllChr-EurOnly/13-ImmunoChipMergeJobs/";
		Gpio.createDir(jobdir);
		String vcfsort = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/vcftools/vcftools_0.1.12b/bin/vcf-sort";
		String vcfconcat = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/vcftools/vcftools_0.1.12b/bin/vcf-concat";
		String immunochipdir = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/";

		String imputeoutdir = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/";
		String outputprefix = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/";
		for (int dataset = 0; dataset < 2; dataset++) {
			for (int reference = 0; reference < 3; reference++) {
				for (int chromosome = 1; chromosome < 23; chromosome++) {

					String jobout = "";
					String datasetName = "";
					String refprefix = "";
					String imputeout = "";
					String famfile = "";

					if (dataset == 0) {
						datasetName = "T1D";
						jobout = jobdir + "T1D";
					} else {
						jobout = jobdir + "RA";
						datasetName = "RA";
					}

					String outputdir = "";
					String immunochipinput = ""; //immunochipmergeddir + datasetName + "/merged/merged-Chr" + chromosome + ".vcf.gz";
					immunochipinput = immunochipdir + datasetName + "/matched/kg/test-matched-sorted-Chr" + chromosome + ".vcf.gz";
					if (reference == 0) {
						refprefix = "kg";
						imputeout = imputeoutdir + datasetName + "/imputed/kg/test-matched-imputed-sorted-Chr" + chromosome + ".vcf.gz";
						outputdir = outputprefix + datasetName + "/imputedmerged/kg/test-imputed-merged";
						jobout += "-kg-" + chromosome + ".sh";
					} else if (reference == 1) {
						refprefix = "seq";
						imputeout = imputeoutdir + datasetName + "/imputed/seq/test-matched-imputed-sorted-Chr" + chromosome + ".vcf.gz";
						outputdir = outputprefix + datasetName + "/imputedmerged/seq/test-imputed-merged";
						jobout += "-seq-" + chromosome + ".sh";
					} else {
						refprefix = "kgseq";
						imputeout = imputeoutdir + datasetName + "/imputed/kgseq/test-matched-imputed-sorted-Chr" + chromosome + ".vcf.gz";
						outputdir = outputprefix + datasetName + "/imputedmerged/kgseq/test-imputed-merged";
						jobout += "-kgseq-" + chromosome + ".sh";
					}
					String outdir = outputprefix + datasetName + "/imputedmerged/" + refprefix + "/";

					String xmx = "-Xmx40g";
					int threads = 6;
					if (chromosome == 2) {
						xmx = "-Xmx35g";
					}

					String command = "mkdir -p " + outdir + "\n" +
							"nice -n20 java -Xmx5g -jar /medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/panelSortAndConcat.jar " + chromosome + " "
							+ imputeout + " "
							+ outputdir + " "
							+ immunochipinput + " "
							+ vcfsort + " "
							+ vcfconcat;
					String bashout = jobout;

					makeJob(command, bashout);

				}
			}
		}


	}


	public void panelSortAndConcat(int chrint, String imputedprefix, String outprefix, String unimputedprefix, String vcfsort, String vcfConcat) throws IOException {
		Chromosome chr = Chromosome.parseChr("" + chrint);
		GenotypeTools t = new GenotypeTools();
		VCFFunctions v = new VCFFunctions();

		if (!Gpio.exists(unimputedprefix + "-" + chr.getName() + ".vcf.gz")) {
			System.out.println(unimputedprefix + "-" + chr.getName() + ".vcf.gz does not exist!");
		} else {
			// filter unimputed phased matched dataset for imputed variants (remove the ones that are imputed)


			String unimputedseq = unimputedprefix + "-" + chr.getName() + ".vcf.gz";
			String imputedseq = imputedprefix + "-" + chr.getName() + ".vcf.gz";
			String unimputedseqout = unimputedprefix + "-imputedVariantsRemoved-" + chr.getName() + ".vcf.gz";

			System.out.println("unimp: " + unimputedseq);
			System.out.println("imp: " + imputedseq);
			System.out.println("out: " + unimputedseqout);
			v.filterVCFVariants(unimputedseq, imputedseq, unimputedseqout); // remove imputed variants from original vcf if there is any overlap

			// concatenate unimputed phased variants with those that are imputed.
			String files = unimputedseqout + " " + imputedseq;
			String mergedout = outprefix + "-unimputedConcatenated-" + chr.getName() + ".vcf.gz";
			String mergedoutsorted = outprefix + "-unimputedConcatenated-sorted-" + chr.getName() + ".vcf.gz";
			String bashfilename = outprefix + "-sortmerge-" + chr.getName() + ".sh";
			System.out.println();
			System.out.println("concat: " + files);
			System.out.println("merged: " + mergedout);
			System.out.println("merged sorted: " + mergedoutsorted);

			t.concatVCF(vcfConcat, vcfsort, files, mergedout, mergedoutsorted, bashfilename); // concatenate and sort

		}


	}

	public void phaseJobs() throws IOException {
		for (int i = 1; i < 23; i++) {

			String command = "nice -n20 java -Xmx10g -jar /medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/beagle.r1399.jar "
					+ "ped=/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/mergefam.txt "
					+ "nthreads=32 "
					+ "gt=/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/reference/6-PanelsMatched/test-matched-sorted-Chr" + i + ".vcf.gz "
					+ "out=/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/reference/6-PanelsMatched/test-matched-sorted-phased-Chr" + i;
			String bashfilename = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/6-PanelsMatchedJobs/phase-chr" + i + ".sh";
			makeJob(command, bashfilename);
		}

	}

	public void makeMergeJobs() throws IOException {
		for (int i = 1; i < 23; i++) {
			String command = "cp /medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/reference/5-1KG-SequencedRegions/filtered-Chr" + i + ".vcf.gz /dev/shm/filtered-Chr" + i + ".vcf.gz\n" +
					"cp /medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/reference/4-ICAndSeqVariantMerged/merged-Chr" + i + ".vcf.gz /dev/shm/merged-Chr" + i + ".vcf.gz\n" +
					"nice -n20 java -Xmx10g -jar /medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/mergeAndIntersectVCF.jar "
					+ i + " /medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/vcftools/vcftools_0.1.12b/bin/vcf-sort"
					+ " /dev/shm/filtered-Chr" + i + ".vcf.gz"
					+ " /dev/shm/merged-Chr" + i + ".vcf.gz"
					+ " /medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/reference/6-PanelsMatched/\n"
					+ "rm /dev/shm/filtered-Chr" + i + ".vcf.gz\n"
					+ "rm /dev/shm/merged-Chr" + i + ".vcf.gz";
			String bashfilename = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/6-PanelsMatchedJobs/chr" + i + ".sh";
			makeJob(command, bashfilename);

		}
	}

	public void makeJob(String command, String bashoutfilename) throws IOException {
		TextFile tf = new TextFile(bashoutfilename, TextFile.W);
		tf.writeln("#!/bin/bash");
		tf.writeln(command);
		tf.close();
	}

	public void makePlot() {
		String[] vcfFiles = new String[]{"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/3-ICRegionFiltered/filtered.vcf",
				"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/2-SeqVariantRegionFilter/filtered.vcf",
				"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/5-1KG-SequencedRegions/filtered.vcf"
		};
		String[] vcfNames = new String[]{"IC",
				"SequencePanel",
				"1000Genomes"};
		String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/cd28region.bed";
		String output = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/regionplot.pdf";
		String sequencedRegions = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/2015-03-30-allRegionsWith50PercentOfSamplesAbove20x.bed";
		try {

//			String ped = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/3-ICRegionFiltered/filtered";
//			String vcf = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/3-ICRegionFiltered/filtered.vcf";
//			VCFFunctions v = new VCFFunctions();
//			v.convertPEDToVCF(ped, vcf);

			//plotRegion(vcfFiles, vcfNames, regionFile, sequencedRegions, output);
		} catch (Exception e) {

		}
	}

	public void regionR2Plots(String regionsFile, String outputdir, String before, String after, String sequencingregionfile, double maf) throws IOException {

		try {


			String gtf = "";
			Gpio.createDir(outputdir);

			TextFile tf = new TextFile(regionsFile, TextFile.R);

			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				if (elems.length >= 3) {
					Feature region = new Feature();
					region.setChromosome(Chromosome.parseChr(elems[0]));
					region.setStart(Integer.parseInt(elems[1]));
					region.setStop(Integer.parseInt(elems[2]));

					int width = region.getStop() - region.getStart();

					if (width > 1000000) {
						width /= 1000;
					}
					System.out.println("region width: " + region.toString() + " -- " + width);

					VariantPlot plot = new VariantPlot(outputdir + region.toString() + "-" + maf + ".pdf", width + 200, 500);
					plot.setMargin(100);
//					plot.plotImputationRSquared(before, after, sequencingregionfile, region, maf, gtf);


				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();


		} catch (Exception e1) {
			e1.printStackTrace();
		}
	}


//	public void compare() {
//		try {
//			String vcf1 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/5-1KG-SequencedRegions/filtered.vcf";
//			String vcf1out = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/test/1kg.vcf";
//			String vcf2 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/4-ICAndSeqVariantMerged/merged.vcf";
//			String vcf2out = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/test/seq.vcf";
//			String log = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/test/log.txt";
//			VCFFunctions v = new VCFFunctions();
//			v.compareAndCorrectVCFVariants(
//					vcf1, vcf1out, vcf2, vcf2out, log, false, true);
//
//			vcf2 = vcf2out;
//			vcf2out = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/test/seq2.vcf";
//			log = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/test/log2.txt";
//
//			v.compareAndCorrectVCFVariants(
//					vcf1, vcf1out, vcf2, vcf2out, log, false, true);
//
//
//		} catch (Exception e) {
//
//		}
//	}

	public void prepareImputationPanels() throws IOException {

		GenotypeTools t = new GenotypeTools();
		VCFFunctions v = new VCFFunctions();
		PedAndMapFunctions p = new PedAndMapFunctions();

		// broad
//		String outputPath = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/reference/";
//
//		// refilter the reference VCF.
//		String startVCF = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/merged-ICIds-MixupsFixed.vcf";
//		String regionFile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/allLoci.bed";
//		String mergedhg19immunochipPed = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/immunochipsequencedsamples/merge";
//		String mergedhg19immunochipMap = mergedhg19immunochipPed + ".map";
//		String mergedhg19immunochipFAM = mergedhg19immunochipPed + ".fam";
//		String plink = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/plink/plink";
//		String beagle = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/beagle.r1399.jar";
//		String merged1kg = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/1kg.phase3.v5-filtered.merged.vcf.gz";
//		String vcfConcat = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/vcftools/vcftools_0.1.12b/bin/vcf-concat";
//		String vcfsort = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/vcftools/vcftools_0.1.12b/bin/vcf-sort";
//		String mergeFam = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/mergefam.txt";


		// # local
		String outputPath = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-30-AllChr-EurOnly/";

		// refilter the reference VCF.
		String startVCF = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-11-MixupsFixed/merged-ICIds-MixupsFixed.vcf";
		String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";
		String mergedhg19immunochipPed = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19variantsinsequencedregions/merge";
		String mergedhg19immunochipMap = mergedhg19immunochipPed + ".map";
		String mergedhg19immunochipFAM = mergedhg19immunochipPed + ".fam";
		String plink = "/Data/Tools/plink-1.07-mac-intel/plink1.9";
		String beagle = "/Data/Tools/beagle/beagle.r1399.jar";
		String merged1kg = "/Data/Ref/BeagleRef/1kg.phase3.v5-filtered.eur.merged.vcf.gz"; // EUR ONLY!
		String vcfConcat = "/Data/Tools/vcftools/bin/vcf-concat";
		String vcfsort = "/Data/Tools/vcftools/bin/vcf-sort";
		String mergeFam = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/mergefam.txt";
		boolean skipimpute = true;
		boolean linux = false;
		boolean skipphasing = false;

		try {
			// filter bad variants
			String filteredVCFOut = outputPath + "1-SeqVariantFilter/";
			Gpio.createDir(filteredVCFOut);
			v.summarizeVCF("/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/6-PanelsMatched/seq-matched.vcf", "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/6-PanelsMatched/seq-matchedsummary.txt");
			v.filterLowFrequencyVariants(startVCF, filteredVCFOut, true, 10, 30, 0.90, 5);
			String mendelianFilteredVCF = filteredVCFOut + "filtered-mendelianerrorsremoved.vcf";

			v.filterMendelianErrors(filteredVCFOut + "filtered.vcf.gz", mergeFam, mendelianFilteredVCF, filteredVCFOut);
			v.summarizeVCF(mendelianFilteredVCF, mendelianFilteredVCF + "-summary.txt");
//			System.exit(-1);

//			mendelianFilteredVCF = filteredVCFOut + "filtered.vcf.gz";
			//v.splitMultipleAllelicVariants(filteredVCFOut + "filtered.vcf", mendelianFilteredVCF);

			// compare against IC genotypes
			//	t.compareVCFGenotypesToPedAndMap(filteredVCFOut + "filtered.vcf.gz", mergedhg19immunochipPed, filteredVCFOut, true);

			// filter for CD28 region
			String filteredVariantOut = outputPath + "2-SeqVariantRegionFilter/";
			Gpio.createDir(filteredVariantOut);

			v.filterVCFForBedRegions(mendelianFilteredVCF, filteredVariantOut + "filtered.vcf.gz", regionFile);
			v.summarizeVCF(filteredVariantOut + "filtered.vcf.gz", filteredVariantOut + "summary.txt");

			// filter immunochip for sequenced regions and samples
			String regionFilteredMapOut = outputPath + "3-ICRegionFiltered/";
			Gpio.createDir(regionFilteredMapOut);
			String variantsToKeep = regionFilteredMapOut + "variantsOverlappingSequencedRegions.txt";
			p.filterMap(mergedhg19immunochipMap, regionFile, variantsToKeep);
			String variantsToExclude = regionFilteredMapOut + "ICvariantsOverlappingSequencedVariants.txt";
			p.filterMapForVCFVariants(mergedhg19immunochipMap, filteredVariantOut + "filtered.vcf.gz", variantsToExclude);


			String regionfilteredpedOut = regionFilteredMapOut + "filtered";
			ProcessBuilder pb = new ProcessBuilder(plink,
					"--extract", variantsToKeep,
					"--file", mergedhg19immunochipPed,
					"--exclude", variantsToExclude,
					"--recode",
					"--out", regionfilteredpedOut,
					"--keep", mergeFam
			);
			t.run(pb);


			// merge with immunochip
			String ICAndSeqVariantMerged = outputPath + "4-ICAndSeqVariantMerged/";
			Gpio.createDir(ICAndSeqVariantMerged);


			v.mergeWithPed(regionfilteredpedOut, filteredVariantOut + "filtered.vcf.gz", ICAndSeqVariantMerged + "merged.vcf.gz");
			v.summarizeVCF(ICAndSeqVariantMerged + "merged.vcf.gz", ICAndSeqVariantMerged + "merged-summary.txt");

//			System.exit(-1);
//
			// filter 1000 genomes for sequenced regions
			String kgSeqRegions = outputPath + "5-1KG-SequencedRegions/";
			Gpio.createDir(kgSeqRegions);
			System.out.println("Filtering VCF for regions: " + merged1kg);
			v.filterVCFForBedRegions(merged1kg, kgSeqRegions + "filtered.vcf.gz", regionFile);

			String kgSeqRegionsHigherfreq = outputPath + "5-1KG-SequencedRegions-HighFreq/";
			Gpio.createDir(kgSeqRegions);
			v.filterLowFrequencyVariants(kgSeqRegions + "filtered.vcf.gz", kgSeqRegionsHigherfreq, true, 0, 0, 0, 1);
			System.out.println("Summarizing: " + kgSeqRegionsHigherfreq + "filtered.vcf.gz");
			v.summarizeVCF(kgSeqRegionsHigherfreq + "filtered.vcf.gz", kgSeqRegionsHigherfreq + "filtered-summary.txt");
			kgSeqRegions = kgSeqRegionsHigherfreq;


			// compare reference and 1kg

			// split per chromosome...
			System.out.println("Splitting ref per chr");
			v.splitPerChromosome(kgSeqRegions + "filtered.vcf.gz", kgSeqRegions + "filtered");
			System.out.println("Splitting test per chr");
			v.splitPerChromosome(ICAndSeqVariantMerged + "merged.vcf.gz", ICAndSeqVariantMerged + "merged");

			for (Chromosome chr : Chromosome.values()) {
				System.out.println(chr.getName() + " imputing..");

				String refVCFChr = kgSeqRegions + "filtered-" + chr.getName() + ".vcf.gz";
				String testVCFChr = ICAndSeqVariantMerged + "merged-" + chr.getName() + ".vcf.gz";
				if (Gpio.exists(refVCFChr) && Gpio.exists(testVCFChr)) {

					String matchedPanelsOut = outputPath + "6-PanelsMatched/";
					Gpio.createDir(matchedPanelsOut);

					v.mergeAndIntersectVCFVariants(
							refVCFChr,
							testVCFChr,
							matchedPanelsOut + "1kg-matched-" + chr.getName() + ".vcf.gz",
							matchedPanelsOut + "seq-matched-" + chr.getName() + ".vcf.gz",
							matchedPanelsOut + "1kg-seq-merged-" + chr.getName() + ".vcf.gz",
							"/",
							matchedPanelsOut + "mergelog-" + chr.getName() + ".txt",
							true);

					// sort vcfs
					t.sortVCF(linux, vcfsort, matchedPanelsOut + "1kg-matched-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "1kg-matched-sorted-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "sort.sh");
					t.sortVCF(linux, vcfsort, matchedPanelsOut + "seq-matched-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "seq-matched-sorted-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "sort.sh");

					// phase sequencing data
					if (!skipphasing) {
						pb = new ProcessBuilder("java",
								"-Xmx6g",
								"-jar", beagle,
								"ped=" + mergeFam,
								"nthreads=32",
								"gt=" + matchedPanelsOut + "seq-matched-sorted-" + chr.getName() + ".vcf.gz",
								"out=" + matchedPanelsOut + "seqpanel-matched-sorted-phased-" + chr.getName() + ""
						);

						t.run(pb);
						t.sortVCF(linux, vcfsort, matchedPanelsOut + "seqpanel-matched-sorted-phased-" + chr.getName() + ".vcf.gz",
								matchedPanelsOut + "seqpanel-matched-sorted-phased-sorted-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "sort.sh");
					}


					if (!skipimpute) {
						// impute 1kg into reference panel
						String imputedPanelsOut = outputPath + "7-PanelsImputed/";
						Gpio.createDir(imputedPanelsOut);

						// phased data
						pb = new ProcessBuilder("java",
								"-Xmx6g",
								"-jar", beagle,
								"ref=" + matchedPanelsOut + "1kg-matched-sorted-" + chr.getName() + ".vcf.gz",
								//"ped=" + mergeFam,
								"usephase=true", "burnin-its=0", "phase-its=0",
								"gt=" + matchedPanelsOut + "seqpanel-matched-sorted-phased-sorted-" + chr.getName() + ".vcf.gz",
								"nthreads=32",
								"out=" + imputedPanelsOut + "seqpanel-1kgimputed-" + chr.getName() + ""
						);
						t.run(pb);

						pb = new ProcessBuilder("java",
								"-Xmx6g",
								"-jar", beagle,
								"ref=" + matchedPanelsOut + "1kg-matched-sorted-" + chr.getName() + ".vcf.gz",
								"ped=" + mergeFam,
								"gt=" + matchedPanelsOut + "seq-matched-sorted-" + chr.getName() + ".vcf.gz",
								"nthreads=32",
								"out=" + imputedPanelsOut + "seqpanel-unphased-1kgimputed-" + chr.getName() + ""
						);
						t.run(pb);

						t.sortVCF(linux, vcfsort, imputedPanelsOut + "seqpanel-1kgimputed-" + chr.getName() + ".vcf.gz",
								imputedPanelsOut + "seqpanel-1kgimputed-sorted-" + chr.getName() + ".vcf.gz", imputedPanelsOut + "sort.sh");


						// impute reference panel into 1kg
						pb = new ProcessBuilder("java",
								"-Xmx6g",
								"-jar", beagle,
								"ref=" + matchedPanelsOut + "seqpanel-matched-sorted-phased-sorted-" + chr.getName() + ".vcf.gz",
								// "ped=" + mergeFam,
								"usephase=true", "burnin-its=0", "phase-its=0",
								"gt=" + matchedPanelsOut + "1kg-matched-sorted-" + chr.getName() + ".vcf.gz",
								"nthreads=32",
								"out=" + imputedPanelsOut + "1kg-seqpanelimputed-" + chr.getName()
						);
						t.run(pb);
						t.sortVCF(linux, vcfsort, imputedPanelsOut + "1kg-seqpanelimputed-" + chr.getName() + ".vcf.gz",
								imputedPanelsOut + "1kg-seqpanelimputed-sorted-" + chr.getName() + ".vcf.gz", imputedPanelsOut + "sort.sh");


						// merge with original phased datasets..
						String imputedPanelsMergedOut = outputPath + "8-PanelsImputedMerged/";
						Gpio.createDir(imputedPanelsMergedOut);
						String unimputed1kg = matchedPanelsOut + "1kg-matched-sorted-" + chr.getName() + ".vcf.gz";
						String imputed1kg = imputedPanelsOut + "1kg-seqpanelimputed-sorted-" + chr.getName() + ".vcf.gz";
						String unimputed1kgout = imputedPanelsMergedOut + "1kg-imputedVariantsFiltered-" + chr.getName() + ".vcf.gz";
						v.filterVCFVariants(unimputed1kg, imputed1kg, unimputed1kgout); // remove imputed variants from the original vcf if there is any overlap

						String files = unimputed1kgout + " " + imputed1kg;
						String mergedout = imputedPanelsMergedOut + "1kg-seqpanelimputed-mergedWithUnimputed-" + chr.getName() + ".vcf.gz";
						String mergedoutsorted = imputedPanelsMergedOut + "1kg-seqpanelimputed-seqpanelimputed-mergedWithUnimputed-sorted-" + chr.getName() + ".vcf.gz";
						String bashfilename = imputedPanelsMergedOut + "sortmerge.sh";
						System.out.println("concat 1");
						System.out.println(files);
						System.out.println(mergedout);
						System.out.println(mergedoutsorted);
						t.concatVCF(vcfConcat, vcfsort, files, mergedout, mergedoutsorted, bashfilename); // concatenate vcf files

						// repeat on the sequenced data..
						t.sortVCF(linux, vcfsort, imputedPanelsOut + "seqpanel-1kgimputed-" + chr.getName() + ".vcf.gz",
								imputedPanelsOut + "seqpanel-1kgimputed-sorted-" + chr.getName() + ".vcf.gz", imputedPanelsOut + "sort.sh");

						String unimputedseq = matchedPanelsOut + "seqpanel-matched-sorted-phased-sorted-" + chr.getName() + ".vcf.gz";
						String imputedseq = imputedPanelsOut + "seqpanel-1kgimputed-" + chr.getName() + ".vcf.gz";
						String unimputedseqout = imputedPanelsMergedOut + "seqpanel-imputedVariantsFiltered-" + chr.getName() + ".vcf.gz";
						v.filterVCFVariants(unimputedseq, imputedseq, unimputedseqout); // remove imputed variants from original vcf if there is any overlap

						files = unimputedseqout + " " + imputedseq;
						mergedout = imputedPanelsMergedOut + "seqpanel-1kgimputed-mergedWithUnimputed-" + chr.getName() + ".vcf.gz";
						mergedoutsorted = imputedPanelsMergedOut + "seqpanel-1kgimputed-1kgimputed-mergedWithUnimputed-sorted-" + chr.getName() + ".vcf.gz";
						bashfilename = imputedPanelsMergedOut + "sortmerge.sh";
						System.out.println("concat 2");
						System.out.println(files);
						System.out.println(mergedout);
						System.out.println(mergedoutsorted);
						t.concatVCF(vcfConcat, vcfsort, files, mergedout, mergedoutsorted, bashfilename); // concatenate

						// bgzip and index
						System.out.println("bgzip and index");
						v.replaceHeader(imputedPanelsMergedOut + "seqpanel-1kgimputed-1kgimputed-mergedWithUnimputed-sorted-" + chr.getName() + ".vcf.gz",
								imputed1kg,
								imputedPanelsMergedOut + "seqpanel-1kgimputed-1kgimputed-mergedWithUnimputed-sorted-header-" + chr.getName() + ".vcf.gz");

						// merge
						System.out.println("merge");
						String finalImputationPanelsOut = outputPath + "9-ImputationPanels/";
						Gpio.createDir(finalImputationPanelsOut);

						v.mergeAndIntersectVCFVariants(
								imputedPanelsMergedOut + "seqpanel-1kgimputed-1kgimputed-mergedWithUnimputed-sorted-header-" + chr.getName() + ".vcf.gz",
								imputedPanelsMergedOut + "seqpanel-1kgimputed-1kgimputed-mergedWithUnimputed-sorted-" + chr.getName() + ".vcf.gz",
								finalImputationPanelsOut + "1kg-" + chr.getName() + ".vcf.gz",
								finalImputationPanelsOut + "seq-" + chr.getName() + ".vcf.gz",
								finalImputationPanelsOut + "1kg-seq-merged-" + chr.getName() + ".vcf.gz",
								"|",
								finalImputationPanelsOut + "mergelog-" + chr.getName() + ".txt",
								true);


						v.summarizeVCF(finalImputationPanelsOut + "1kg-seq-merged-" + chr.getName() + ".vcf.gz", finalImputationPanelsOut + "1kg-seq-merged-summarystats-" + chr.getName() + ".txt");
					}

				}

			}


		} catch (Exception e) {

			e.printStackTrace();
		}

	}


	public void mergeAndIntersect(boolean linux, int chrint, String vcfsort, String refVCF, String testVCF, String matchedPanelsOut, String separator) throws IOException {
		Chromosome chr = Chromosome.parseChr("" + chrint);
		VCFFunctions v = new VCFFunctions();
		v.mergeAndIntersectVCFVariants(
				refVCF,
				testVCF,
				matchedPanelsOut + "ref-matched-" + chr.getName() + ".vcf.gz",
				matchedPanelsOut + "test-matched-" + chr.getName() + ".vcf.gz",
				matchedPanelsOut + "ref-test-merged-" + chr.getName() + ".vcf.gz",
				separator,
				matchedPanelsOut + "mergelog-" + chr.getName() + ".txt",
				true);

		GenotypeTools t = new GenotypeTools();
		t.sortVCF(linux, vcfsort, matchedPanelsOut + "ref-matched-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "ref-matched-sorted-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "ref-sort-" + chr.getName() + ".sh");
		t.sortVCF(linux, vcfsort, matchedPanelsOut + "test-matched-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "test-matched-sorted-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "test-sort-" + chr.getName() + ".sh");
		t.sortVCF(linux, vcfsort, matchedPanelsOut + "ref-test-merged-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "ref-test-merged-sorted-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "test-sort-" + chr.getName() + ".sh");
	}


	public void prepareGWASDatasets() {

		try {

			String removeTheseIndividuals = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/mergefam.txt";

			String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";
			String[] datasets = new String[]{
//					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_ES_QCgp",
//					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_NL_QCgp",
//					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_SE-E_QCgp",
//					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_SE-U_QCgp",
//					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_UK_QCgp",
//					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_US_QCgp",
			};

			String[] datasetOutput = new String[]{
//					"/Data/ImmunoChip/RA/AllRegions/ES/",
//					"/Data/ImmunoChip/RA/AllRegions/NL/",
//					"/Data/ImmunoChip/RA/AllRegions/SEE/",
//					"/Data/ImmunoChip/RA/AllRegions/SEU/",
//					"/Data/ImmunoChip/RA/AllRegions/UK/",
//					"/Data/ImmunoChip/RA/AllRegions/US/",
			};

			String[] refNames = new String[]{"1kg", "seq", "1kgseq"};
			String[] refSets = new String[]{
					"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/5-1KG-SequencedRegions/filtered.vcf",
					"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/6-PanelsMatched/seqpanel-phased-sorted.vcf",
					"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/9-PanelsMerged/merged.vcfsorted.vcf"};


			for (int i = 0; i < datasets.length; i++) {
				String outdir = datasetOutput[i];
				Gpio.createDir(outdir);
//				liftOverAndFilterPlinkDataset(datasets[i], removeTheseIndividuals, outdir, regionFile);
			}


			datasets = new String[]{
					"/Data/ImmunoChip/T1D/eur",
					"/Data/ImmunoChip/T1D/uk"
			};

			datasetOutput = new String[]{
					"/Data/ImmunoChip/T1D/AllRegions/eur/",
					"/Data/ImmunoChip/T1D/AllRegions/uk/"
			};

			for (int i = 0; i < datasets.length; i++) {
				String outdir = datasetOutput[i];
				Gpio.createDir(outdir);
				liftOverAndFilterPlinkDataset(datasets[i], removeTheseIndividuals, outdir, regionFile);
			}

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public void rsquaredfilter(String vcfIn, String vcfOut) throws IOException {
		TextFile tfin = new TextFile(vcfIn, TextFile.R);
		TextFile tfout = new TextFile(vcfOut, TextFile.W);
		System.out.println();
		System.out.println("Filtering:");
		System.out.println(vcfIn);
		System.out.println(vcfOut);
		String ln = tfin.readLine();
		while (ln != null) {
			if (ln.startsWith("#")) {
				tfout.writeln(ln);
			} else {
				VCFVariant var = new VCFVariant(ln);
				Double ar2 = var.getInfo().get("AR2");
				if (ar2 != null && ar2 >= 0.8) {
					tfout.writeln(ln);
				}
			}
			ln = tfin.readLine();
		}
		tfin.close();
		tfout.close();
	}

	public void datasetmerger() throws IOException {
		String[] datasetOutput = new String[]{
				"/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip/RA/ES/",
				"/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip/RA/NL/",
				"/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip/RA/SEE/",
				"/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip/RA/SEU/",
				"/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip/RA/UK/",
				"/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip/RA/US/",

		};
		String outputdir = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip/RA/merged/";
		mergeGWASDatasets(datasetOutput, outputdir);
		datasetOutput = new String[]{
				"/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip/T1D/eur/",
				"/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip/T1D/uk/"
		};
		outputdir = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip/T1D/merged/";
		mergeGWASDatasets(datasetOutput, outputdir);
	}

	public void mergeGWASDatasets(String[] datasetOutput, String outputdir) throws IOException {


		Gpio.createDir(outputdir);
		for (int chr = 1; chr < 23; chr++) {

			System.out.println("Processing Chr " + chr);
			ArrayList<String> files = new ArrayList<String>();
			for (int d = 0; d < datasetOutput.length; d++) {
				String file = datasetOutput[d] + "genotypes-filtered-sorted-Chr" + chr + ".vcf.gz";
				if (Gpio.exists(file)) {
					files.add(file);
				}
			}

			VCFFunctions v = new VCFFunctions();
			if (files.size() != datasetOutput.length) {
				System.out.println("Found only: " + files.size() + " files for chr: " + chr);
			} else {
				String file1 = files.get(0);
				String file2 = files.get(1);
				v.mergeAndIntersectVCFVariants(file1,
						file2,
						outputdir + "tmp1-Chr" + chr + ".vcf.gz",
						outputdir + "tmp2-Chr" + chr + ".vcf.gz",
						outputdir + "mergedtmp-Chr" + chr + ".vcf.gz",
						"/",
						outputdir + "mergelog-Chr" + chr + "-" + 0 + ".txt.gz",
						false);
				if (files.size() > 2) {
					for (int i = 2; i < files.size(); i++) {
						file2 = files.get(i);
						v.mergeAndIntersectVCFVariants(outputdir + "mergedtmp-Chr" + chr + ".vcf.gz",
								file2,
								outputdir + "tmp1-Chr" + chr + ".vcf.gz",
								outputdir + "tmp2-Chr" + chr + ".vcf.gz",
								outputdir + "merged-Chr" + chr + ".vcf.gz",
								"/",
								outputdir + "mergelog-Chr" + chr + "-" + i + ".txt.gz",
								false);
						Gpio.moveFile(outputdir + "merged-Chr" + chr + ".vcf.gz", outputdir + "mergedtmp-Chr" + chr + ".vcf.gz");
					}
				}
				Gpio.moveFile(outputdir + "mergedtmp-Chr" + chr + ".vcf.gz", outputdir + "merged-Chr" + chr + ".vcf.gz");
			}


		}
	}

	public void liftOverAndFilterPlinkDataset(String plinkDataset, String removeFile, String outdir, String regionFile) throws IOException {

		PedAndMapFunctions pedAndMapFunctions = new PedAndMapFunctions();
		VCFFunctions vcfFunctions = new VCFFunctions();
		ProcessBuilder pb = null;
		GenotypeTools t = new GenotypeTools();

		String map = plinkDataset + ".map";
		String refmap = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC-RsIDsFromT1DStudy.map";
		String chrupdate = plinkDataset + "-raciChrNames.map";
		pedAndMapFunctions.rewriteMapFileChromosomeNames(refmap, map, chrupdate);

		String rsnameref = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC-RsIDsFromT1DStudy.map";
		String rsupdate = plinkDataset + "-raciChrNames-updatedRS.map";
		pedAndMapFunctions.updateMapFileRsIdsUsingMapFile(rsnameref, chrupdate, rsupdate);
//		Gpio.delete(chrupdate);

		String rsupdatededup = plinkDataset + "-raciChrNames-updatedRS-dedup.map";
		pedAndMapFunctions.deduplicateMAP(rsupdate, rsupdatededup);
//		Gpio.delete(rsupdate);
//
		String bedout = plinkDataset + "-raciChrNames-dedup.bed";
		pedAndMapFunctions.rewriteMapToBed(rsupdatededup, bedout);
//

		// liftover
		String lifted = plinkDataset + "-lifted.bed";
		String unlifted = plinkDataset + "-unlifted.bed";
//
		pb = new ProcessBuilder("/Data/Projects/2014-FR-Reseq/ImmunoChip/liftOver", bedout,
				"/Data/Projects/2014-FR-Reseq/ImmunoChip/hg18ToHg19.over.chain.gz", lifted, unlifted);
		t.run(pb);
		System.out.println("Lifted over: " + lifted + " | " + unlifted);
//
//
		String hg19map = plinkDataset + "-raciChrNames-updatedRS-dedup-hg19.map";
		pedAndMapFunctions.convertPostLiftOverMAP(rsupdatededup, hg19map, lifted);
//		Gpio.delete(rsupdatededup);
//

//		Gpio.delete(hg19map);
//
		String hg19mapupd = plinkDataset + "-raciChrNames-hg19-updRS.map";
		String dbsnpvcf = "/Data/Projects/2014-FR-Reseq/GATKResources/dbsnp_138.b37.vcf.gz";
		pedAndMapFunctions.updateRSNames(dbsnpvcf, hg19map, hg19mapupd);
//		Gpio.delete(hg19mapdedup);
//

		// dedup
		String hg19mapdedup = plinkDataset + "-raciChrNames-updatedRS-dedup-hg19-updRS-dedup.map";
		pedAndMapFunctions.deduplicateMAP(hg19mapupd, hg19mapdedup);

		String variantSelect = plinkDataset + "-raciChrNames-hg19-dedup-selectVariants.txt";
//		String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/loci.txt";
		pedAndMapFunctions.filterMap(hg19mapdedup, regionFile, variantSelect);
//
		Gpio.copyFile(plinkDataset + ".map", plinkDataset + ".hg18map");
		Gpio.copyFile(hg19mapdedup, plinkDataset + ".map");

		if (removeFile != null) {
			pb = new ProcessBuilder("/Data/Tools/plink-1.07-mac-intel/plink1.9",
					"--extract", variantSelect,
					"--file", plinkDataset,
					"--recode", "--out", outdir + "genotypes-filtered",
					"--remove", removeFile);
		} else {
			pb = new ProcessBuilder("/Data/Tools/plink-1.07-mac-intel/plink1.9",
					"--extract", variantSelect,
					"--file", plinkDataset,
					"--recode", "--out", outdir + "genotypes-filtered");
		}
		t.run(pb);

		Gpio.copyFile(plinkDataset + ".hg18map", plinkDataset + ".map");
//
		String ped = outdir + "genotypes-filtered";
		String vcf = outdir + "genotypes-filtered.vcf.gz";
		vcfFunctions.convertPEDToVCF(ped, vcf);

		Gpio.delete(ped + ".ped");
		Gpio.delete(ped + ".map");

		// sort
		String sortedvcf = outdir + "genotypes-filtered-sorted.vcf";
		String bashCommand = "gzcat " + vcf + " | /Data/Tools/vcftools/bin/vcf-sort > " + sortedvcf + "\nrm " + sortedvcf + ".gz\ngzip " + sortedvcf;
		String bashfilename = plinkDataset + "-filter-sort.sh";
		TextFile tf = new TextFile(bashfilename, TextFile.W);
		tf.writeln("#!/bin/bash\n" + bashCommand);
		tf.close();
		pb = new ProcessBuilder("bash", bashfilename);
		t.run(pb);

//		Gpio.delete(vcf);
		vcfFunctions.splitPerChromosome(sortedvcf + ".gz", outdir + "genotypes-filtered-sorted");


	}

//	public String preparePlinkDataset(String referenceVCF) throws IOException {
//		outdir += "gt";

	// String referenceVCF = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/MergedWithIC/merged-lowfreqvariantsremoved-sorted-phased-sorted.vcf";
//

//


//
//		// select variants from plink


	// split per chromosome


//		Chromosome[] chromosomes = Chromosome.values();
//		for (Chromosome chr : chromosomes) {
//			// compare VCF to reference dataset
//			vcfFunctions.mergeAndIntersectVCFVariants(
//					referenceVCF + "-" + chr.getName() + ".vcf",
//					outdir + "-filtered-sorted-" + chr.getName() + ".vcf",
//					outdir + "-reference-sorted-matched.vcf",
//					outdir + "-filtered-sorted-" + chr.getName() + "-matched.vcf",
//					outdir + "-filtered-sorted-" + chr.getName() + "-merged.vcf",
//					"/",
//					outdir + "-filtered-sorted-matched-mergelog-" + chr.getName() + ".txt",
//					true);
//		}


//		// next up: beagle compare
//		// impute unfiltered
//		String sortedvcfphased = origmap + "-filtered-sorted-phasedAndImputed";
//		pb = new ProcessBuilder("java",
//				"-jar",
//				"/Data/Tools/beagle/beagle.r1399.jar",
//				"ref=" + referenceVCF,
//				"gt=" + sortedvcf,
//				"out=" + sortedvcfphased,
//				"nthreads=4");
//		run(pb);
//		return null;
//	}

	public void createCovariateMatrixAndMergeFamFiles(String[] famfiles,
													  String[] covariatefiles,
													  String[] cohortnames,
													  String output) throws IOException {

		/*
		 Family ID
     Individual ID
     Paternal ID
     Maternal ID
     Sex (1=male; 2=female; other=unknown)
     Phenotype
		 */

		TextFile mergedFamOut = new TextFile(output + "mergedfam.fam", TextFile.W);
		TextFile mergedDiseaseOut = new TextFile(output + "mergeddisease.txt", TextFile.W);

		HashMap<String, Integer> sampleToSex = new HashMap<String, Integer>();

		HashMap<String, Integer> availableCovariates = new HashMap<String, Integer>();
		ArrayList<String> availableCovariatesList = new ArrayList<String>();


		String mergedHeader = "Sample";
		int covariateCtr = 0;
		for (int cohort = 0; cohort < famfiles.length; cohort++) {
			TextFile famIn = new TextFile(famfiles[cohort], TextFile.R);

			String[] elems = famIn.readLineElems(Strings.whitespace);
			while (elems != null) {

				mergedFamOut.writeln(Strings.concat(elems, Strings.tab));

				String sample = elems[1];
				mergedDiseaseOut.writeln(sample + "\t" + elems[5]);


				int sex = Integer.parseInt(elems[4]);
				if (sex != 1 && sex != 2) {
					sampleToSex.put(sample, -1);
				} else {
					sampleToSex.put(sample, sex);
				}

				elems = famIn.readLineElems(Strings.whitespace);
			}
			famIn.close();

			TextFile covariatefile = new TextFile(covariatefiles[cohort], TextFile.R);

			String[] covariateHeader = covariatefile.readLineElems(Strings.whitespace);


			// fid sid cov1 cov2...
			for (int c = 2; c < covariateHeader.length; c++) {
				String s = covariateHeader[c];
				if (!availableCovariates.containsKey(s)) {
					availableCovariates.put(s, covariateCtr);
					availableCovariatesList.add(s);

					System.out.println("remapping: " + s + "\t" + c + "\t" + covariateCtr);
					mergedHeader += "\t" + s;
					covariateCtr++;
				}
			}
			covariatefile.close();

		}

		mergedFamOut.close();
		mergedDiseaseOut.close();
		mergedHeader += "\tSex";

		for (int cohort = 0; cohort < famfiles.length; cohort++) {
			mergedHeader += "\t" + cohortnames[cohort];
		}


		TextFile mergedCovariatesOut = new TextFile(output + "mergedCovariates.txt", TextFile.W);
		mergedCovariatesOut.writeln(mergedHeader);
		for (int cohort = 0; cohort < famfiles.length; cohort++) {

			String cohortstr = "";
			for (int cohort2 = 0; cohort2 < famfiles.length; cohort2++) {
				if (cohort2 == cohort) {
					cohortstr += "\t1";
				} else {
					cohortstr += "\t0";
				}
			}

			TextFile covariatefile = new TextFile(covariatefiles[cohort], TextFile.R);

			// remap covariates...
			String[] covariateHeader = covariatefile.readLineElems(Strings.whitespace);
			int[] covariateIndex = new int[covariateHeader.length]; // available covariates + sex + cohort will be independently added

			for (int c = 2; c < covariateHeader.length; c++) {
				String covariate = covariateHeader[c];
				covariateIndex[c] = availableCovariates.get(covariate);

			}

			HashMap<String, String> covariatesStrings = new HashMap<String, String>();
			String[] elems = covariatefile.readLineElems(Strings.whitespace);
			HashSet<String> excludeSample = new HashSet<String>();
			while (elems != null) {

				String sample = elems[1];
				String lnOut = sample; // sample

				String[] newCovariateOrder = new String[availableCovariates.size()];
				boolean missingCovariates = false;
				for (int c = 2; c < elems.length; c++) {
					if (elems[c].equals("-9")) {
						missingCovariates = true;
						excludeSample.add(sample);
					}
					newCovariateOrder[covariateIndex[c]] = elems[c];
				}

				if (!missingCovariates) {

					Integer sex = sampleToSex.get(sample);
					String sexStr = "";
					if (sex == null) {
						System.err.println("No gender for: " + sample);
						excludeSample.add(sample);
					} else {
//					if (sex == 1) {
//						sexStr = "1\t0";
//					} else if (sex == 2) {
//						sexStr = "0\t1";
//					} else {
//						sexStr = "0\t0";
//					}


						sexStr = "" + sex;
						lnOut += "\t" + Strings.concat(newCovariateOrder, Strings.tab, "" + 0) + "\t" + sexStr + cohortstr;
						covariatesStrings.put(sample, lnOut);

					}
				}

				elems = covariatefile.readLineElems(Strings.whitespace);
			}
			covariatefile.close();

			String famfile = famfiles[cohort];
			TextFile tffam = new TextFile(famfile, TextFile.R);

			String[] famelems = tffam.readLineElems(Strings.whitespace);
			while (famelems != null) {
				if (famelems.length > 2) {
					String sample = famelems[1];
					Integer sex = sampleToSex.get(sample);
					String covariateStr = covariatesStrings.get(sample);

					if (sex != null) {
						if (covariateStr == null) {
							int[] newCovariateOrder = new int[availableCovariates.size()];
							String sexStr = "" + sex;
							covariateStr = sample + "\t" + Strings.concat(newCovariateOrder, Strings.tab) + "\t" + sexStr + cohortstr;
						}
						if (!excludeSample.contains(sample)) {
							mergedCovariatesOut.writeln(covariateStr);
						}
					}

				} else {
					System.out.println("Could not split line: " + Strings.concat(famelems, Strings.tab));
				}
				famelems = tffam.readLineElems(Strings.whitespace);
			}
			tffam.close();

			//mergedCovariatesOut.writeln(lnOut);
		}

		mergedCovariatesOut.close();
	}

	public void immunoChipRSquaredPlots() throws IOException, DocumentException {

		String sequencedRegionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/2015-03-30-allRegionsWith50PercentOfSamplesAbove20x.bed";
		String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";
		String gtf = "/Data/Annotation/UCSC/genes.gtf";


		String[] refs = new String[]{"seqvaronly", "kg", "seq"};
		String[] datasetIds = new String[]{"RA", "T1D"};
		for (int dataset = 0; dataset < 2; dataset++) {

			for (int chromosome = 1; chromosome < 23; chromosome++) {
				String[] filesForPlotting = new String[refs.length];
				for (int ref = 0; ref < refs.length; ref++) {
					String afterImputation = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/" + datasetIds[dataset] + "/" + refs[ref] + "/merged-variantinfo-Chr" + chromosome + ".txt";
					if (Gpio.exists(afterImputation)) {
						filesForPlotting[ref] = afterImputation;
					}

				}

				String outdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/RsquaredPlots/" + datasetIds[dataset] + "/";
				Gpio.createDir(outdir);
				Chromosome chrin = Chromosome.parseChr("" + chromosome);

				TextFile tf2 = new TextFile(regionFile, TextFile.R);
				String[] elems = tf2.readLineElems(TextFile.tab);
				while (elems != null) {
					Feature region = new Feature();
					region.setChromosome(Chromosome.parseChr(elems[0]));
					region.setStart(Integer.parseInt(elems[1]));
					region.setStop(Integer.parseInt(elems[2]));

					if (region.getChromosome().equals(chrin)) {
						VariantPlot plot = new VariantPlot(outdir + region.toString() + ".pdf", 1400, 1400);
						plot.setMargin(200);
						plot.plotImputationRSquared(filesForPlotting, refs, sequencedRegionFile, region, gtf);

					}
					elems = tf2.readLineElems(TextFile.tab);
				}
				tf2.close();

			}


		}


	}

	public void mergeAssociationResults(String[] names, String[] files, String[] impQualFiles, String outputfilename) throws IOException {


		ArrayList<HashMap<Feature, Pair<String, String>>> dataPVal = new ArrayList<HashMap<Feature, Pair<String, String>>>();
		ArrayList<HashMap<Feature, Pair<String, String>>> dataQualVal = new ArrayList<HashMap<Feature, Pair<String, String>>>();

		HashSet<Feature> uniqueFeats = new HashSet<Feature>();

		// pvalfiles
		for (int f = 0; f < files.length; f++) {

			HashMap<Feature, Pair<String, String>> p = new HashMap<Feature, Pair<String, String>>();
			String file = files[f];

			if (Gpio.exists(file)) {
				boolean isTab = false;
				if (file.endsWith(".tab")) {
					//
					isTab = true;
				}

				TextFile tf = new TextFile(file, TextFile.R);
				tf.readLine();


				String[] elems = tf.readLineElems(Strings.whitespace);
				int ln = 0;
				while (elems != null) {

					if (isTab) {

						String maf = "NaN";
						Double pvald = Double.parseDouble(elems[3]);
						double log10p = -Math.log10(pvald);

						Feature feat = new Feature();
						feat.setChromosome(Chromosome.parseChr(elems[1]));
						Integer pos = Integer.parseInt(elems[2]);
						String name = elems[0];

						feat.setName(name);
						feat.setStart(pos);
						feat.setStop(pos);

						Pair<String, String> pair = new Pair<String, String>(maf, "" + log10p);
						uniqueFeats.add(feat);
						p.put(feat, pair);

					} else {
						String maf = elems[3];
						String pval = elems[elems.length - 1];

						Feature feat = new Feature();
						feat.setChromosome(Chromosome.parseChr(elems[1]));
						Integer pos = Integer.parseInt(elems[2]);
						String name = elems[0];

						feat.setName(name);
						feat.setStart(pos);
						feat.setStop(pos);

						if (name.equals("rs1028182")) {
							System.out.println(Strings.concat(elems, Strings.tab) + ": " + file);

						}


						Pair<String, String> pair = new Pair<String, String>(maf, pval);

						p.put(feat, pair);
						uniqueFeats.add(feat);
					}

					elems = tf.readLineElems(Strings.whitespace);
					ln++;
				}
				tf.close();

				System.out.println(ln + " lines and " + p.size() + " variants parsed.");
			}


			dataPVal.add(p);
		}

		System.out.println(uniqueFeats.size() + " unique variants tested for association...");


		// qual files

		for (int f = 0; f < impQualFiles.length; f++) {
			String file = impQualFiles[f];
			HashMap<Feature, Pair<String, String>> p = new HashMap<Feature, Pair<String, String>>();
			if (file != null && Gpio.exists(file)) {
				TextFile tf = new TextFile(file, TextFile.R);
				String[] elems = tf.readLineElems(Strings.whitespace);
				int ln = 0;
				while (elems != null) {
					if (!elems[0].startsWith("#")) {

						String allele1 = elems[3];
						String allele2 = elems[4];
						String alleles = allele1 + "/" + allele2;

						String arstr = elems[elems.length - 1];
						String[] arstrelems = arstr.split("=");

						if (arstrelems.length == 1) {
							System.out.println("Could not split: " + arstrelems[0] + " for line: " + elems[2]);
						}

						String ar2 = arstrelems[1];
						Feature feat = new Feature();
						feat.setChromosome(Chromosome.parseChr(elems[0]));
						Integer pos = Integer.parseInt(elems[1]);
						String name = elems[2];

						feat.setName(name);
						feat.setStart(pos);
						feat.setStop(pos);


						Pair<String, String> pair = new Pair<String, String>(ar2, alleles);

						p.put(feat, pair);
					}


					elems = tf.readLineElems(Strings.whitespace);
					ln++;
				}
				tf.close();

				System.out.println(ln + " lines and " + p.size() + " variants parsed.");
			}

			dataQualVal.add(p);
		}

		// merge

		TextFile out = new TextFile(outputfilename, TextFile.W);
		// make some kind of header

		String header = "Chr\t" +
				"Pos\t" +
				"Id";

		for (int f = 0; f < files.length; f++) {
			String name = names[f];
			header += "\t" + name + "-MAF\t" + name + "-AR2\t" + name + "-Alleles\t" + name + "-AssocP";
		}

		out.writeln(header);

		for (Feature feat : uniqueFeats) {
			String outln = feat.getChromosome().toString() + "\t" + feat.getStart() + "\t" + feat.getName();
			for (int f = 0; f < files.length; f++) {
				Pair<String, String> quals = dataQualVal.get(f).get(feat);
				Pair<String, String> assoc = dataPVal.get(f).get(feat);
				String maf = "NaN";
				String ar2 = "NaN";
				String alleles = "NaN";
				String assocp = "NaN";
				if (quals != null) {
					alleles = quals.getRight();
					ar2 = quals.getLeft();
				}
				if (assoc != null) {
					maf = assoc.getLeft();
					assocp = assoc.getRight();
				}
				outln += "\t" + maf + "\t" + ar2 + "\t" + alleles + "\t" + assocp;
			}
			out.writeln(outln);
		}
		out.close();


	}


}
