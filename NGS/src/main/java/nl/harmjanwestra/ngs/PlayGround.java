package nl.harmjanwestra.ngs;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.genotypes.GenotypeTools;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.plink.PedAndMapFunctions;
import nl.harmjanwestra.utilities.shell.ProcessRunner;
import nl.harmjanwestra.utilities.vcf.VCFFunctions;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.*;

/**
 * Created by hwestra on 4/21/15.
 */
public class PlayGround {

	public static void main(String[] args) {

		PlayGround g = new PlayGround();



//		g.prepareImputationPanels();
		// g.prepareGWASDatasets();
//		try {
////			PedAndMapFunctions p = new PedAndMapFunctions();
////			p.filterFAM("/Data/Projects/2014-FR-Reseq/2015-finalRun/AllSequencedImmunoChipIDs.txt",
////					"/Data/ImmunoChip/2015-11-28-hg19/Merged/RA.fam",
////					"/Data/ImmunoChip/2015-11-28-hg19/Merged/RA.fam.sequencedIndividuals.txt");
////			g.prepareImputationPanels();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

		PlayGround pg = new PlayGround();
		String plink = "/Data/Homozygosity/EGCUT/EGCUT";
		String outdir = "/Data/Homozygosity/EGCUT/hg19/liftover/EGCUT";
		try {
			pg.liftOver(plink, outdir);
		} catch (IOException e) {
			e.printStackTrace();
		}

		System.exit(-1);

//		Integer chrint = Integer.parseInt(args[0]);
//		String vcfsort = args[1];
//		String refvcf = args[2];
//		String testvcf = args[3];
//		String outdir = args[4];
//		boolean phased = Boolean.parseBoolean(args[5]);
//		try {
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


//
//			g.prepareGWASDatasets();

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

//			g.filterEuropean("/Data/Ref/BeagleRef/1kg.phase3.v5-filtered.merged.vcf.gz",
// "/Data/Ref/BeagleRef/1kg.phase3.v5-filtered.eur.merged.vcf.gz",
// "/Data/Ref/BeagleRef/europeanpopulations.txt");

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

//			String assocdir1 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-22-Assoc/T1D/seq-pseudo/";
//			String assocdir2 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-22-Assoc/T1D/kg-pseudo/";
//			String out1 = "/Data/tmp/assocs.txt";
//			g.matchAssoc(assocdir1, assocdir2, out1);
//
//			String rsqdir1 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/T1D/seqvaronly/";
//			String rsqdir2 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/T1D/kg/";
//			// String out = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/T1D/rsqcompare.txt";
//			String out = "/Data/tmp/rsqcompare.txt";
//			g.matchRSq(rsqdir1, rsqdir2, out);

//
//


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


//			g.associationPlots();
//			g.makeRSquaredPlots();

//			g.mergeAssociationResults(names,files,impqualfiles,out);
//
//			if (args.length < 2) {
//				System.out.println("FilterBadStuff.jar bamin bamout");
//			} else {
//				Resequencing r = new Resequencing();
//				r.filterBAMForImproperMates(args[0], args[1]);
//			}

//			g.readAndCompareXinliPseudoFiles();

//			g.testSingleVCFFileLocalPseudoControls();
//			g.runTest(args);
//			g.testSingleVCFFileLocalPseudoControls();
//
//			g.getPFromTab();
//			g.mergeAssocFiles();
//			g.runTest(args);

		// g.makeRSquaredPlots();
//			g.associationPlots();

//			g.runTest(args);
//			g.testSingleVCFFileLocalPseudoControls();
//			g.getPFromTab();
//			g.associationPlots();

		// g.mergeAssocFiles();
//			g.determineProperlyImputedVariants();
//
//
//			g.mergeAssocFiles();
//			g.august192015();
//			g.datasetmerger();
//			g.august192015();

//			if (args.length == 1) {
//				int threads = Integer.parseInt(args[0]);
//				g.testVCF(threads);
//			} else {
//				System.out.println("Usage: threads");
//			}

//			PosteriorPvalues p = new PosteriorPvalues();
//
//
//
//		} catch (IOException e) {
//			e.printStackTrace();
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
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

	public void august192015() throws IOException {


		// merge ChrX
		VCFFunctions f = new VCFFunctions();
		String in = "/Data/Ref/BeagleRef/chrX.1kg.phase3.v5.vcf.gz";
		String out = "/Data/Ref/BeagleRef/chrX.1kg.phase3.v5-filtered.vcf.gz";
		String regionfile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";

		// filter X chrom in 1kg for regions
//		f.filterVCFForBedRegions(in, out, regionfile);


		GenotypeTools t = new GenotypeTools();
		String vcfsort = "/Data/Tools/vcftools/bin/vcf-sort";
		String vcfconcat = "/Data/Tools/vcftools/bin/vcf-concat";
		String files = "/Data/Ref/BeagleRef/1kg.phase3.v5-filtered.merged.vcf.gz /Data/Ref/BeagleRef/chrX.1kg.phase3.v5-filtered.vcf.gz";
		String merged = "/Data/Ref/BeagleRef/1kg.phase3.v5-filtered.merged-withX.vcf.gz";
		String mergedsorted = "/Data/Ref/BeagleRef/1kg.phase3.v5-filtered.merged-withX-sorted.vcf.gz";
		String bashfilename = "/Data/Ref/BeagleRef/chrXMerge.sh";
//		t.concatVCF(vcfconcat, vcfsort, files, merged, mergedsorted, bashfilename);
//
//		filterEuropean("/Data/Ref/BeagleRef/1kg.phase3.v5-filtered.merged-withX-sorted.vcf.gz",
//				"/Data/Ref/BeagleRef/1kg.phase3.v5-filtered.eur.merged-withX-sorted.vcf.gz",
//				"/Data/Ref/BeagleRef/europeanpopulations.txt");

		// filter european samples


		// prepareGWASDatasets();
//
//
//		for (int i = 1; i < 24; i++) {
//			System.out.println("filtering: " + i);
//			if (i != 23) {
//
//				filterEuropean("/Data/Ref/BeagleRef/chr" + i + ".1kg.phase3.v5.vcf.gz",
//						"/Data/Ref/BeagleRef/chr" + i + ".1kg.phase3.v5.eur.vcf.gz",
//						"/Data/Ref/BeagleRef/europeanpopulations.txt");
//			} else {
//				filterEuropean("/Data/Ref/BeagleRef/chrX.1kg.phase3.v5.vcf.gz",
//						"/Data/Ref/BeagleRef/chrX.1kg.phase3.v5.eur.vcf.gz",
//						"/Data/Ref/BeagleRef/europeanpopulations.txt");
//			}
//
//		}

//		prepareImputationPanels();
		splitsScripts();
		imputationScripts();
	}

	private void splitsScripts() throws IOException {
		String scriptout = "/Data/tmp/split/";
		Gpio.createDir(scriptout);
		String[] dsS = new String[]{"T1D", "RA"};
		String[] chrS = new String[]{"1", "2", "4", "5",
				"6", "7", "9", "10", "11", "12", "13",
				"14", "15", "16", "17", "18", "19", "20", "21", "22"};
		String[] refS = new String[]{"1kg", "seq", "1kg-seq-merged"};

		int nrBatches = 11;

		for (String ref : refS) {
			for (String ds : dsS) {


				for (String chr : chrS) {
					String scriptsh = scriptout + ds + "-" + ref + "-chr" + chr + ".sh";
					String workdir = "/medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/";
					String script = "mkdir -p " + workdir + "" + ds + "/matched/" + ref + "/\n";


					String logdir = "" + workdir + "logs/";

					String logprefix = logdir + ds + "-" + ref + "-chr" + chr + ".log";

					script += "mkdir -p " + logdir + "\n";


					script += "\n\n# merge and check vcf\n";
					script += "nice -n20 java -Xmx10g -jar /medpop/srlab/hwestra/tools/VCFMerge.jar \\\n"
							+ "\ttrue \\\n"
							+ "\t" + chr + " \\\n"
							+ "\t/medpop/srlab/hwestra/tools/vcftools/vcftools_0.1.12b/bin/vcf-sort \\\n"
							+ "\t" + workdir + "ref/" + ref + "-Chr" + chr + ".vcf.gz \\\n"
							+ "\t" + workdir + "" + ds + "/input/merged-Chr" + chr + ".vcf.gz \\\n"
							+ "\t" + workdir + "" + ds + "/matched/" + ref + "/ \\\n"
							+ "\ttab\n";

					script += "\n\n# merge and check vcf\n";
					script += "mkdir -p /medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/" + ds + "/matched/" + ref + "/batches/\n";
					script += "nice -n20 java -Xmx10g -jar /medpop/srlab/hwestra/tools/VCFBatchSplitter.jar \\\n"
							+ "\t/medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/" + ds + "/matched/" + ref + "/test-matched-sorted-Chr" + chr + ".vcf.gz \\\n"
							+ "\t" + workdir + "" + ds + ".fam \\\n"
							+ "\t/medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/" + ds + "/matched/" + ref + "/batches/chr" + chr + "\\\n"
							+ "\t2500\n";
					makeJob(script, scriptsh);
				}

			}
		}
	}

	private void imputationScripts() throws IOException {


		String scriptout = "/Data/tmp/";
		String[] dsS = new String[]{"T1D", "RA"};
		String[] chrS = new String[]{"1", "2", "4", "5",
				"6", "7", "9", "10", "11", "12", "13",
				"14", "15", "16", "17", "18", "19", "20", "21", "22"};
		String[] refS = new String[]{"1kg", "seq", "1kg-seq-merged"};

		int nrBatches = 11;

		for (String ref : refS) {
			for (String ds : dsS) {
				if (ds.equals("T1D")) {
					nrBatches = 11;
				} else {
					nrBatches = 10;
				}
				for (int b = 0; b < nrBatches; b++) {
					for (String chr : chrS) {
						String scriptsh = scriptout + ds + "-" + ref + "-chr" + chr + "-" + b + ".sh";
						String workdir = "/medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/";
						String script = "mkdir -p " + workdir + "" + ds + "/matched/" + ref + "/\n";


						String logdir = "" + workdir + "logs/";

						String logprefix = logdir + ds + "-" + ref + "-chr" + chr + "-batch" + b + ".log";

						script += "mkdir -p " + logdir + "\n";
						script += "rm " + logprefix + " \n";

//						script += "\n\n# merge and check vcf\n";
//						script += "nice -n20 java -Xmx10g -jar /medpop/srlab/hwestra/tools/VCFMerge.jar \\\n"
//								+ "\ttrue \\\n"
//								+ "\t" + chr + " \\\n"
//								+ "\t/medpop/srlab/hwestra/tools/vcftools/vcftools_0.1.12b/bin/vcf-sort \\\n"
//								+ "\t" + workdir + "ref/" + ref + "-Chr" + chr + ".vcf.gz \\\n"
//								+ "\t" + workdir + "" + ds + "/input/merged-Chr" + chr + ".vcf.gz \\\n"
//								+ "\t" + workdir + "" + ds + "/matched/" + ref + "/ \\\n"
//								+ "\ttab >> " + logprefix + " 2>&1";

						// convert to plink format
						script += "\n\n# convert to plink format\n";
						script += "rm " + workdir + "" + ds + "/matched/" + ref + "/chr" + chr + "-batch-" + b + ".ped\n";
						script += "rm " + workdir + "" + ds + "/matched/" + ref + "/chr" + chr + "-batch-" + b + ".map\n";
						script += "nice -n20 /medpop/srlab/hwestra/tools/plink2/plink \\\n"
								+ "\t--vcf " + workdir + "" + ds + "/matched/" + ref + "/batches/chr" + chr + "-batch-" + b + ".vcf.gz \\\n"
								+ "\t--out " + workdir + "" + ds + "/matched/" + ref + "/batches/chr" + chr + "-batch-" + b + " \\\n"
								+ "\t--recode \\\n"
								+ "\t--double-id \\\n"
								+ "\t>> " + logprefix + " 2>&1\n";

						// replace ped headers
						script += "\n\n# replace ped headers\n";
						script += "nice -n20 java -Xmx10g -jar /medpop/srlab/hwestra/tools/replacePEDHeaders.jar \\\n"
								+ "\t" + workdir + "" + ds + "/matched/" + ref + "/batches/chr" + chr + "-batch-" + b + ".ped \\\n"
								+ "\t" + workdir + "" + ds + ".fam \\\n"
								+ "\t>> " + logprefix + " 2>&1\n";

						// check with shapeit
						script += "\n\n# check with shapeit\n";
						script += "rm " + workdir + "" + ds + "/matched/" + ref + "/batches/chr" + chr + "-batch-" + b + ".shapeit.alignments*\n";
						script += "nice -n20 /medpop/srlab/hwestra/tools/shapeit/shapeit \\\n"
								+ "\t-check \\\n"
								+ "\t-P " + workdir + "" + ds + "/matched/" + ref + "/batches/chr" + chr + "-batch-" + b + " \\\n"
								+ "\t-M " + workdir + "imputeref/1000GP_Phase3/genetic_map_chr" + chr + "_combined_b37.txt \\\n"
								+ "\t--input-ref " + workdir + "imputeref/1000GP_Phase3/1000GP_Phase3_chr" + chr + ".hap.gz \\\n"
								+ "\t" + workdir + "imputeref/1000GP_Phase3/1000GP_Phase3_chr" + chr + ".legend.gz \\\n"
								+ "\t" + workdir + "imputeref/1000GP_Phase3/1000GP_Phase3.sample \\\n"
								+ "\t--include-grp " + workdir + "imputeref/1000GP_Phase3/group.list \\\n"
								+ "\t-T 5 \\\n"
								+ "\t--output-log " + workdir + "" + ds + "/matched/" + ref + "/batches/chr" + chr + "-batch-" + b + ".shapeit.alignments \\\n"
								+ "\t>> " + logprefix + " 2>&1\n";

						// phase with reference and pedigree info
						script += "\n\n# phase with reference and pedigree info\n";
						script += "mkdir -p " + workdir + "" + ds + "/phased/" + ref + "\n";
						script += "rm " + workdir + "" + ds + "/phased/" + ref + "/chr" + chr + "-batch-" + b + "*\n";

						script += "if [ -f " + workdir + "" + ds + "/matched/" + ref + "/batches/chr" + chr + "-batch-" + b + ".shapeit.alignments.snp.strand.exclude ]; then \n";

						script += "nice -n20 /medpop/srlab/hwestra/tools/shapeit/shapeit \\\n"
								+ "\t-P " + workdir + "" + ds + "/matched/" + ref + "/batches/chr" + chr + "-batch-" + b + " \\\n"
								+ "\t-M " + workdir + "imputeref/1000GP_Phase3/genetic_map_chr" + chr + "_combined_b37.txt \\\n"
								+ "\t--input-ref " + workdir + "imputeref/1000GP_Phase3/1000GP_Phase3_chr" + chr + ".hap.gz \\\n"
								+ "\t" + workdir + "imputeref/1000GP_Phase3/1000GP_Phase3_chr" + chr + ".legend.gz \\\n"
								+ "\t" + workdir + "imputeref/1000GP_Phase3/1000GP_Phase3.sample \\\n"
								+ "\t--include-grp " + workdir + "imputeref/1000GP_Phase3/group.list \\\n"
								+ "\t--exclude-snp " + workdir + "" + ds + "/matched/" + ref + "/batches/chr" + chr + "-batch-" + b + ".shapeit.alignments.snp.strand.exclude \\\n"
								+ "\t--duohmm \\\n"
								+ "\t-O " + workdir + "" + ds + "/phased/" + ref + "/chr" + chr + "-batch-" + b + " \\\n"
								+ "\t-T 5 \\\n"
								+ "\t--output-log " + workdir + "" + ds + "/phased/" + ref + "/chr" + chr + "-batch-" + b + ".shapeit.phasing \\\n"
								+ "\t>> " + logprefix + " 2>&1\n";

						script += " else \n";

						script += "nice -n20 /medpop/srlab/hwestra/tools/shapeit/shapeit \\\n"
								+ "\t-P " + workdir + "" + ds + "/matched/" + ref + "/batches/chr" + chr + "-batch-" + b + " \\\n"
								+ "\t-M " + workdir + "imputeref/1000GP_Phase3/genetic_map_chr" + chr + "_combined_b37.txt \\\n"
								+ "\t--input-ref " + workdir + "imputeref/1000GP_Phase3/1000GP_Phase3_chr" + chr + ".hap.gz \\\n"
								+ "\t" + workdir + "imputeref/1000GP_Phase3/1000GP_Phase3_chr" + chr + ".legend.gz \\\n"
								+ "\t" + workdir + "imputeref/1000GP_Phase3/1000GP_Phase3.sample \\\n"
								+ "\t--include-grp " + workdir + "imputeref/1000GP_Phase3/group.list \\\n"
								+ "\t--duohmm \\\n"
								+ "\t-O " + workdir + "" + ds + "/phased/" + ref + "/chr" + chr + "-batch-" + b + " \\\n"
								+ "\t-T 5 \\\n"
								+ "\t--output-log " + workdir + "" + ds + "/phased/" + ref + "/chr" + chr + "-batch-" + b + ".shapeit.phasing \\\n"
								+ "\t>> " + logprefix + " 2>&1\n";

						script += " fi \n";

						// convert back to VCF
						script += "\n\n# convert back to VCF\n";
						script += "nice -n20 java -Xmx10g -jar /medpop/srlab/hwestra/tools/convertHapsToVCF.jar \\\n"
								+ "\t" + workdir + "" + ds + "/phased/" + ref + "/chr" + chr + "-batch-" + b + " \\\n"
								+ "\t" + workdir + "" + ds + "/phased/" + ref + "/chr" + chr + "-batch-" + b + " \\\n"
								+ "\ttrue \\\n"
								+ "\t/medpop/srlab/hwestra/tools/vcftools/vcftools_0.1.12b/bin/vcf-sort \\\n"
								+ "\t>> " + logprefix + " 2>&1\n";

//					// impute-ah!
						script += "\n\n# impute\n";
						script += "mkdir -p /medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/" + ds + "/imputed/" + ref + "\n";
						script += "rm /medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/" + ds + "/imputed/" + ref + "/chr" + chr + "-batch-" + b + "*\n";
						script += "nice -n20 java -Xmx10g -jar /medpop/srlab/hwestra/tools/beagle.r1399.jar \\\n"
								+ "\tref=/medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/ref/" + ref + "-Chr" + chr + ".vcf.gz \\\n"
								+ "\tusephase=true \\\n"
								+ "\tburnin-its=0 \\\n"
								+ "\tphase-its=0 \\\n"
								+ "\tgt=/medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/" + ds + "/phased/" + ref + "/chr" + chr + "-batch-" + b + "-sorted.vcf.gz \\\n"
								+ "\tnthreads=5 \\\n"
								+ "\tout=/medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/" + ds + "/imputed/" + ref + "/chr" + chr + "-batch-" + b + " \\\n"
								+ "\t>> " + logprefix + " 2>&1\n";

						makeJob(script, scriptsh);
					}
				}

			}
		}


	}

	public void determineProperlyImputedVariants() throws IOException {
		String startvcf1 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-30-AllChr-EurOnly/2-SeqVariantRegionFilter/filtered.vcf.gz";

		VCFFunctions functions = new VCFFunctions();
		ArrayList<Feature> variants = functions.getVariantsFromVCF(startvcf1, true, false);
		HashSet<Feature> featureHash = new HashSet<Feature>();
		featureHash.addAll(variants);

		TextFile out = new TextFile("/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-22-Assoc/ImputedVariants/SeqPanelInto1KG.txt", TextFile.W);

		for (int i = 1; i < 23; i++) {
			String imputationResults = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/ImputationPanels/ref-impquals-Chr" + i + ".txt";

			String preimputation = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/ImputationPanels/1kg-preimputation-input-Chr" + i + ".txt";
			HashSet<Feature> alreadyindataset = new HashSet<Feature>();
			TextFile tfin = new TextFile(preimputation, TextFile.R);
			String[] preimpElems = tfin.readLineElems(Strings.whitespace);
			while (preimpElems != null) {
				if (!preimpElems[0].startsWith("#")) {
					Feature f = new Feature();
					// 1 1955373 rs4648791 C T . . .
					Chromosome chr = Chromosome.parseChr(preimpElems[0]);
					Integer stat = Integer.parseInt(preimpElems[1]);
					String name = preimpElems[2];
					f.setStart(stat);
					f.setStop(stat);
					f.setChromosome(chr);
					f.setName(name);
					alreadyindataset.add(f);

				}
				preimpElems = tfin.readLineElems(Strings.whitespace);
			}
			tfin.close();

			TextFile tf = new TextFile(imputationResults, TextFile.R);
			String[] elems = tf.readLineElems(Strings.whitespace);
			int ln = 0;
			while (elems != null) {
				if (!elems[0].startsWith("#")) {

					String allele1 = elems[3];
					String allele2 = elems[4];
					String alleles = allele1 + "/" + allele2;

					String info = elems[elems.length - 1];
					String[] infoelems = info.split(";");
					double ar2 = 0;
					for (String s : infoelems) {
						String[] selems = s.split("=");
						if (selems[0].equals("AR2")) {
							ar2 = Double.parseDouble(selems[1]);
						}
					}

					Feature feat = new Feature();
					feat.setChromosome(Chromosome.parseChr(elems[0]));
					Integer pos = Integer.parseInt(elems[1]);
					String name = elems[2];

					feat.setName(name);
					feat.setStart(pos);
					feat.setStop(pos);


					if (featureHash.contains(feat)) {

						if (ar2 >= 0.8) {
							boolean alreadyInDataset = alreadyindataset.contains(feat);
							out.writeln(Strings.concat(elems, Strings.tab) + "\t" + ar2 + "\t" + alreadyInDataset);
						}
					}


				}


				elems = tf.readLineElems(Strings.whitespace);
				ln++;
			}
			tf.close();
		}
		out.close();

	}

	public void getPFromTab() throws IOException {
		String onen = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/T1D-Onengut/hg19_gwas_ic_t1d_onengut_meta_4_18_1.tab";
		TextFile tf = new TextFile(onen, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, Double> d = new HashMap<String, Double>();
		while (elems != null) {

			String snp = elems[0];
			String p = elems[elems.length - 2];
			try {
				Double d2 = Double.parseDouble(p);
				d2 = -(Math.log10(d2));//-Math.log(d2);

				d.put(snp, d2);
			} catch (NumberFormatException e) {
				System.out.println(p + " is not a number for " + snp);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		String list = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/T1D-Onengut/test.txt";
		String outf = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/T1D-Onengut/test-p.txt";
		TextFile tf2 = new TextFile(list, TextFile.R);
		TextFile out = new TextFile(outf, TextFile.W);
		String ln = tf2.readLine();
		while (ln != null) {
			out.writeln("" + d.get(ln.trim()));

			ln = tf2.readLine();
		}


		out.close();
		tf2.close();

	}


	public void readAndCompareXinliPseudoFiles() throws IOException {
		String hjfile = "/Data/tmp/pseudotest/-pseudosamplelist.txt";
		String xinxinfile = "/Data/tmp/xinli/pseudocontrols.txt";

		String famfile = "/Data/ImmunoChip/T1D/binary/eur.fam";
		TextFile tf3 = new TextFile(famfile, TextFile.R);
		String[] elems3 = tf3.readLineElems(Strings.whitespace);
		HashSet<String> allfams = new HashSet<String>();
		while (elems3 != null) {
			String family = elems3[0];
			allfams.add(family);
			elems3 = tf3.readLineElems(Strings.whitespace);
		}
		tf3.close();

		System.out.println(allfams.size() + " all families");

		HashSet<String> families = new HashSet<String>();


		TextFile tf = new TextFile(xinxinfile, TextFile.R);
		String[] elems = tf.readLineElems(Strings.whitespace);

		while (elems != null) {

			String sample1 = elems[0];
			String sample2 = elems[1];
			if (sample2.endsWith("_ctrl")) {
				families.add(sample1);
			}


			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();
		System.out.println(families.size() + " families");

		TextFile tf2 = new TextFile(hjfile, TextFile.R);


		HashSet<String> familieshj = new HashSet<String>();

		String[] elems2 = tf2.readLineElems(TextFile.tab);
		while (elems2 != null) {
			String sample = elems2[0];
			if (sample.endsWith("-PseudoControl")) {

				String samplefam = sample.replaceAll("-PseudoControl", "");
				samplefam = samplefam.substring(0, samplefam.length() - 2);

				familieshj.add(samplefam);


			}
			elems2 = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();

		for (String s : families) {
			if (!familieshj.contains(s)) {
				System.out.println(s + " not found");
			}
		}
	}

	public void mergeAssocFiles() throws IOException {

		String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/2015-03-30-allRegionsWith50PercentOfSamplesAbove20xPlus1K.bed"; // "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";
		boolean targetregionsonly = false;

		for (int chr = 1; chr < 23; chr++) {
			Chromosome chrin = Chromosome.parseChr("" + chr);

			String[] files = null;
			String[] impquals = null;
			String[] names = null;

			for (int d = 0; d < 2; d++) {
				if (d == 0) {
					files = new String[]{
							"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/T1D-Onengut/hg19_gwas_ic_t1d_onengut_meta_4_18_1.tab",
							"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-22-Assoc/T1D/before-UK/Chr" + chr + "-gwas.txt",
							"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-22-Assoc/T1D/kg-UK/Chr" + chr + "-gwas.txt",
							"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-22-Assoc/T1D/seq-UK/Chr" + chr + "-gwas.txt"
					};
					names = new String[]{"Onengut", "ImmunoChip", "1000Genomes", "Sequencing"};
					impquals = new String[]{null,
							null,
							"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/T1D/kg/merged-variantinfo-Chr" + chr + ".txt",
							"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/T1D/seq/merged-variantinfo-Chr" + chr + ".txt"};
				} else {
					files = new String[]{
							"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/RA-Eyre/hg19_gwas_ic_ra_eyre_4_18_0.tab",
							"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-22-Assoc/RA/before/Chr" + chr + "-gwas.txt",
							"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-22-Assoc/RA/kg/Chr" + chr + "-gwas.txt",
							"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-22-Assoc/RA/seq/Chr" + chr + "-gwas.txt"
					};
					names = new String[]{"Eyre", "ImmunoChip", "1000Genomes", "Sequencing"};
					impquals = new String[]{null,
							null,
							"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/RA/kg/merged-variantinfo-Chr" + chr + ".txt",
							"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/RA/seq/merged-variantinfo-Chr" + chr + ".txt"};
				}


				Feature chrfeat = new Feature();
				chrfeat.setChromosome(chrin);
				chrfeat.setStart(0);
				chrfeat.setStop(chrin.getLength());
				if (d == 0) {
					String chrout = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-22-Assoc/AssociationOutput/AssociationFiles/T1D/";
					mergeAssociationResults(chrfeat, names, files, impquals, chrout);
				} else {
					String chrout = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-22-Assoc/AssociationOutput/AssociationFiles/RA/";
					mergeAssociationResults(chrfeat, names, files, impquals, chrout);
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
							String out = null;
							if (d == 0) {
								out = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-22-Assoc/AssociationOutput/AssociationFiles/T1D/T1D-UK-merged-Chr" + chr + "_" + region.getStart() + "-" + region.getStop() + ".txt";
							} else {
								out = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-22-Assoc/AssociationOutput/AssociationFiles/RA/RA-merged-Chr" + chr + "_" + region.getStart() + "-" + region.getStop() + ".txt";
							}

							mergeAssociationResults(region, names, files, impquals, out);
						}
					}
					elems = tf2.readLineElems(TextFile.tab);
				}
				tf2.close();
			}


		}
	}

	public void runTest(String[] args) throws IOException {
		PlayGround g = new PlayGround();
		if (args.length < 4) {
			System.out.println("usage: d r theads pseudo [start] [stop]");
			System.out.println("d0 = T1D");
			System.out.println("d1 = RA");
			System.out.println("before kg seq kgseq seqvaronly");
		} else {


			int start = 1;
			int stop = 23;
			if (args.length >= 6) {
				start = Integer.parseInt(args[4]);
				stop = Integer.parseInt(args[5]);
			}

			// g.testVCF(Integer.parseInt(args[0]), Integer.parseInt(args[1]), Integer.parseInt(args[2]), Boolean.parseBoolean(args[3]), start, stop);
			//"before", "kg", "seq", "kgseq", "seqvaronly"
		}
	}

	public void runCoverage(String[] args) throws IOException {
		if (args.length < 5) {
			System.out.println("Usage: java -Xmx10g -jar Coverage.jar listoffiles.txt outdir targets.bed makebedgraps nrthreads\n" +
					"makebedgraphs should be a boolean: true/false" +
					"if makebedgraphs is enabled, overlapping regions will be merged into one single region" +
					"nrthreads should be an integer [1,)" +
					"other variables should be paths" +
					"listoffiles.txt format: samplename\\tpath/to/bamfile.bam");
		} else {
			Coverage c = new Coverage();
			String listfile = args[0];
			String outdir = args[1];
			String targets = args[2];
			Boolean b = Boolean.parseBoolean(args[3]);
			int nrthreads = Integer.parseInt(args[4]);
			if (nrthreads <= 0) {
				nrthreads = 1;
			}

			c.bamToBedWithinRegionsForList(listfile, outdir, targets, b, nrthreads);
		}

		System.out.println("Have a nice day.");
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
			String f = dir + "Chr" + chr + "-gwas.txt";
			if (Gpio.exists(f)) {

				System.out.println("loading: " + f);

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
		System.out.println(assocValues.size() + " assocs loaded");
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

	private void matchAssoc(String assocdir1, String assocdir2, String out) throws IOException {

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

		TextFile outf = new TextFile(out, TextFile.W);

		for (int i = 0; i < xarr.length; i++) {
			outf.writeln(xarr[i] + "\t" + yarr[i]);
		}

		outf.close();

		Grid g = new Grid(500, 500, 1, 1, 100, 100);


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


		TextFile outfq = new TextFile(out + "-rsquareds.txt", TextFile.W);

		for (int i = 0; i < xarr.length; i++) {
			outfq.writeln(xarr[i] + "\t" + yarr[i]);
		}

		outfq.close();

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

				TextFile associn = new TextFile(assoc + "Chr" + chr + "-gwas.txt", TextFile.R);
				TextFile assocout = new TextFile(assocOut + "Chr" + chr + "-gwas.txt", TextFile.W);
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

	void mergeImputedBatches(String prefix, String outfilename, int nrbatches, int chr) throws IOException {
		VCFFunctions f = new VCFFunctions();

//		f.mergeImputationBatches(prefix, outfilename, nrbatches, chr);
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

//	public void splitVCFs() throws IOException {
//
//		VCFFunctions f = new VCFFunctions();
//
//		for (int reference = 3; reference < 4; reference++) {
//			String refStr = "kg";
//			if (reference == 1) {
//				refStr = "kgseq";
//			} else if (reference == 2) {
//				refStr = "seq";
//			} else if (reference == 3) {
//				refStr = "seqvaronly";
//			}
//			for (int dataset = 0; dataset < 2; dataset++) {
//
//				String datasetStr = "RA";
//				String famfile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/RA.fam";
//				if (dataset == 1) {
//					datasetStr = "T1D";
//					famfile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/T1D.fam";
//				}
//
//				for (int chromosome = 1; chromosome < 23; chromosome++) {
//
//					String input = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/" + datasetStr + "/matched/" + refStr + "/test-matched-sorted-Chr" + chromosome + ".vcf.gz";
//					if (reference == 3) {
//						input = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/" + datasetStr + "/matched/" + refStr + "/matched-test-matched-sorted-Chr" + chromosome + ".vcf.gz";
//					}
//					if (Gpio.exists(input)) {
//						String outputprefix = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/" + datasetStr + "/matched/" + refStr + "/split/Chr" + chromosome + "/";
//						Gpio.createDir(outputprefix);
//						outputprefix += "test-match-sorted-Chr" + chromosome;
//						f.splitVCFOverRandomBatches(input, famfile, outputprefix, 1000);
//					} else {
//						System.out.println("cannot find: " + input);
//					}
//				}
//
//			}
//		}
//	}


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

//	private void associationPlots() throws IOException {
//		try {
//			String sequencedRegionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/2015-03-30-allRegionsWith50PercentOfSamplesAbove20x.bed";
//			String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/2015-03-30-allRegionsWith50PercentOfSamplesAbove20xPlus1K.bed"; // "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";
////			String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/interesting.bed";
//			String gtf = "/Data/Annotation/UCSC/genes.gtf";
//
//			String assocInputDir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-22-Assoc/";
//			double mafthreshold = 0.005;
//			boolean logtransform = false;
//			boolean targetregionsonly = false;
//			if (targetregionsonly) {
//				regionFile = sequencedRegionFile;
//			}
//
//			String[] referenceSets = new String[]{"ICStudy", "before", "kg", "seq"};
//			String[] referenceNames = new String[]{"ImmunoBase", "ImmunoChip", "1000Genomes", "Sequencing"};
//			for (int dataset = 0; dataset < 2; dataset++) {
//				String datasetStr = "RA";
//				if (dataset == 1) {
//					datasetStr = "T1D";
//				}
//
//				for (int chr = 1; chr < 23; chr++) {
//					Chromosome chrin = Chromosome.parseChr("" + chr);
//					String[] datasetFiles = new String[referenceSets.length];
//
////					if (datasetStr.equals("T1D")) {
////						// datasetFiles = new String[(referenceSets.length * 2) - 1];
////
////						referenceNames = new String[(referenceSets.length * 2) - 1];
////						for (int r = 0; r < referenceSets.length; r++) {
////							if (r == referenceSets.length - 1) {
////								referenceNames[referenceNames.length - 1] = referenceSets[r];
////							} else {
////								referenceNames[r * 2] = referenceSets[r];
////								referenceNames[r] = referenceSets[r] + "-pseudo";
////							}
////						}
////					} else {
////
////					}
//					String[] refNames = new String[referenceNames.length];
//					for (int i = 0; i < refNames.length; i++) {
//						refNames[i] = new String(referenceNames[i]);
//					}
//
//
//					String output = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-22-Assoc/plots/association/" + datasetStr + "/";
//					Gpio.createDir(output);
//
//					for (int reference = 0; reference < referenceNames.length; reference++) {
//						String refStr = referenceSets[reference];
//						if (dataset == 0 && refStr.equals("ICStudy")) {
//							datasetFiles[reference] = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/RA-Eyre/hg19_gwas_ic_ra_eyre_4_18_0.tab";
//							refNames[reference] = "Eyre et al";
//						} else if (refStr.equals("ICStudy")) {
//							datasetFiles[reference] = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/T1D-Onengut/hg19_gwas_ic_t1d_onengut_meta_4_18_1.tab";
//							refNames[reference] = "Onengut et al";
//						} else if (dataset == 1) {
//							String afterImputation = assocInputDir + datasetStr + "/" + refStr + "-UK/" + chrin.toString() + "-gwas.txt";
//							datasetFiles[reference] = afterImputation;
//						} else {
//							String afterImputation = assocInputDir + datasetStr + "/" + refStr + "/" + chrin.toString() + "-gwas.txt";
//							datasetFiles[reference] = afterImputation;
//						}
//					}
//
//					for (int i = 0; i < datasetFiles.length; i++) {
//						System.out.println(dataset + "\t" + datasetFiles[i]);
//					}
//
//
//					TextFile tf2 = new TextFile(regionFile, TextFile.R);
//					String[] elems = tf2.readLineElems(TextFile.tab);
//
//					while (elems != null) {
//						Feature region = new Feature();
//						if (elems.length > 2) {
//							if (targetregionsonly) {
//								region.setChromosome(Chromosome.parseChr(elems[0]));
//								region.setStart(Integer.parseInt(elems[1]) - 10000);
//								region.setStop(Integer.parseInt(elems[2]) + 10000);
//							} else {
//								region.setChromosome(Chromosome.parseChr(elems[0]));
//								region.setStart(Integer.parseInt(elems[1]));
//								region.setStop(Integer.parseInt(elems[2]));
//							}
//							if (region.getChromosome().equals(chrin)) {
//
//								String plotout = output + region.toString() + ".pdf";
//
//								int height = 500 + (referenceNames.length * 250);
//
//
//								VariantPlot plot = new VariantPlot(plotout, 1400, height);
//								plot.setMargin(200);
//
//								plot.plotAssocPvalue(gtf,
//										datasetFiles,
//										refNames,
//										sequencedRegionFile,
//										region,
//										mafthreshold,
//										logtransform);
//							}
//						}
//						elems = tf2.readLineElems(TextFile.tab);
//					}
//
//					tf2.close();
//				}
//			}
//
//
//			// String[] variantFiles, String[] variantFileNames, String sequencedRegionFile, String regionFile
//
//
////			plot.close();
//
//		} catch (DocumentException e) {
//			e.printStackTrace();
//		}
//
//	}

//	public void mergeImmunoChipVCF() throws IOException {
//
//		String data1 = "/Data/ImmunoChip/T1D/AllRegions/eur/genotypes-filtered-sorted-";
//		String data2 = "/Data/ImmunoChip/T1D/AllRegions/uk/genotypes-filtered-sorted-";
//
//		String outdir = "/Data/ImmunoChip/T1D/AllRegions/merged/genotypes-filtered-sorted-";
//
//		VCFFunctions v = new VCFFunctions();
//		for (int i = 1; i < 23; i++) {
//
//			String ref = data1 + "Chr" + i + ".vcf.gz";
//			String test = data2 + "Chr" + i + ".vcf.gz";
//			if (Gpio.exists(ref)) {
//				String refout = outdir + "-ref.vcf.gz";
//				String testout = outdir + "-ref.vcf.gz";
//				String mergedout = outdir + "Chr" + i + ".vcf.gz";
//				String log = outdir + "Chr" + i + ".log.gz";
//				v.mergeAndIntersectVCFVariants(ref, test, refout, testout, mergedout, "/", log, false);
//			}
//
//		}
//
//	}

//	public void testSingleVCFFileLocalPseudoControls() {
//
//		String vcf = "/Data/tmp/pseudotest/merged-filtered.Chr7-first200.vcf";
//		String out = "/Data/tmp/pseudotest/R/";
//
//
//		String diseasestatus = "/Data/ImmunoChip/T1D/binary/covarmerged.txtmergeddisease.txt";
//		String covariatefile = "/Data/ImmunoChip/T1D/binary/covarmerged.txtmergedCovariatesWoM.txt";
//		HashSet<String> covariatestoinclude = null;
//		HashSet<String> snpLimit = null;
//		String samplestoexclude = null;
//		String famfile = "/Data/ImmunoChip/T1D/binary/covarmerged.txtmergedfam.fam";
//
//		boolean filterimputationquality = true;
//		double imputationqualitythreshold = 0.8;
//		int mafthreshold = 5;
//		int threadnum = 1;
//
//		LogitTestR testObj = new LogitTestR(vcf,
//				out,
//				diseasestatus,
//				covariatefile,
//				covariatestoinclude,
//				snpLimit,
//				samplestoexclude,
//				null,
//				filterimputationquality,
//				imputationqualitythreshold,
//				mafthreshold, threadnum);
//
//		out = out + "-pseudo";
//
//		LogitTestR testObj2 = new LogitTestR(vcf,
//				out,
//				diseasestatus,
//				covariatefile,
//				covariatestoinclude,
//				snpLimit,
//				samplestoexclude,
//				famfile,
//				filterimputationquality,
//				imputationqualitythreshold,
//				mafthreshold, threadnum);
//
//		try {
//			testObj.call();
////			testObj2.call();
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
//
//
//	}

//	public void testVCFLocal(int nrThreads) throws IOException {
//		String covariatefile = "";
//		String outputdir = "/Data/tmp/pseudocontroltest/";
//		String diseasestatus = "";
//
//		HashSet<String> covariatestoinclude = null;
//		String samplestoexclude = null;
//		boolean filterimputationquality = false;
//		double imputationqualitythreshold = 0.8;
//		int mafthreshold = 5;
//
//		System.out.println("Opening threadpool for " + nrThreads + " threads.");
//		ExecutorService threadPool = Executors.newFixedThreadPool(nrThreads);
//		CompletionService<Boolean> pool = new ExecutorCompletionService<Boolean>(threadPool);
//
//		HashSet<String> snpLimit = null;
//
//		int submit = 0;
//		int start = 1;
//		int stop = 23;
//		String famfile = null;
//		for (int d = 0; d < 1; d++) {
//
////			if (d == 0) {
////				start = 1;
////				stop = 23;
////			} else {
////				start = 1;
////				stop = 23;
////			}
//
//			for (int chr = start; chr < stop; chr++) {
//				String vcf = "";
//				String out = "";
//				for (int reference = 2; reference < 3; reference++) {
//					String ref = "beforeImputation";
//
//
//					if (d == 0) {
//						covariatefile = "/Data/ImmunoChip/T1D/binary/covarmerged.txtmergedCovariates.txt";
//						//vcf = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationOutput/T1D/merged-filtered-Chr" + chr + ".vcf.gz";
//
//						vcf = "/Data/ImmunoChip/T1D/merged/merged-Chr" + chr + ".vcf.gz";
//						diseasestatus = "/Data/ImmunoChip/T1D/binary/covarmerged.txtmergeddisease.txt";
//						out = outputdir + "T1D/" + ref + "/";
//						famfile = "/Data/ImmunoChip/T1D/binary/covarmerged.txtmergedfam.fam";
//					} else {
//						covariatefile = "/Data/ImmunoChip/RA/covar/covarmerged.txtmergedCovariates.txt";
//						vcf = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationOutput/RA/merged-filtered-Chr" + chr + ".vcf.gz";
//						diseasestatus = "/Data/ImmunoChip/RA/covar/covarmerged.txtmergeddisease.txt";
//						out = outputdir + "RA/" + ref + "/";
//						famfile = "/Data/ImmunoChip/RA/covar/covarmerged.txtmergedfam.fam";
//					}
//					Gpio.createDir(out);
//
//
//					String pseudoOut = "";
//					if (d == 0) {
//						pseudoOut = out + "/pseudo/";
//						Gpio.createDir(pseudoOut);
//						pseudoOut += "Chr" + chr + "-";
//					}
//					out += "Chr" + chr + "-";
//
//					if (chr != 8 && chr != 3) {
//						if (Gpio.exists(vcf)) {
//							if (d == 0) {
//								LogitTestR testObj = new LogitTestR(vcf,
//										out,
//										diseasestatus,
//										covariatefile,
//										covariatestoinclude,
//										snpLimit,
//										samplestoexclude,
//										null,
//										filterimputationquality,
//										imputationqualitythreshold,
//										mafthreshold, submit);
//
//								pool.submit(testObj);
//								submit++;
//								testObj = new LogitTestR(vcf,
//										pseudoOut,
//										diseasestatus,
//										covariatefile,
//										covariatestoinclude,
//										snpLimit,
//										samplestoexclude,
//										famfile,
//										filterimputationquality,
//										imputationqualitythreshold,
//										mafthreshold, submit);
//
//								pool.submit(testObj);
//								submit++;
//							} else {
//								LogitTestR testObj = new LogitTestR(vcf,
//										out,
//										diseasestatus,
//										covariatefile,
//										covariatestoinclude,
//										snpLimit,
//										samplestoexclude,
//										famfile,
//										filterimputationquality,
//										imputationqualitythreshold,
//										mafthreshold, submit);
//
//								pool.submit(testObj);
//								submit++;
//							}
//
//						} else {
//							System.err.println("file does not exist: " + vcf);
//						}
//					}
//				}
//
//			}
//
//		}
//
//		int returned = 0;
//		while (returned < submit) {
//
//			try {
//
//
//				Boolean result = pool.take().get();
//
//				if (result) {
//					returned++;
//				} else {
//
//					System.exit(-1);
//				}
//
//
//				System.out.println(returned + " / " + submit + " returned");
//			} catch (Exception e) {
//				e.printStackTrace();
//			}
//
//		}
//
//		threadPool.shutdown();
//
//	}

	//	public void testVCF(int d, int reference, int nrThreads, boolean pseudo, int start, int stop) throws IOException {
//	public void testVCF(int nrThreads) throws IOException {
//
////		String vcf = "/Data/tmp/test-matched-Chr10.vcf.gz";
////		String covariatefile = "/Data/ImmunoChip/RA/covar/covarmerged.txtmergedCovariates.txt";
////		String outputdir = "/Data/tmp/";
////		String diseasestatus = "/Data/ImmunoChip/RA/covar/covarmerged.txtmergeddisease.txt";
////
////		HashSet<String> covariatestoinclude = null;
////		String samplestoexclude = null;
////		boolean filterimputationquality = false;
////		double imputationqualitythreshold = 0.8;
////		double mafthreshold = 0.005;
////
////		test.logitTest(vcf, outputdir, diseasestatus, covariatefile, covariatestoinclude, samplestoexclude, filterimputationquality, imputationqualitythreshold, mafthreshold);
////"/Data/ImmunoChip/T1D/binary/plinkvcf/plink.vcf";
//
//		String covariatefile = "";
//		String outputdir = "/medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/Assoc/";
//		String diseasestatus = "";
//
//		HashSet<String> covariatestoinclude = null;
//		String samplestoexclude = null;
//		boolean filterimputationquality = true;
//		double imputationqualitythreshold = 0.8;
//		int mafthreshold = 5;
//
//		System.out.println("Opening threadpool for " + nrThreads + " threads.");
//		ExecutorService threadPool = Executors.newFixedThreadPool(nrThreads);
//		CompletionService<Boolean> pool = new ExecutorCompletionService<Boolean>(threadPool);
//
//
//		String bedfile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/allLoci.bed";
//
//		BedFileReader reader = new BedFileReader();
//		ArrayList<Feature> bedregions = reader.readAsList(bedfile);
//		System.out.println(bedregions.size() + " bed regions loaded");
//
//
//		int submit = 0;
//		String famfile = null;
////		String[] refs = new String[]{"before", "kg", "seq", "kgseq", "seqvaronly"};
//		String[] refs = new String[]{"1kg", "seq", "1kg-seq-merged"};
////		String[] refs = new String[]{"1kg", "1kg-seq-merged"};
//
//		int endD = 2;
////		if (pseudo) {
////			endD = 1;
////		}
//
////		int[] chrs = new int[]{9};
//		for (int reference = 0; reference < refs.length; reference++) {
//
//
//			for (int d = 0; d < endD; d++) {
//
//				for (int chr = 1; chr < 23; chr++) {
////				for (int chr : chrs) {
//					String vcf = "";
//					String out = "";
////					for (int reference = 0; reference < refs.length; reference++) {
//					String ref = refs[reference];
//
//					if (ref.equals("before")) {
//						filterimputationquality = false;
//					} else {
//						filterimputationquality = true;
//					}
//
////			if (d == 0) {
////				covariatefile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/T1D-covarmerged.txtmergedCovariatesUK.txt";
////				vcf = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/15-EverythingImputedAgain/T1D/" + ref + "/merged-filtered.Chr" + chr + ".vcf.gz";
////				if (ref.equals("before")) {
////					vcf = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip2/T1D/merged/merged-Chr" + chr + ".vcf.gz";
////				}
////				diseasestatus = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/T1D-covarmerged.txtmergeddisease.txt";
////				if (pseudo) {
////					famfile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/T1D.fam";
////					out = outputdir + "T1D/" + ref + "-pseudo/";
////				} else {
////					out = outputdir + "T1D/" + ref + "-UK/";
////				}
////
////			} else {
////				covariatefile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/RA-covarmerged.txtmergedCovariates.txt";
////				vcf = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/15-EverythingImputedAgain/RA/" + ref + "/merged-filtered.Chr" + chr + ".vcf.gz";
////				diseasestatus = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/RA-covarmerged.txtmergeddisease.txt";
////				if (ref.equals("before")) {
////					vcf = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip2/RA/merged/merged-Chr" + chr + ".vcf.gz";
////				}
////				out = outputdir + "RA/" + ref + "/";
////			}
//
//					if (d == 0) {
//						covariatefile = "/medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/T1D-covarmerged.txtmergedCovariates-withPseudos.txt";
//						vcf = "/medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/T1D/imputed-merged/" + ref + "/merged-Chr" + chr + ".vcf.gz";
//						samplestoexclude = "/medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/T1D-parentsToExclude.txt";
//						diseasestatus = "/medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/T1D-diseaseStatus.txt";
//						famfile = "/medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/T1D.fam";
//						out = "/medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/T1D/gwas-iter/" + ref + "/";
//
//					} else {
//						covariatefile = "/medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/RA-covarmerged.txtmergedCovariates.txt";
//						vcf = "/medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/RA/imputed-merged/" + ref + "/merged-Chr" + chr + ".vcf.gz";
//						famfile = null;
//						samplestoexclude = null;
//						diseasestatus = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/RA-covarmerged.txtmergeddisease.txt";
//						out = "/medpop/srlab/hwestra/fr-reseq/2015-09-25-WithX/RA/gwas-iter/" + ref + "/";
//					}
//					Gpio.createDir(out);
//					out += "Chr" + chr + "-";
//
//					HashSet<String> snpLimit = null;
//
//					if (chr != 8 && chr != 3) {
//						if (Gpio.exists(vcf)) {
//							LogitTestR testObj = new LogitTestR(vcf,
//									out,
//									diseasestatus,
//									covariatefile,
//									covariatestoinclude,
//									snpLimit,
//									samplestoexclude,
//									famfile,
//									filterimputationquality,
//									imputationqualitythreshold,
//									mafthreshold, submit);
//							testObj.skipMakingPseudoControls(true);
//
//
//							testObj.setRunIterative(true, getBedRegionsForChr(chr, bedregions));
//							pool.submit(testObj);
//
//							submit++;
//						}
//
//					}
//
////				if (d == 0) {
////					for (int set = 0; set < 2; set++) {
////						famfile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/T1D.fam";
////						if (set == 0) {
////							covariatefile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/T1D-covarmerged.txtmergedCovariatesUS.txt";
////							out = outputdir + "T1D/" + ref + "-US-pseudo/";
////						} else {
////							covariatefile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/T1D-covarmerged.txtmergedCovariatesEUR.txt";
////							out = outputdir + "T1D/" + ref + "-EUR-pseudo/";
////						}
////
////						vcf = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/15-EverythingImputedAgain/T1D/" + ref + "/merged-filtered.Chr" + chr + ".vcf.gz";
////						if (ref.equals("before")) {
////							vcf = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/ImmunoChip2/T1D/merged/merged-Chr" + chr + ".vcf.gz";
////						}
////						diseasestatus = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/2015-05-30-Eur-WoMendelianErrors/T1D-covarmerged.txtmergeddisease.txt";
////
////						if (Gpio.exists(vcf)) {
////							LogitTestR testObj = new LogitTestR(vcf,
////									out,
////									diseasestatus,
////									covariatefile,
////									covariatestoinclude,
////									snpLimit,
////									samplestoexclude,
////									famfile,
////									filterimputationquality,
////									imputationqualitythreshold,
////									mafthreshold, submit);
////							pool.submit(testObj);
////							submit++;
////						}
////					}
////				}
//
//				}
//
//			}
//		}
////		}
//
//		int returned = 0;
//		while (returned < submit) {
//
//			try {
//				Boolean result = pool.take().get();
//				if (result) {
//					returned++;
//				} else {
//
//					System.exit(-1);
//				}
//
//				System.out.println(returned + " / " + submit + " returned");
//			} catch (Exception e) {
//				e.printStackTrace();
//			}
//
//		}
//
//		threadPool.shutdown();
//
//
//	}

	private ArrayList<Feature> getBedRegionsForChr(int chr, ArrayList<Feature> bedregions) {
		ArrayList<Feature> output = new ArrayList<>();
		Chromosome chrom = Chromosome.parseChr("" + chr);
		for (Feature f : bedregions) {
			if (f.getChromosome().equals(chrom)) {
				output.add(f);
			}
		}
		return output;
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

//			t.concatVCF(vcfConcat, vcfsort, files, mergedout, mergedoutsorted, bashfilename); // concatenate and sort

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

//	public void regionR2Plots(String regionsFile, String outputdir, String before, String after, String sequencingregionfile, double maf) throws IOException {
//
//		try {
//
//
//			String gtf = "";
//			Gpio.createDir(outputdir);
//
//			TextFile tf = new TextFile(regionsFile, TextFile.R);
//
//			String[] elems = tf.readLineElems(TextFile.tab);
//			while (elems != null) {
//				if (elems.length >= 3) {
//					Feature region = new Feature();
//					region.setChromosome(Chromosome.parseChr(elems[0]));
//					region.setStart(Integer.parseInt(elems[1]));
//					region.setStop(Integer.parseInt(elems[2]));
//
//					int width = region.getStop() - region.getStart();
//
//					if (width > 1000000) {
//						width /= 1000;
//					}
//					System.out.println("region width: " + region.toString() + " -- " + width);
//
//					VariantPlot plot = new VariantPlot(outputdir + region.toString() + "-" + maf + ".pdf", width + 200, 500);
//					plot.setMargin(100);
////					plot.plotImputationRSquared(before, after, sequencingregionfile, region, maf, gtf);
//
//
//				}
//				elems = tf.readLineElems(TextFile.tab);
//			}
//			tf.close();
//
//
//		} catch (Exception e1) {
//			e1.printStackTrace();
//		}
//	}


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

//	public void prepareImputationPanels() throws IOException {
//
//		GenotypeTools t = new GenotypeTools();
//
//		VCFFunctions v = new VCFFunctions();
//		PedAndMapFunctions p = new PedAndMapFunctions();
//
//		// broad
////		String outputPath = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/reference/";
////
////		// refilter the reference VCF.
////		String startVCF = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/merged-ICIds-MixupsFixed.vcf";
////		String regionFile = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/allLoci.bed";
////		String mergedhg19immunochipPed = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/immunochipsequencedsamples/merge";
////		String mergedhg19immunochipMap = mergedhg19immunochipPed + ".map";
////		String mergedhg19immunochipFAM = mergedhg19immunochipPed + ".fam";
////		String plink = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/plink/plink";
////		String beagle = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/beagle.r1399.jar";
////		String merged1kg = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/1kg.phase3.v5-filtered.merged.vcf.gz";
////		String vcfConcat = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/vcftools/vcftools_0.1.12b/bin/vcf-concat";
////		String vcfsort = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/vcftools/vcftools_0.1.12b/bin/vcf-sort";
////		String mergeFam = "/medpop/srlab/hwestra/fr-reseq/2015-05-19-Imputation/mergefam.txt";
//
//		// # local
//		String outputPath = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-11-29-AllChr-EurOnly-WithX/";
//
//		// refilter the reference VCF.
//		String startVCF = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-11-MixupsFixed/merged-ICIds-MixupsFixed.vcf";
//		String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";
//		String mergedhg19immunochipPed = "/Data/ImmunoChip/SequencingProject/2015-11-29/Merged/merge2-sequencedSamples";
//		String mergedhg19immunochipMap = mergedhg19immunochipPed + ".map";
//		String startVCFFAM = mergedhg19immunochipPed + ".fam";
//		String plink = "/Data/Tools/plink-1.07-mac-intel/plink1.9";
//		String beagle = "/Data/Tools/beagle/beagle.r1399.jar";
//		String merged1kg = "/Data/Ref/BeagleRef/1kg.phase3.v5-filtered.eur.merged-withX-sorted.vcf.gz"; // EUR ONLY!
//		String kgdir = "/Data/Ref/BeagleRef/";
//
//		String merged1kgFam = "/Data/Ref/BeagleRef/20130606_g1k.ped.fam";
//
//		String shapeit = "/Data/Tools/shapeit/shapeit";
//
//		String vcfConcat = "/Data/Tools/vcftools/bin/vcf-concat";
//		String vcfsort = "/Data/Tools/vcftools/bin/vcf-sort";
//		String mergeFam = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/mergefam.txt";
//
//		String shapeitgeneticmap = "/Data/Ref/ImputeRef/geneticmaps/genetic_map_CHROMOSOME_combined_b37.txt";
//		String shapeitref = "/Data/Ref/ImputeRef/haplotypes/1000GP_Phase3_CHROMOSOME.hap.gz";
//		String shapeitreflegend = "/Data/Ref/ImputeRef/haplotypes/1000GP_Phase3_CHROMOSOME.legend.gz";
//		String shapeitrefsample = "/Data/Ref/ImputeRef/1000GP_Phase3.sample";
//		String grouplistfile = "/Data/Ref/ImputeRef/group.list";
//
//
//		boolean skipimpute = false;
//		boolean onlymerge = false;
//		boolean linux = false;
//		boolean phasingWithBeagle = false;
//		boolean phasingWithShapeit = true;
//		boolean skipPhasing = false;
//
//		try {
//			// filter bad variants
//			String filteredVCFOut = outputPath + "1-SeqVariantFilter/";
//			Gpio.createDir(filteredVCFOut);
//
//			v.filterLowFrequencyVariants(startVCF, filteredVCFOut, startVCFFAM, true, 10, 30, 0.90, 5);
//			String mendelianFilteredVCF = filteredVCFOut + "filtered-mendelianerrorsremoved.vcf";
//
//			v.filterMendelianErrors(filteredVCFOut + "filtered.vcf.gz", mergeFam, mendelianFilteredVCF, filteredVCFOut);
//			v.summarizeVCF(mendelianFilteredVCF, mendelianFilteredVCF + "-summary.txt");
//
////			mendelianFilteredVCF = filteredVCFOut + "filtered.vcf.gz";
////			v.splitMultipleAllelicVariants(filteredVCFOut + "filtered.vcf", mendelianFilteredVCF);
//
//			// compare against IC genotypes
//			//	t.compareVCFGenotypesToPedAndMap(filteredVCFOut + "filtered.vcf.gz", mergedhg19immunochipPed, filteredVCFOut, true);
//
//			// filter for CD28 region
//			String filteredVariantOut = outputPath + "2-SeqVariantRegionFilter/";
//			Gpio.createDir(filteredVariantOut);
//
//			v.filterVCFForBedRegions(mendelianFilteredVCF, filteredVariantOut + "filtered.vcf.gz", regionFile);
////			v.summarizeVCF(filteredVariantOut + "filtered.vcf.gz", filteredVariantOut + "summary.txt");
//
//			// filter immunochip for sequenced regions and samples
//			String regionFilteredMapOut = outputPath + "3-ICRegionFiltered/";
//			Gpio.createDir(regionFilteredMapOut);
//			String variantsToKeep = regionFilteredMapOut + "variantsOverlappingSequencedRegions.txt";
//			p.filterMap(mergedhg19immunochipMap, regionFile, variantsToKeep);
//			String variantsToExclude = regionFilteredMapOut + "ICvariantsOverlappingSequencedVariants.txt";
//			p.filterMapForVCFVariants(mergedhg19immunochipMap, filteredVariantOut + "filtered.vcf.gz", variantsToExclude);
//
//
//			String regionfilteredpedOut = regionFilteredMapOut + "filtered";
//			ProcessBuilder pb = new ProcessBuilder(plink,
//					"--extract", variantsToKeep,
//					"--file", mergedhg19immunochipPed,
//					"--exclude", variantsToExclude,
//					"--recode",
//					"--out", regionfilteredpedOut,
//					"--keep", mergeFam
//			);
//			t.run(pb);
//
//
//			// merge with immunochip
//			String ICAndSeqVariantMerged = outputPath + "4-ICAndSeqVariantMerged/";
//			Gpio.createDir(ICAndSeqVariantMerged);
//
//
//			v.mergeWithPed(regionfilteredpedOut, filteredVariantOut + "filtered.vcf.gz", ICAndSeqVariantMerged + "merged.vcf.gz");
//			v.summarizeVCF(ICAndSeqVariantMerged + "merged.vcf.gz", ICAndSeqVariantMerged + "merged-summary.txt");
//
//			// filter 1000 genomes for sequenced regions
////			String kgSeqRegions = outputPath + "5-1KG-SequencedRegions/";
////			Gpio.createDir(kgSeqRegions);
////			System.out.println("Filtering VCF for regions: " + merged1kg);
////			v.filterVCFForBedRegions(merged1kg, kgSeqRegions + "filtered.vcf.gz", regionFile);
////
////			String kgSeqRegionsHigherfreq = outputPath + "5-1KG-SequencedRegions-HighFreq/";
////			Gpio.createDir(kgSeqRegionsHigherfreq);
////			v.filterLowFrequencyVariants(kgSeqRegions + "filtered.vcf.gz", kgSeqRegionsHigherfreq, merged1kgFam, true, 0, 0, 0, 1);
////			System.out.println("Summarizing: " + kgSeqRegionsHigherfreq + "filtered.vcf.gz");
////			v.summarizeVCF(kgSeqRegionsHigherfreq + "filtered.vcf.gz", kgSeqRegionsHigherfreq + "filtered-summary.txt");
////			kgSeqRegions = kgSeqRegionsHigherfreq;
//
//
//			// compare reference and 1kg
//
////			// split per chromosome...
////			System.out.println("Splitting ref per chr");
////			v.splitPerChromosome(kgSeqRegions + "filtered.vcf.gz", kgSeqRegions + "filtered");
////			System.out.println("Splitting test per chr");
//			v.splitPerChromosome(ICAndSeqVariantMerged + "merged.vcf.gz", ICAndSeqVariantMerged + "merged");
//
//
//			for (Chromosome chr : Chromosome.values()) {
//
//				System.out.println(chr.getName() + " imputing..");
//
//				// filter 1000 genomes for sequenced regions
//				String kgSeqRegions = outputPath + "5-1KG-SequencedRegions/";
//				Gpio.createDir(kgSeqRegions);
//				String kgChrName = kgdir + chr.getName().toLowerCase() + ".1kg.phase3.v5.eur.vcf.gz";
//				String kgChrFilteredName = kgdir + chr.getName().toLowerCase() + ".filtered.1kg.phase3.v5.eur.vcf.gz";
//				String kgSeqRegionsHigherfreq = outputPath + "5-1KG-SequencedRegions-HighFreq/";
//				Gpio.createDir(kgSeqRegionsHigherfreq);
//				String kgChrHighFreq = kgSeqRegionsHigherfreq + chr.getName().toLowerCase() + "-";
//				if (chr.equals(Chromosome.X)) {
//					kgChrName = kgdir + "chrX.1kg.phase3.v5.eur.vcf.gz";
//					kgChrFilteredName = kgdir + "chrX.filtered.1kg.phase3.v5.eur.vcf.gz";
//					kgChrHighFreq = kgSeqRegionsHigherfreq + "chrX-";
//				}
//				System.out.println("Filtering VCF for regions: " + kgChrName);
//				v.filterVCFForBedRegions(kgChrName, kgChrFilteredName, regionFile);
//
//
//				v.filterLowFrequencyVariants(kgChrFilteredName, kgChrHighFreq, merged1kgFam, true, 0, 0, 0, 1);
//				kgSeqRegions = kgSeqRegionsHigherfreq;
//
//
//				String refVCFChr = kgChrHighFreq + "filtered.vcf.gz";//kgSeqRegions + "filtered-" + chr.getName() + ".vcf.gz";
//				System.out.println("Result: " + refVCFChr);
//
//				String testVCFChr = ICAndSeqVariantMerged + "merged-" + chr.getName() + ".vcf.gz";
//				if (Gpio.exists(refVCFChr) && Gpio.exists(testVCFChr)) {
//
//					String matchedPanelsOut = outputPath + "6-PanelsMatched/";
//					Gpio.createDir(matchedPanelsOut);
//
//					v.mergeAndIntersectVCFVariants(
//							refVCFChr,
//							testVCFChr,
//							matchedPanelsOut + "1kg-matched-" + chr.getName() + ".vcf.gz",
//							matchedPanelsOut + "seq-matched-" + chr.getName() + ".vcf.gz",
//							matchedPanelsOut + "1kg-seq-merged-" + chr.getName() + ".vcf.gz",
//							"/",
//							matchedPanelsOut + "mergelog-" + chr.getName() + ".txt",
//							true);
//
//					// sort vcfs
//					t.sortVCF(linux, vcfsort, matchedPanelsOut + "1kg-matched-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "1kg-matched-sorted-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "sort.sh");
//					t.sortVCF(linux, vcfsort, matchedPanelsOut + "seq-matched-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "seq-matched-sorted-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "sort.sh");
//
//					// phase sequencing data
//					if (phasingWithBeagle) {
//						pb = new ProcessBuilder("java",
//								"-Xmx6g",
//								"-jar", beagle,
//								"ped=" + mergeFam,
//								"nthreads=32",
//								"gt=" + matchedPanelsOut + "seq-matched-sorted-" + chr.getName() + ".vcf.gz",
//								"out=" + matchedPanelsOut + "seqpanel-matched-sorted-phased-" + chr.getName() + ""
//						);
//
//						t.run(pb);
//						t.sortVCF(linux, vcfsort, matchedPanelsOut + "seqpanel-matched-sorted-phased-" + chr.getName() + ".vcf.gz",
//								matchedPanelsOut + "seqpanel-matched-sorted-phased-sorted-" + chr.getName() + ".vcf.gz", matchedPanelsOut + "sort.sh");
//
//						if (!skipimpute) {
//							// impute 1kg into reference panel
//							String imputedPanelsOut = outputPath + "7-PanelsImputed/";
//							Gpio.createDir(imputedPanelsOut);
//
//							// phased data
//							pb = new ProcessBuilder("java",
//									"-Xmx6g",
//									"-jar", beagle,
//									"ref=" + matchedPanelsOut + "1kg-matched-sorted-" + chr.getName() + ".vcf.gz",
//									//"ped=" + mergeFam,
//									"usephase=true", "burnin-its=0", "phase-its=0",
//									"gt=" + matchedPanelsOut + "seqpanel-matched-sorted-phased-sorted-" + chr.getName() + ".vcf.gz",
//									"nthreads=32",
//									"out=" + imputedPanelsOut + "seqpanel-1kgimputed-" + chr.getName() + ""
//							);
//							t.run(pb);
//
//
//							pb = new ProcessBuilder("java",
//									"-Xmx6g",
//									"-jar", beagle,
//									"ref=" + matchedPanelsOut + "1kg-matched-sorted-" + chr.getName() + ".vcf.gz",
//									"ped=" + mergeFam,
//									"gt=" + matchedPanelsOut + "seq-matched-sorted-" + chr.getName() + ".vcf.gz",
//									"nthreads=32",
//									"out=" + imputedPanelsOut + "seqpanel-unphased-1kgimputed-" + chr.getName() + ""
//							);
//							t.run(pb);
//
//							t.sortVCF(linux, vcfsort, imputedPanelsOut + "seqpanel-1kgimputed-" + chr.getName() + ".vcf.gz",
//									imputedPanelsOut + "seqpanel-1kgimputed-sorted-" + chr.getName() + ".vcf.gz", imputedPanelsOut + "sort.sh");
//
//
//							// impute reference panel into 1kg
//							pb = new ProcessBuilder("java",
//									"-Xmx6g",
//									"-jar", beagle,
//									"ref=" + matchedPanelsOut + "seqpanel-matched-sorted-phased-sorted-" + chr.getName() + ".vcf.gz",
//									// "ped=" + mergeFam,
//									"usephase=true", "burnin-its=0", "phase-its=0",
//									"gt=" + matchedPanelsOut + "1kg-matched-sorted-" + chr.getName() + ".vcf.gz",
//									"nthreads=32",
//									"out=" + imputedPanelsOut + "1kg-seqpanelimputed-" + chr.getName()
//							);
//							t.run(pb);
//							t.sortVCF(linux, vcfsort, imputedPanelsOut + "1kg-seqpanelimputed-" + chr.getName() + ".vcf.gz",
//									imputedPanelsOut + "1kg-seqpanelimputed-sorted-" + chr.getName() + ".vcf.gz", imputedPanelsOut + "sort.sh");
//
//
//							// merge with original phased datasets..
//							String imputedPanelsMergedOut = outputPath + "8-PanelsImputedMerged/";
//							Gpio.createDir(imputedPanelsMergedOut);
//							String unimputed1kg = matchedPanelsOut + "1kg-matched-sorted-" + chr.getName() + ".vcf.gz";
//							String imputed1kg = imputedPanelsOut + "1kg-seqpanelimputed-sorted-" + chr.getName() + ".vcf.gz";
//							String unimputed1kgout = imputedPanelsMergedOut + "1kg-imputedVariantsFiltered-" + chr.getName() + ".vcf.gz";
//							v.filterVCFVariants(unimputed1kg, imputed1kg, unimputed1kgout); // remove imputed variants from the original vcf if there is any overlap
//
//							String files = unimputed1kgout + " " + imputed1kg;
//							String mergedout = imputedPanelsMergedOut + "1kg-seqpanelimputed-mergedWithUnimputed-" + chr.getName() + ".vcf.gz";
//							String mergedoutsorted = imputedPanelsMergedOut + "1kg-seqpanelimputed-seqpanelimputed-mergedWithUnimputed-sorted-" + chr.getName() + ".vcf.gz";
//							String bashfilename = imputedPanelsMergedOut + "sortmerge.sh";
//							System.out.println("concat 1");
//							System.out.println(files);
//							System.out.println(mergedout);
//							System.out.println(mergedoutsorted);
//							t.concatVCF(vcfConcat, vcfsort, files, mergedout, mergedoutsorted, bashfilename); // concatenate vcf files
//
//							// repeat on the sequenced data..
//							t.sortVCF(linux, vcfsort, imputedPanelsOut + "seqpanel-1kgimputed-" + chr.getName() + ".vcf.gz",
//									imputedPanelsOut + "seqpanel-1kgimputed-sorted-" + chr.getName() + ".vcf.gz", imputedPanelsOut + "sort.sh");
//
//							String unimputedseq = matchedPanelsOut + "seqpanel-matched-sorted-phased-sorted-" + chr.getName() + ".vcf.gz";
//							String imputedseq = imputedPanelsOut + "seqpanel-1kgimputed-" + chr.getName() + ".vcf.gz";
//							String unimputedseqout = imputedPanelsMergedOut + "seqpanel-imputedVariantsFiltered-" + chr.getName() + ".vcf.gz";
//							v.filterVCFVariants(unimputedseq, imputedseq, unimputedseqout); // remove imputed variants from original vcf if there is any overlap
//
//							files = unimputedseqout + " " + imputedseq;
//							mergedout = imputedPanelsMergedOut + "seqpanel-1kgimputed-mergedWithUnimputed-" + chr.getName() + ".vcf.gz";
//							mergedoutsorted = imputedPanelsMergedOut + "seqpanel-1kgimputed-1kgimputed-mergedWithUnimputed-sorted-" + chr.getName() + ".vcf.gz";
//							bashfilename = imputedPanelsMergedOut + "sortmerge.sh";
//							System.out.println("concat 2");
//							System.out.println(files);
//							System.out.println(mergedout);
//							System.out.println(mergedoutsorted);
//							t.concatVCF(vcfConcat, vcfsort, files, mergedout, mergedoutsorted, bashfilename); // concatenate
//
//							// bgzip and index
//							System.out.println("bgzip and index");
//							v.replaceHeader(imputedPanelsMergedOut + "seqpanel-1kgimputed-1kgimputed-mergedWithUnimputed-sorted-" + chr.getName() + ".vcf.gz",
//									imputed1kg,
//									imputedPanelsMergedOut + "seqpanel-1kgimputed-1kgimputed-mergedWithUnimputed-sorted-header-" + chr.getName() + ".vcf.gz");
//
//							// merge
//							System.out.println("merge");
//							String finalImputationPanelsOut = outputPath + "9-ImputationPanels/";
//							Gpio.createDir(finalImputationPanelsOut);
//
//							v.mergeAndIntersectVCFVariants(
//									imputedPanelsMergedOut + "seqpanel-1kgimputed-1kgimputed-mergedWithUnimputed-sorted-header-" + chr.getName() + ".vcf.gz",
//									imputedPanelsMergedOut + "seqpanel-1kgimputed-1kgimputed-mergedWithUnimputed-sorted-" + chr.getName() + ".vcf.gz",
//									finalImputationPanelsOut + "1kg-" + chr.getName() + ".vcf.gz",
//									finalImputationPanelsOut + "seq-" + chr.getName() + ".vcf.gz",
//									finalImputationPanelsOut + "1kg-seq-merged-" + chr.getName() + ".vcf.gz",
//									"|",
//									finalImputationPanelsOut + "mergelog-" + chr.getName() + ".txt",
//									true);
//
//
//							v.summarizeVCF(finalImputationPanelsOut + "1kg-seq-merged-" + chr.getName() + ".vcf.gz", finalImputationPanelsOut + "1kg-seq-merged-summarystats-" + chr.getName() + ".txt");
//						}
//					}
//
//					if (phasingWithShapeit) {
//						// shapeit
//						// convert to ped/map/fam --> use plink
//
//						String phasetmpdir = matchedPanelsOut + "shapeit/";
//						System.out.println("phase tmp dir: " + phasetmpdir);
//
//
//						Gpio.createDir(phasetmpdir);
//						pb = new ProcessBuilder(plink,
//								"--vcf", matchedPanelsOut + "seq-matched-sorted-" + chr.getName() + ".vcf.gz",
//								"--out", phasetmpdir + "seq-matched-sorted-" + chr.getName(),
//								"--recode"
//						);
//						t.run(pb);
//
//
//						if (!skipPhasing) {
//
//							if (chr.equals(Chromosome.X)) {
//
//							} else {
//
//
////							// need to check the 'matched' data first
////							// remove non-matching variants
////
////
////							// update family info in ped file
//								p.updatePedWithFAM(phasetmpdir + "seq-matched-sorted-" + chr.getName() + ".ped", mergeFam);
////
//
//								String chrMap = shapeitgeneticmap.replaceAll("CHROMOSOME", chr.getName().toLowerCase());
//								String chrRef = shapeitref.replaceAll("CHROMOSOME", chr.getName().toLowerCase());
//								String chrLegend = shapeitreflegend.replaceAll("CHROMOSOME", chr.getName().toLowerCase());
////
////							/*
////							shapeit -check \
////        -B gwas \
////        -M genetic_map.txt \
////        --input-ref reference.haplotypes.gz reference.legend.gz reference.sample \
////        --output-log gwas.alignments
////							 */
////
////							// convert ped/map to gen/sample
////							// gtool -P --ped myPlinkTextData.ped --map myPlinkTextData.map --og myGtoolTextData.gen --os myGtoolTextData.sample
////							p.convertToGenSample(phasetmpdir + "seq-matched-sorted-" + chr.getName());
//////
//////							System.exit(-1);
////
//
//								// phase without a reference
//								pb = new ProcessBuilder(shapeit,
//										"-P", phasetmpdir + "seq-matched-sorted-" + chr.getName(),
//										"-M", chrMap,
//										"-T", "4",
//										"--duohmm",
//										"-O", phasetmpdir + "seq-matched-sorted-shapeitRound1-" + chr.getName()
//
//								);
//								t.run(pb);
//
//								// remove variants that phasing with reference would choke on
//								pb = new ProcessBuilder(shapeit, "-check",
//										"-P", phasetmpdir + "seq-matched-sorted-" + chr.getName(),
//										"-M", chrMap,
//										"--input-ref", chrRef, chrLegend, shapeitrefsample,
//										"--include-grp", grouplistfile,
//										"-T", "4",
//										"--output-log", phasetmpdir + "seq-matched-sorted-shapeit.alignments." + chr.getName()
//								);
//								t.run(pb);
////
//								// phase using the 1kg as a reference
//								pb = new ProcessBuilder(shapeit,
//										"-P", phasetmpdir + "seq-matched-sorted-" + chr.getName(),
//										"-M", chrMap,
//										"--input-ref", chrRef, chrLegend, shapeitrefsample,
//										"--include-grp", grouplistfile,
//										"-T", "4",
//										"--duohmm",
//										"--exclude-snp", phasetmpdir + "seq-matched-sorted-shapeit.alignments." + chr.getName() + ".snp.strand.exclude",
//										"-O", phasetmpdir + "seq-matched-sorted-shapeitRound2-" + chr.getName()
//
//								);
//								t.run(pb);
//
////							// reintroduce variants unique to dataset: put them back into the gen/sample file
//								reintroduceVariants(phasetmpdir + "seq-matched-sorted-shapeit.alignments." + chr.getName() + ".snp.strand.exclude",
//										phasetmpdir + "seq-matched-sorted-shapeitRound2-" + chr.getName(),
//										phasetmpdir + "seq-matched-sorted-shapeitRound1-" + chr.getName(),
//										phasetmpdir + "seq-matched-sorted-shapeitRound2-WithNovelVariants-" + chr.getName()
//								);
//
////							 System.exit(-1);
//
//								// convert back to VCF
//
//							}
//						}
//
//						if (Gpio.exists(phasetmpdir + "seq-matched-sorted-shapeitRound2-WithNovelVariants-" + chr.getName() + ".haps")) {
//							String phasedSeqData = covertHapsSampleToVCF(phasetmpdir + "seq-matched-sorted-shapeitRound2-WithNovelVariants-" + chr.getName(),
//									phasetmpdir + "seq-matched-sorted-shapeitRound2-WithNovelVariants-" + chr.getName(), linux, vcfsort);
//
//							if (!skipimpute) {
//// impute 1kg into reference panel
//								String imputedPanelsOut = outputPath + "7-PanelsImputed/";
//								Gpio.createDir(imputedPanelsOut);
//
//
//								// impute reference panel into 1kg
//								pb = new ProcessBuilder("java",
//										"-Xmx7g",
//										"-jar", beagle,
//										"ref=" + phasedSeqData,
//										// "ped=" + mergeFam,
//										"usephase=true", "burnin-its=0", "phase-its=0",
//										"gt=" + matchedPanelsOut + "1kg-matched-sorted-" + chr.getName() + ".vcf.gz",
//										"nthreads=32",
//										"out=" + imputedPanelsOut + "1kg-seqpanelimputed-" + chr.getName()
//								);
//								t.run(pb);
//								t.sortVCF(linux, vcfsort, imputedPanelsOut + "1kg-seqpanelimputed-" + chr.getName() + ".vcf.gz",
//										imputedPanelsOut + "1kg-seqpanelimputed-sorted-" + chr.getName() + ".vcf.gz", imputedPanelsOut + "sort.sh");
//
//								String imputedPanelsMergedOut = outputPath + "8-PanelsImputedMerged/";
//								Gpio.createDir(imputedPanelsMergedOut);
//								String unimputed1kg = matchedPanelsOut + "1kg-matched-sorted-" + chr.getName() + ".vcf.gz";
//								String imputed1kg = imputedPanelsOut + "1kg-seqpanelimputed-sorted-" + chr.getName() + ".vcf.gz";
//								String unimputed1kgout = imputedPanelsMergedOut + "1kg-imputedVariantsFiltered-" + chr.getName() + ".vcf.gz";
//								v.filterVCFVariants(unimputed1kg, imputed1kg, unimputed1kgout); // remove imputed variants from the original vcf if there is any overlap
//
//								String files = unimputed1kgout + " " + imputed1kg;
//								String mergedout = imputedPanelsMergedOut + "1kg-seqpanelimputed-mergedWithUnimputed-" + chr.getName() + ".vcf.gz";
//								String mergedoutsorted = imputedPanelsMergedOut + "1kg-seqpanelimputed-seqpanelimputed-mergedWithUnimputed-sorted-" + chr.getName() + ".vcf.gz";
//								String bashfilename = imputedPanelsMergedOut + "sortmerge.sh";
//								System.out.println("concat 1");
//								System.out.println(files);
//								System.out.println(mergedout);
//								System.out.println(mergedoutsorted);
//								t.concatVCF(vcfConcat, vcfsort, files, mergedout, mergedoutsorted, bashfilename); // concatenate vcf files
//
//
//								// repeat on the sequenced data..
//								// phased data
//								pb = new ProcessBuilder("java",
//										"-Xmx7g",
//										"-jar", beagle,
//										"ref=" + matchedPanelsOut + "1kg-matched-sorted-" + chr.getName() + ".vcf.gz",
//										//"ped=" + mergeFam,
//										"usephase=true", "burnin-its=0", "phase-its=0",
//										"gt=" + phasedSeqData,
//										"nthreads=32",
//										"out=" + imputedPanelsOut + "seqpanel-1kgimputed-" + chr.getName() + ""
//								);
//								t.run(pb);
//								t.sortVCF(linux, vcfsort, imputedPanelsOut + "seqpanel-1kgimputed-" + chr.getName() + ".vcf.gz",
//										imputedPanelsOut + "seqpanel-1kgimputed-sorted-" + chr.getName() + ".vcf.gz", imputedPanelsOut + "sort.sh");
//
//								String unimputedseq = phasedSeqData;
//								String imputedseq = imputedPanelsOut + "seqpanel-1kgimputed-" + chr.getName() + ".vcf.gz";
//								String unimputedseqout = imputedPanelsMergedOut + "seqpanel-imputedVariantsFiltered-" + chr.getName() + ".vcf.gz";
//								v.filterVCFVariants(unimputedseq, imputedseq, unimputedseqout); // remove imputed variants from original vcf if there is any overlap
//
//								files = unimputedseqout + " " + imputedseq;
//								mergedout = imputedPanelsMergedOut + "seqpanel-1kgimputed-mergedWithUnimputed-" + chr.getName() + ".vcf.gz";
//								mergedoutsorted = imputedPanelsMergedOut + "seqpanel-1kgimputed-1kgimputed-mergedWithUnimputed-sorted-" + chr.getName() + ".vcf.gz";
//								bashfilename = imputedPanelsMergedOut + "sortmerge.sh";
//								System.out.println("concat 2");
//								System.out.println(files);
//								System.out.println(mergedout);
//								System.out.println(mergedoutsorted);
//								t.concatVCF(vcfConcat, vcfsort, files, mergedout, mergedoutsorted, bashfilename); // concatenate
//
//								// bgzip and index
//								System.out.println("bgzip and index");
//								v.replaceHeader(imputedPanelsMergedOut + "seqpanel-1kgimputed-1kgimputed-mergedWithUnimputed-sorted-" + chr.getName() + ".vcf.gz",
//										imputed1kg,
//										imputedPanelsMergedOut + "seqpanel-1kgimputed-1kgimputed-mergedWithUnimputed-sorted-header-" + chr.getName() + ".vcf.gz");
//
//								// merge
//								System.out.println("merge");
//								String finalImputationPanelsOut = outputPath + "9-ImputationPanels/";
//								Gpio.createDir(finalImputationPanelsOut);
//
//								v.mergeAndIntersectVCFVariants(
//										imputedPanelsMergedOut + "1kg-seqpanelimputed-seqpanelimputed-mergedWithUnimputed-sorted-" + chr.getName() + ".vcf.gz",
//										imputedPanelsMergedOut + "seqpanel-1kgimputed-1kgimputed-mergedWithUnimputed-sorted-" + chr.getName() + ".vcf.gz",
//										finalImputationPanelsOut + "1kg-" + chr.getName() + ".vcf.gz",
//										finalImputationPanelsOut + "seq-" + chr.getName() + ".vcf.gz",
//										finalImputationPanelsOut + "1kg-seq-merged-" + chr.getName() + ".vcf.gz",
//										"|",
//										finalImputationPanelsOut + "mergelog-" + chr.getName() + ".txt",
//										true);
//
//
//								v.summarizeVCF(finalImputationPanelsOut + "1kg-seq-merged-" + chr.getName() + ".vcf.gz", finalImputationPanelsOut + "1kg-seq-merged-summarystats-" + chr.getName() + ".txt");
//							}
//
//							if (onlymerge) {
//								String imputedPanelsMergedOut = outputPath + "8-PanelsImputedMerged/";
//								String finalImputationPanelsOut = outputPath + "9-ImputationPanels/";
//								VCFFunctions func = new VCFFunctions();
//								func.mergeAndIntersectVCFVariants(
//										imputedPanelsMergedOut + "1kg-seqpanelimputed-seqpanelimputed-mergedWithUnimputed-sorted-" + chr.getName() + ".vcf.gz",
//										imputedPanelsMergedOut + "seqpanel-1kgimputed-1kgimputed-mergedWithUnimputed-sorted-" + chr.getName() + ".vcf.gz",
//										finalImputationPanelsOut + "1kg-" + chr.getName() + ".vcf.gz",
//										finalImputationPanelsOut + "seq-" + chr.getName() + ".vcf.gz",
//										finalImputationPanelsOut + "1kg-seq-merged-" + chr.getName() + ".vcf.gz",
//										"|",
//										finalImputationPanelsOut + "mergelog-" + chr.getName() + ".txt",
//										true);
//
//								v.filterVariantsWithWeirdEncoding(finalImputationPanelsOut + "seq-" + chr.getName() + ".vcf.gz");
//								v.filterVariantsWithWeirdEncoding(finalImputationPanelsOut + "1kg-" + chr.getName() + ".vcf.gz");
//								v.filterVariantsWithWeirdEncoding(finalImputationPanelsOut + "1kg-seq-merged-" + chr.getName() + ".vcf.gz");
//							}
//						}
//
//
//					}
//
//
//				}
//
//			}
//
//
//		} catch (Exception e) {
//
//			e.printStackTrace();
//		}
//
//	}


	/*
	reintroduceVariants(phasetmpdir + "seq-matched-sorted-shapeit.alignments." + chr.getName() + ".snp.strand.exclude",
									phasetmpdir + "seq-matched-sorted-shapeit-" + chr.getName()+".haps",
									phasetmpdir + "seq-matched-sorted-" + chr.getName(),
									phasetmpdir + "seq-matched-sorted-shapeit-withNovelVariants" + chr.getName()
									);
	 */
	private void reintroduceVariants(String excludefile, String phasedRound2File, String phasedRound1File, String outfile) throws IOException {

		// read list of SNPs that were excluded
		TextFile tf = new TextFile(excludefile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		HashSet<String> idsToAdd = new HashSet<String>();
		while (elems != null) {
/*
			type	pos	main_id	main_A	main_B	main_flippable	ref_id	ref_A	ref_B	ref_flippable
Missing	Missing	2487762	rs2227313	G	A	1	NA	NA	NA	1
*/
			if (elems.length > 3) {
				idsToAdd.add(elems[3] + ":" + elems[2]);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		TextFile output = new TextFile(outfile + ".haps", TextFile.W);
		TextFile hapsin = new TextFile(phasedRound2File + ".haps", TextFile.R);

		String ln = hapsin.readLine();
		while (ln != null) {
			output.writeln(ln);
			ln = hapsin.readLine();
		}
		hapsin.close();

		// read variants that were actually excluded
		// read map file first to identify the column to read.
		TextFile gen = new TextFile(phasedRound1File + ".haps", TextFile.R);


		// 1 rs10910095 2510755 G A 0 0
		String[] genelems = gen.readLineElems(Strings.whitespace);
		while (genelems != null) {

			String id = genelems[1];
			String position = genelems[3];

			if (idsToAdd.contains(id + ":" + position)) {
				// write to output
				output.writeln(Strings.concat(genelems, Strings.whitespace));
			}

			genelems = gen.readLineElems(Strings.whitespace);
		}
		gen.close();
		output.close();

		System.out.println("copying: " + phasedRound2File + ".sample to " + outfile + ".sample");
		Gpio.copyFile(phasedRound2File + ".sample", outfile + ".sample");

	}


//	public void prepareGWASDatasets() {
//
//		try {
//
//			String removeTheseIndividuals = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/mergefam.txt";
//
//			String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";
//			String[] datasets = new String[]{
//					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_ES_QCgp",
//					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_NL_QCgp",
//					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_SE-E_QCgp",
//					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_SE-U_QCgp",
//					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_UK_QCgp",
//					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_US_QCgp",
//			};
//
//			String[] datasetOutput = new String[]{
//					"/Data/ImmunoChip/RA/AllRegions/ES/2015-09-23-IncX/",
//					"/Data/ImmunoChip/RA/AllRegions/NL/2015-09-23-IncX/",
//					"/Data/ImmunoChip/RA/AllRegions/SEE/2015-09-23-IncX/",
//					"/Data/ImmunoChip/RA/AllRegions/SEU/2015-09-23-IncX/",
//					"/Data/ImmunoChip/RA/AllRegions/UK/2015-09-23-IncX/",
//					"/Data/ImmunoChip/RA/AllRegions/US/2015-09-23-IncX/",
//			};
//
//			String[] refNames = new String[]{"1kg", "seq", "1kgseq"};
//			String[] refSets = new String[]{
//					"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/5-1KG-SequencedRegions/filtered.vcf",
//					"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/6-PanelsMatched/seqpanel-phased-sorted.vcf",
//					"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/9-PanelsMerged/merged.vcfsorted.vcf"};
//
//
//			for (int i = 0; i < datasets.length; i++) {
//				String outdir = datasetOutput[i];
//				Gpio.createDir(outdir);
//				liftOverAndFilterPlinkDataset(datasets[i], removeTheseIndividuals, outdir, regionFile);
//			}
//
//
//			datasets = new String[]{
//					"/Data/ImmunoChip/T1D/eur",
//					"/Data/ImmunoChip/T1D/uk"
//			};
//
//			datasetOutput = new String[]{
//					"/Data/ImmunoChip/T1D/AllRegions/eur/2015-09-23-IncX/",
//					"/Data/ImmunoChip/T1D/AllRegions/uk/2015-09-23-IncX/"
//			};
//
//			for (int i = 0; i < datasets.length; i++) {
//				String outdir = datasetOutput[i];
//				Gpio.createDir(outdir);
//				liftOverAndFilterPlinkDataset(datasets[i], removeTheseIndividuals, outdir, regionFile);
//			}
//
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
//
//	}

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
//
//	public void datasetmerger() throws IOException {
//		String[] datasetOutput = new String[]{
//				"/Data/ImmunoChip/RA/2015-09-23-IncX/ES/",
//				"/Data/ImmunoChip/RA/2015-09-23-IncX/NL/",
//				"/Data/ImmunoChip/RA/2015-09-23-IncX/SEE/",
//				"/Data/ImmunoChip/RA/2015-09-23-IncX/SEU/",
//				"/Data/ImmunoChip/RA/2015-09-23-IncX/UK/",
//				"/Data/ImmunoChip/RA/2015-09-23-IncX/US/",
//		};
//
//
//		String outputdir = "/Data/ImmunoChip/RA/2015-09-23-IncX/merged/";
//		mergeGWASDatasets(datasetOutput, outputdir);
//		datasetOutput = new String[]{
//				"/Data/ImmunoChip/T1D/2015-09-23-IncX/eur/",
//				"/Data/ImmunoChip/T1D/2015-09-23-IncX/uk/"
//		};
//		outputdir = "/Data/ImmunoChip/T1D/2015-09-23-IncX/merged/";
//		mergeGWASDatasets(datasetOutput, outputdir);
//	}

//	public void mergeGWASDatasets(String[] datasetOutput, String outputdir) throws IOException {
//
//
//		Gpio.createDir(outputdir);
//		for (int chr = 1; chr < 23; chr++) {
//
//			System.out.println("Processing Chr " + chr);
//			ArrayList<String> files = new ArrayList<String>();
//			for (int d = 0; d < datasetOutput.length; d++) {
//				String file = datasetOutput[d] + "genotypes-filtered-sorted-Chr" + chr + ".vcf.gz";
//				if (Gpio.exists(file)) {
//					files.add(file);
//				}
//			}
//
//			VCFFunctions v = new VCFFunctions();
//			if (files.size() != datasetOutput.length) {
//				System.out.println("Found only: " + files.size() + " files for chr: " + chr);
//			} else {
//				String file1 = files.get(0);
//				String file2 = files.get(1);
//				v.mergeAndIntersectVCFVariants(file1,
//						file2,
//						outputdir + "tmp1-Chr" + chr + ".vcf.gz",
//						outputdir + "tmp2-Chr" + chr + ".vcf.gz",
//						outputdir + "mergedtmp-Chr" + chr + ".vcf.gz",
//						"/",
//						outputdir + "mergelog-Chr" + chr + "-" + 0 + ".txt.gz",
//						false);
//				if (files.size() > 2) {
//					for (int i = 2; i < files.size(); i++) {
//						file2 = files.get(i);
//						v.mergeAndIntersectVCFVariants(outputdir + "mergedtmp-Chr" + chr + ".vcf.gz",
//								file2,
//								outputdir + "tmp1-Chr" + chr + ".vcf.gz",
//								outputdir + "tmp2-Chr" + chr + ".vcf.gz",
//								outputdir + "merged-Chr" + chr + ".vcf.gz",
//								"/",
//								outputdir + "mergelog-Chr" + chr + "-" + i + ".txt.gz",
//								false);
//						Gpio.moveFile(outputdir + "merged-Chr" + chr + ".vcf.gz", outputdir + "mergedtmp-Chr" + chr + ".vcf.gz");
//					}
//				}
//				Gpio.moveFile(outputdir + "mergedtmp-Chr" + chr + ".vcf.gz", outputdir + "merged-Chr" + chr + ".vcf.gz");
//			}
//
//
//		}
//	}

	public void liftOver(String plinkDataset, String outfile) throws IOException {
		PedAndMapFunctions pedAndMapFunctions = new PedAndMapFunctions();
		VCFFunctions vcfFunctions = new VCFFunctions();
		ProcessBuilder pb = null;
		GenotypeTools t = new GenotypeTools();


		String bedout = outfile + "_hg18.bed";
		pedAndMapFunctions.rewriteMapToBed(plinkDataset + ".map", bedout);

		String lifted = bedout + "-lifted.bed";
		String unlifted = bedout + "-unlifted.bed";
//
		pb = new ProcessBuilder("/Data/ImmunoChip/SequencingProject/liftOver", bedout,
				"/Data/ImmunoChip/SequencingProject/hg18ToHg19.over.chain.gz", lifted, unlifted);
		ProcessRunner.run(pb);

		String hg19map = outfile + "_hg19.map";
		pedAndMapFunctions.convertPostLiftOverMAP(plinkDataset + ".map", hg19map, lifted);

		String hg19mapupd = outfile + "_hg19_dbSNP138.map";
		String dbsnpvcf = "/Data/Projects/2014-FR-Reseq/GATKResources/dbsnp_138.b37.vcf.gz";
		pedAndMapFunctions.updateRSNames(dbsnpvcf, hg19map, hg19mapupd);

		String hg19mapupddedup = outfile + "_hg19_dbSNP138_dedup.map";
		pedAndMapFunctions.deduplicateMAP(hg19mapupd, hg19mapupddedup);
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
		ProcessRunner.run(pb);

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
		ProcessRunner.run(pb);


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
		ProcessRunner.run(pb);

//		Gpio.delete(vcf);
		vcfFunctions.splitPerChromosome(sortedvcf + ".gz", outdir + "genotypes-filtered-sorted");


	}
//
////	public String preparePlinkDataset(String referenceVCF) throws IOException {
////		outdir += "gt";
//
//	// String referenceVCF = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/MergedWithIC/merged-lowfreqvariantsremoved-sorted-phased-sorted.vcf";
////
//
////
//
//
////
////		// select variants from plink
//
//
//	// split per chromosome
//
//
////		Chromosome[] chromosomes = Chromosome.values();
////		for (Chromosome chr : chromosomes) {
////			// compare VCF to reference dataset
////			vcfFunctions.mergeAndIntersectVCFVariants(
////					referenceVCF + "-" + chr.getName() + ".vcf",
////					outdir + "-filtered-sorted-" + chr.getName() + ".vcf",
////					outdir + "-reference-sorted-matched.vcf",
////					outdir + "-filtered-sorted-" + chr.getName() + "-matched.vcf",
////					outdir + "-filtered-sorted-" + chr.getName() + "-merged.vcf",
////					"/",
////					outdir + "-filtered-sorted-matched-mergelog-" + chr.getName() + ".txt",
////					true);
////		}
//
//
////		// next up: beagle compare
////		// impute unfiltered
////		String sortedvcfphased = origmap + "-filtered-sorted-phasedAndImputed";
////		pb = new ProcessBuilder("java",
////				"-jar",
////				"/Data/Tools/beagle/beagle.r1399.jar",
////				"ref=" + referenceVCF,
////				"gt=" + sortedvcf,
////				"out=" + sortedvcfphased,
////				"nthreads=4");
////		run(pb);
////		return null;
////	}

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

//	public void makeRSquaredPlots() throws IOException, DocumentException {
//
//		String sequencedRegionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/2015-03-30-allRegionsWith50PercentOfSamplesAbove20x.bed";
//		String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/2015-03-30-allRegionsWith50PercentOfSamplesAbove20xPlus1K.bed"; // "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";
//		String gtf = "/Data/Annotation/UCSC/genes.gtf";
//
//
//		String[] refs = new String[]{"seqvaronly", "seq", "kg", "kgseq"};
//		String[] datasetIds = new String[]{"RA", "T1D"};
//		for (int dataset = 0; dataset < 2; dataset++) {
//
//			for (int chromosome = 1; chromosome < 23; chromosome++) {
//				String[] filesForPlotting = new String[refs.length];
//				String[] preimputation = new String[refs.length];
//				for (int ref = 0; ref < refs.length; ref++) {
//					String afterImputation = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/" + datasetIds[dataset] + "/" + refs[ref] + "/merged-variantinfo-Chr" + chromosome + ".txt";
//					String beforImputation = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/" + datasetIds[dataset] + "/" + refs[ref] + "-preimputation-input-Chr" + chromosome + ".txt";
//					if (Gpio.exists(afterImputation)) {
//						filesForPlotting[ref] = afterImputation;
//						preimputation[ref] = beforImputation;
//					}
//
//				}
//
//				String outdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-22-Assoc/RsquaredPlots/" + datasetIds[dataset] + "/";
//				Gpio.createDir(outdir);
//				Chromosome chrin = Chromosome.parseChr("" + chromosome);
//
//				TextFile tf2 = new TextFile(regionFile, TextFile.R);
//				String[] elems = tf2.readLineElems(TextFile.tab);
//				while (elems != null) {
//					Feature region = new Feature();
//					region.setChromosome(Chromosome.parseChr(elems[0]));
//					region.setStart(Integer.parseInt(elems[1]));
//					region.setStop(Integer.parseInt(elems[2]));
//
//					if (region.getChromosome().equals(chrin)) {
//						VariantPlot plot = new VariantPlot(outdir + region.toString() + ".pdf", 1400, 1400);
//						plot.setMargin(200);
//						plot.plotImputationRSquared(filesForPlotting, refs, preimputation, sequencedRegionFile, region, gtf);
//
//					}
//					elems = tf2.readLineElems(TextFile.tab);
//				}
//				tf2.close();
//			}
//		}
//
//
//		for (int chromosome = 1; chromosome < 23; chromosome++) {
//			String[] filesForPlotting = new String[2];
//			filesForPlotting[0] = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/ImputationPanels/ref-impquals-Chr" + chromosome + ".txt";
//			filesForPlotting[1] = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/ImputationPanels/test-impquals-Chr" + chromosome + ".txt";
//
//			String[] preimputation = new String[2];
//			preimputation[0] = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/ImputationPanels/1kg-preimputation-input-Chr" + chromosome + ".txt";
//			preimputation[1] = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/ImputationPanels/seqpanel-preimputation-input-Chr" + chromosome + ".txt";
//
//			String outdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/ImputationQuals/ImputationPanels/plots/";
//			Gpio.createDir(outdir);
//			TextFile tf2 = new TextFile(regionFile, TextFile.R);
//			String[] elems = tf2.readLineElems(TextFile.tab);
//			refs = new String[]{"kg", "seqpanel"};
//			while (elems != null) {
//				Feature region = new Feature();
//				region.setChromosome(Chromosome.parseChr(elems[0]));
//				region.setStart(Integer.parseInt(elems[1]));
//				region.setStop(Integer.parseInt(elems[2]));
//
//
//				Chromosome chrin = Chromosome.parseChr("" + chromosome);
//				if (region.getChromosome().equals(chrin)) {
//					VariantPlot plot = new VariantPlot(outdir + region.toString() + ".pdf", 1400, 1400);
//					plot.setMargin(200);
//					plot.plotImputationRSquared(filesForPlotting, refs, preimputation, sequencedRegionFile, region, gtf);
//
//				}
//				elems = tf2.readLineElems(TextFile.tab);
//			}
//			tf2.close();
//		}
//
//
//	}

//	public void makeListOfMissingICVariants() throws IOException {
//
//		String regionsFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";
//
//		String rsquareDir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/ImpQScores/";
//		String[] refs = new String[]{"1kg", "seq", "1kg-seq-merged"};
//		String immunoChipT1D = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/T1D-Onengut/hg19_gwas_ic_t1d_onengut_cc_4_18_1.tab";
//		String immunoChipRA = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/RA-Eyre/hg19_gwas_ic_ra_eyre_4_18_0.tab";
//
//		;
//		String outdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/";
//
//		BedFileReader reader = new BedFileReader();
//		ArrayList<Feature> regions = reader.readAsList(regionsFile);
//
//		for (int d = 0; d < 2; d++) {
//			String ic = immunoChipRA;
//			String ds = "RA";
//			String[] subsets = new String[]{"ES", "Merged", "NL", "SEE", "SEU", "UK", "US"};
//			String origDataDir = "/Data/ImmunoChip/RA/2015-09-23-IncX/";
//			if (d == 1) {
//				ds = "T1D";
//				ic = immunoChipT1D;
//				origDataDir = "/Data/ImmunoChip/T1D/2015-09-23-IncX/";
//				subsets = new String[]{"eur", "uk", "merged"};
//			}
//			// outfile
//			// chr\tpos\trsq1\trsq2\trsq3\tnrMissing\tICPval
//
//			TextFile out = new TextFile(outdir + "ComparisonToIChip-" + ds + ".txt", TextFile.W);
//			out.writeln("Region\tChr\tPos\t" + Strings.concat(refs, Strings.tab) + "\tNrMissingPostImputation\t" + Strings.concat(subsets, Strings.tab) + "\tNrMissingPreImputation\tIcPVal");
//
//			for (int r = 0; r < regions.size(); r++) {
//				Feature region = regions.get(r);
//
//				Chromosome chr = region.getChromosome();
//
//				if (!chr.equals(Chromosome.X)) {
//					VariantPlot variantPlot = new VariantPlot();
//					Pair<HashSet<String>, ArrayList<Pair<Integer, Double>>> icPvalsData = variantPlot.readTabFile(ic, region);
//					ArrayList<Pair<Integer, Double>> icPvals = icPvalsData.getRight();
//					String[] outputStr = new String[icPvals.size()];
//					String[] outputStrSubsets = new String[icPvals.size()];
//					int[] nrMissing = new int[icPvals.size()];
//					int[] nrMissingSubsets = new int[icPvals.size()];
//					for (int p = 0; p < icPvals.size(); p++) {
//						Pair<Integer, Double> icp = icPvals.get(p);
//						outputStr[p] = region.toString() + "\t" + chr.getName() + "\t" + icp.getLeft();
//					}
//
//					for (int s = 0; s < subsets.length; s++) {
//						String subsetfile = origDataDir + subsets[s] + "/header-Chr" + chr.getNumber() + ".txt";
//
//						TextFile t = new TextFile(subsetfile, TextFile.R);
//						String ln = t.readLine();
//						HashSet<Integer> positions = new HashSet<Integer>();
//						while (ln != null) {
//
//							if (ln.startsWith("#")) {
//
//							} else {
//								String[] elems = ln.split("\t");
//								if (elems.length > 2) {
//									Integer pos = Integer.parseInt(elems[1]);
//									positions.add(pos);
//								}
//							}
//
//							ln = t.readLine();
//						}
//						t.close();
//
//						for (int p = 0; p < icPvals.size(); p++) {
//							Pair<Integer, Double> icp = icPvals.get(p);
//							Integer pos = icp.getLeft();
//							String outstr = "T";
//							if (!positions.contains(pos)) {
//								outstr = "F";
//								nrMissingSubsets[p]++;
//							}
//							outputStrSubsets[p] += "\t" + outstr;
//						}
//
//					}
//
//					for (int refId = 0; refId < refs.length; refId++) {
//						String ref = refs[refId];
//						String rsquareFile = rsquareDir + ds + "/" + ref + "/merged-variantinfo-Chr" + chr.getNumber() + ".txt";
//						if (Gpio.exists(rsquareFile)) {
//							VCFRSquares q = new VCFRSquares();
//							ArrayList<Pair<Feature, Double>> rsquareds = q.loadRSquareds(rsquareFile, region);
//
//
//							// overlap files
//							HashMap<Integer, Double> posToRSq = new HashMap<Integer, Double>();
//							for (Pair<Feature, Double> p : rsquareds) {
//								Feature f = p.getLeft();
//								Integer pos = f.getStart();
//								Double rsq = p.getRight();
//
//								if (posToRSq.containsKey(pos)) {
//									// System.err.println("Position already used: " + ref + "\t" + f.toString());
//								} else {
//									posToRSq.put(pos, rsq);
//								}
//							}
//
//							System.out.println(posToRSq.size() + " pvals loaded from: " + rsquareFile);
//
//
//							for (int p = 0; p < icPvals.size(); p++) {
//								Pair<Integer, Double> icp = icPvals.get(p);
//								Double rsq = posToRSq.get(icp.getLeft());
//								outputStr[p] = outputStr[p] + "\t" + rsq;
//								if (rsq == null) {
//									nrMissing[p]++;
////									System.out.println("Could not find position: " + icp.getLeft());
//								}
//
//							}
//
//						}
//					}
//					for (int p = 0; p < icPvals.size(); p++) {
//						Pair<Integer, Double> icp = icPvals.get(p);
//						String outStr = outputStr[p] + "\t" + nrMissing[p] + "\t" + outputStrSubsets[p] + "\t" + nrMissingSubsets[p] + "\t" + icp.getRight();
//						out.writeln(outStr);
//					}
//				}
//
//
//			}
//
//			out.close();
//
//
//		}
//	}

	public void mergeAssociationResults(Feature region, String[] names, String[] files, String[] impQualFiles, String outputfilename) throws IOException {


		ArrayList<HashMap<Feature, Pair<String, String>>> dataPVal = new ArrayList<HashMap<Feature, Pair<String, String>>>();
		ArrayList<HashMap<Feature, Pair<String, String>>> dataQualVal = new ArrayList<HashMap<Feature, Pair<String, String>>>();

		HashSet<Feature> uniqueFeats = new HashSet<Feature>();

		System.out.println("Merging into: " + outputfilename);
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

						Double log10p = Double.NaN;
						try {
							Double pvald = Double.parseDouble(elems[3]);
							log10p = -Math.log10(pvald);
						} catch (NumberFormatException e) {

						}
						Feature feat = new Feature();
						feat.setChromosome(Chromosome.parseChr(elems[1]));
						Integer pos = Integer.parseInt(elems[2]);
						String name = elems[0];

						feat.setName(name);
						feat.setStart(pos);
						feat.setStop(pos);

						if (feat.overlaps(region)) {

							Pair<String, String> pair = new Pair<String, String>(maf, "" + log10p);
							uniqueFeats.add(feat);
							p.put(feat, pair);
						}

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

						if (feat.overlaps(region)) {
							Pair<String, String> pair = new Pair<String, String>(maf, pval);
							p.put(feat, pair);
							uniqueFeats.add(feat);
						}
					}

					elems = tf.readLineElems(Strings.whitespace);
					ln++;
				}
				tf.close();

//				System.out.println(ln + " lines and " + p.size() + " variants parsed.");
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

//				System.out.println(ln + " lines and " + p.size() + " variants parsed.");
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

		ArrayList<Feature> featArr = new ArrayList<Feature>();
		featArr.addAll(uniqueFeats);
		Collections.sort(featArr, new FeatureComparator(true));

		for (Feature feat : featArr) {
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
