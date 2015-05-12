package nl.harmjanwestra.ngs;

import com.lowagie.text.DocumentException;
import nl.harmjanwestra.ngs.GenotypeFormats.PedAndMapFunctions;
import nl.harmjanwestra.ngs.GenotypeFormats.VCFFunctions;
import nl.harmjanwestra.ngs.graphics.VariantPlot;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import java.util.regex.Pattern;

/**
 * Created by hwestra on 4/21/15.
 */
public class PlayGround {

	public static void main(String[] args) {

		PlayGround g = new PlayGround();
		g.cd28April21();
//		g.cd28April21PlinkFiles();
		//	g.compare();

		VCFFunctions v = new VCFFunctions();
		String vcfIn = "/Data/Ref/BeagleRef/1kg.phase3.v5-filtered.merged.vcf.gz";
		String vcfOut = "";
		try {
			v.rewriteVariantsAtSamePositionAsMultiAllelic(vcfIn, vcfOut);
		} catch (IOException e) {
			e.printStackTrace();
		}

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

			plotRegion(vcfFiles, vcfNames, regionFile, sequencedRegions, output);
		} catch (Exception e) {

		}
	}

	public void compare() {
		try {
			String vcf1 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/5-1KG-SequencedRegions/filtered.vcf";
			String vcf1out = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/test/1kg.vcf";
			String vcf2 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/4-ICAndSeqVariantMerged/merged.vcf";
			String vcf2out = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/test/seq.vcf";
			String log = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/test/log.txt";
			VCFFunctions v = new VCFFunctions();
			v.compareAndCorrectVCFVariants(
					vcf1, vcf1out, vcf2, vcf2out, log, false, true);

			vcf2 = vcf2out;
			vcf2out = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/test/seq2.vcf";
			log = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/test/log2.txt";

			v.compareAndCorrectVCFVariants(
					vcf1, vcf1out, vcf2, vcf2out, log, false, true);


		} catch (Exception e) {

		}
	}

	public void cd28April21() {

		GenotypeTools t = new GenotypeTools();
		VCFFunctions v = new VCFFunctions();
		PedAndMapFunctions p = new PedAndMapFunctions();

		String outputPath = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/";

		// refilter the reference VCF.
		String startVCF = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-11-MixupsFixed/merged-ICIds-MixupsFixed.vcf";
		String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/cd28region.bed";
		String mergedhg19immunochipPed = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19/merge";
		String mergedhg19immunochipMap = mergedhg19immunochipPed + ".map";
		String mergedhg19immunochipFAM = mergedhg19immunochipPed + ".fam";
		String plink = "/Data/Tools/plink-1.07-mac-intel/plink1.9";
		String beagle = "/Data/Tools/beagle/beagle.r1399.jar";
		String merged1kg = "/Data/Ref/BeagleRef/1kg.phase3.v5-filtered.merged.vcf.gz";
		String bcftools = "/Data/Tools/bcftools/bcftools-1.2/bcftools";
		boolean skipimpute = true;

		try {
			// filter bad variants
			String filteredVCFOut = outputPath + "1-SeqVariantFilter/";
			Gpio.createDir(filteredVCFOut);
			v.filterLowFrequencyVariants(startVCF, filteredVCFOut, true, 10, 30, 0.90, 5);
			String seqvcfwithmultiallelicvariantssplit = filteredVCFOut +"filtered-splitmultiallelic.vcf";
			v.splitMultipleAllelicVariants(filteredVCFOut + "filtered.vcf", seqvcfwithmultiallelicvariantssplit);

			// compare against IC genotypes
			//t.compareVCFGenotypesToPedAndMap(outputPath + "1-SeqVariantFilter/filtered.vcf", mergedhg19immunochipPed, filteredVCFOut, true);

			// filter for CD28 region
			String filteredVariantOut = outputPath + "2-SeqVariantRegionFilter/";
			Gpio.createDir(filteredVariantOut);

			v.filterVCFForBedRegions(seqvcfwithmultiallelicvariantssplit, filteredVariantOut + "filtered.vcf", regionFile);
			v.summarizeVCF(filteredVariantOut + "filtered.vcf", filteredVariantOut + "summary.txt");

			// filter immunochip for sequenced regions and samples
			String regionFilteredMapOut = outputPath + "3-ICRegionFiltered/";
			Gpio.createDir(regionFilteredMapOut);
			String variantsToKeep = regionFilteredMapOut + "variantsOverlappingSequencedRegions.txt";
			p.filterMap(mergedhg19immunochipMap, regionFile, variantsToKeep);
			String variantsToExclude = regionFilteredMapOut + "ICvariantsOverlappingSequencedVariants.txt";
			p.filterMapForVCFVariants(mergedhg19immunochipMap, filteredVariantOut + "filtered.vcf", variantsToExclude);

			String regionfilteredpedOut = regionFilteredMapOut + "filtered";
			ProcessBuilder pb = new ProcessBuilder(plink,
					"--extract", variantsToKeep,
					"--file", mergedhg19immunochipPed,
					"--exclude", variantsToExclude,
					"--recode",
					"--out", regionfilteredpedOut,
					"--keep", "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19variantsinsequencedregions/mergefam.txt"
			);
			t.run(pb);

			// merge with immunochip
			String ICAndSeqVariantMerged = outputPath + "4-ICAndSeqVariantMerged/";
			Gpio.createDir(ICAndSeqVariantMerged);
			v.mergeWithPed(regionfilteredpedOut, filteredVariantOut + "filtered.vcf", ICAndSeqVariantMerged + "merged.vcf");
			v.summarizeVCF(ICAndSeqVariantMerged + "merged.vcf", ICAndSeqVariantMerged + "merged-summary.txt");

			// filter 1000 genomes for sequenced regions
			String kgSeqRegions = outputPath + "5-1KG-SequencedRegions/";
			Gpio.createDir(kgSeqRegions);
			// v.filterVCFForBedRegions(merged1kg, kgSeqRegions + "filtered.vcf", regionFile);
			v.summarizeVCF(kgSeqRegions + "filtered.vcf", kgSeqRegions + "filtered-summary.txt");

			// compare reference and 1kg
			String matchedPanelsOut = outputPath + "6-PanelsMatched/";
			Gpio.createDir(matchedPanelsOut);
			v.compareAndCorrectVCFVariants(
					kgSeqRegions + "filtered.vcf",
					matchedPanelsOut + "1kg.vcf",
					ICAndSeqVariantMerged + "merged.vcf",
					matchedPanelsOut + "seqpanel.vcf",
					matchedPanelsOut + "seqpanelto1kgcomparisonlog.txt",
					true,
					false);

			// check the converted panel against itself.
			v.compareAndCorrectVCFVariants(
					ICAndSeqVariantMerged + "merged.vcf",
					ICAndSeqVariantMerged + "merged-recheck.vcf",
					matchedPanelsOut + "seqpanel.vcf",
					matchedPanelsOut + "seqpanel-recheck.vcf",
					matchedPanelsOut + "seqpaneltoseqpanel-recheck-comparisonlog.txt",
					true,
					false);

			v.compareAndCorrectVCFVariants(
					kgSeqRegions + "filtered.vcf",
					matchedPanelsOut + "1kg-recheck.vcf",
					matchedPanelsOut + "seqpanel.vcf",
					matchedPanelsOut + "seqpanel-recheck1kg.vcf",
					matchedPanelsOut + "seqpaneltoseqpanel-recheck-1kg-comparisonlog.txt",
					true,
					false);

			// sort vcfs
			t.sortVCF(matchedPanelsOut + "seqpanel.vcf", matchedPanelsOut + "seqpanel-sorted.vcf", matchedPanelsOut + "sort.sh");
			t.sortVCF(matchedPanelsOut + "1kg.vcf", matchedPanelsOut + "1kg-sorted.vcf", matchedPanelsOut + "sort.sh");

			// phase sequencing data
			if (!skipimpute) {
				pb = new ProcessBuilder("java",
						"-jar", beagle,
						"ped=/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19variantsinsequencedregions/merge.fam",
						"gt=" + matchedPanelsOut + "seqpanel-sorted.vcf",
						"out=" + matchedPanelsOut + "seqpanel-phased"
				);

				t.run(pb);
				t.sortVCF(matchedPanelsOut + "seqpanel-phased.vcf.gz", matchedPanelsOut + "seqpanel-phased-sorted.vcf", matchedPanelsOut + "sort.sh");
			}

			v.compareAndCorrectVCFVariants(
					kgSeqRegions + "filtered.vcf",
					matchedPanelsOut + "1kg-recheck-tophased.vcf",
					matchedPanelsOut + "seqpanel-phased-sorted.vcf",
					matchedPanelsOut + "seqpanel-phased-sorted-1kgchecked.vcf",
					matchedPanelsOut + "seqpanel-phased-sorted-1kgchecked-log.txt",
					true,
					false);


			// impute 1kg into reference panel
			String imputedPanelsOut = outputPath + "7-PanelsImputed/";
			Gpio.createDir(imputedPanelsOut);
			if (!skipimpute) {

				pb = new ProcessBuilder("java",
						"-jar", beagle,
						"ref=" + matchedPanelsOut + "1kg-sorted.vcf",
						"ped=/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19variantsinsequencedregions/merge.fam",
						"gt=" + matchedPanelsOut + "seqpanel-phased-sorted.vcf",
						"out=" + imputedPanelsOut + "seqpanel-1kgimputed"
				);
				t.run(pb);

				v.compareAndCorrectVCFVariants(
						kgSeqRegions + "filtered.vcf",
						imputedPanelsOut + "1kg-recheck-tophased.vcf",
						imputedPanelsOut + "seqpanel-1kgimputed.vcf.gz",
						imputedPanelsOut + "seqpanel-1kgimputed-compared1kg.vcf",
						imputedPanelsOut + "seqpanel-1kgimputed-compared1kg.txt",
						true,
						false);

				pb = new ProcessBuilder("java",
						"-jar", beagle,
						"ref=" + matchedPanelsOut + "1kg-sorted.vcf",
						"ped=/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19variantsinsequencedregions/merge.fam",
						"gt=" + matchedPanelsOut + "seqpanel-sorted.vcf",
						"out=" + imputedPanelsOut + "seqpanel-unphased-1kgimputed"
				);
				t.run(pb);

				v.compareAndCorrectVCFVariants(
						kgSeqRegions + "filtered.vcf",
						imputedPanelsOut + "1kg-recheck-tophased.vcf",
						imputedPanelsOut + "seqpanel-unphased-1kgimputed.vcf.gz",
						imputedPanelsOut + "seqpanel-unphased-1kgimputed-compared1kg.vcf",
						imputedPanelsOut + "seqpanel-unphased-1kgimputed-compared1kg.txt",
						true,
						false);


				t.sortVCF(imputedPanelsOut + "seqpanel-1kgimputed.vcf.gz", imputedPanelsOut + "seqpanel-1kgimputed-sorted.vcf", imputedPanelsOut + "sort.sh");


				// impute reference panel into 1kg
				pb = new ProcessBuilder("java",
						"-jar", beagle,
						"ref=" + matchedPanelsOut + "seqpanel-phased-sorted.vcf",
						"ped=/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19variantsinsequencedregions/merge.fam",
						"gt=" + matchedPanelsOut + "1kg-sorted.vcf",
						"out=" + imputedPanelsOut + "1kg-seqpanelimputed"
				);
				t.run(pb);
				t.sortVCF(imputedPanelsOut + "1kg-seqpanelimputed.vcf.gz", imputedPanelsOut + "1kg-seqpanelimputed-sorted.vcf", imputedPanelsOut + "sort.sh");
			}

			// merge with original phased datasets..
			String imputedPanelsMergedOut = outputPath + "8-PanelsImputedMerged/";
			Gpio.createDir(imputedPanelsMergedOut);
			String unimputed1kg = matchedPanelsOut + "1kg-sorted.vcf";
			String imputed1kg = imputedPanelsOut + "1kg-seqpanelimputed-sorted.vcf";
			String unimputed1kgout = imputedPanelsMergedOut + "1kg-imputedVariantsFiltered.vcf";
			v.filterVCFVariants(unimputed1kg, imputed1kg, unimputed1kgout); // remove imputed variants from the original vcf if there is any overlap

			String files = unimputed1kgout + " " + imputed1kg;
			String mergedout = imputedPanelsMergedOut + "1kg-seqpanelimputed-merged.vcf";
			String mergedoutsorted = imputedPanelsMergedOut + "1kg-seqpanelimputed-merged-sorted.vcf";
			String bashfilename = imputedPanelsMergedOut + "sortmerge.sh";
			System.out.println("concat 1");
			System.out.println(files);
			System.out.println(mergedout);
			System.out.println(mergedoutsorted);
			t.concatVCF(files, mergedout, mergedoutsorted, bashfilename); // concatenate vcf files


			// repeat on the sequenced data..
			String unimputedseq = matchedPanelsOut + "seqpanel-phased-sorted.vcf";
			String imputedseq = imputedPanelsOut + "seqpanel-1kgimputed-sorted.vcf";
			String unimputedseqout = imputedPanelsMergedOut + "seqpanel-imputedVariantsFiltered.vcf";
			v.filterVCFVariants(unimputedseq, imputedseq, unimputedseqout);

			files = unimputedseqout + " " + imputedseq;
			mergedout = imputedPanelsMergedOut + "seqpanel-1kgimputed-merged.vcf";
			mergedoutsorted = imputedPanelsMergedOut + "seqpanel-1kgimputed-merged-sorted.vcf";
			bashfilename = imputedPanelsMergedOut + "sortmerge.sh";
			System.out.println("concat 2");
			System.out.println(files);
			System.out.println(mergedout);
			System.out.println(mergedoutsorted);
			t.concatVCF(files, mergedout, mergedoutsorted, bashfilename);

			// remove variants that don't intersect
			// TODO: Beagle falsely flips alleles?
			v.compareAndCorrectVCFVariants(
					imputedPanelsMergedOut + "1kg-seqpanelimputed-merged-sorted.vcf",
					imputedPanelsMergedOut + "1kg-seqpanelimputed-merged-sorted-matched.vcf",
					imputedPanelsMergedOut + "seqpanel-1kgimputed-merged-sorted.vcf",
					imputedPanelsMergedOut + "seqpanel-1kgimputed-merged-sorted-matched.vcf",
					imputedPanelsMergedOut + "match-log.txt", false, true);

			// also compare with imputed
			v.compareAndCorrectVCFVariants(
					imputedPanelsMergedOut + "1kg-seqpanelimputed-merged-sorted.vcf",
					imputedPanelsMergedOut + "1kg-seqpanelimputed-merged-sorted-matched.vcf",
					unimputedseq,
					imputedPanelsMergedOut + "seqpanel-unimputed-comparedto-1kgseqpanelimputed.vcf",
					imputedPanelsMergedOut + "match-log-1gseqpanelimputedvsunimputedseq.txt", false, true);

			// bgzip and merge
			System.out.println("bgzip and index");
			v.replaceHeader(imputedPanelsMergedOut + "1kg-seqpanelimputed-merged-sorted-matched.vcf", imputed1kg, imputedPanelsMergedOut + "1kg-seqpanelimputed-merged-sorted-matched-header.vcf");
			t.bgzipAndIndex(imputedPanelsMergedOut + "1kg-seqpanelimputed-merged-sorted-matched-header.vcf", imputedPanelsMergedOut + "bgzipindex.sh");
			t.bgzipAndIndex(imputedPanelsMergedOut + "seqpanel-1kgimputed-merged-sorted-matched.vcf", imputedPanelsMergedOut + "bgzipindex.sh");

			// merge
			System.out.println("merge");
			String imputedPanelsMergedMergedOut = outputPath + "9-PanelsMerged/";
			Gpio.createDir(imputedPanelsMergedMergedOut);
			files = imputedPanelsMergedOut + "1kg-seqpanelimputed-merged-sorted-matched-header.vcf.gz " + imputedPanelsMergedOut + "seqpanel-1kgimputed-merged-sorted-matched.vcf.gz";
			System.out.println(files);
			t.mergeVCF(files, imputedPanelsMergedMergedOut + "merged.vcf", imputedPanelsMergedMergedOut + "merge.sh");

			// compare to original imputed files
			v.compareAndCorrectVCFVariants(
					kgSeqRegions + "filtered.vcf",
					imputedPanelsMergedMergedOut + "1kg-matched.vcf",
					imputedPanelsMergedMergedOut + "merged.vcf",
					imputedPanelsMergedMergedOut + "merged-comparedTo1KgUnimputed.vcf",
					imputedPanelsMergedMergedOut + "merged-comparedTo1KgUnimputed.txt", false, true);

			v.compareAndCorrectVCFVariants(
					matchedPanelsOut + "seqpanel.vcf",
					imputedPanelsMergedMergedOut + "seqpanel-matched.vcf",
					imputedPanelsMergedMergedOut + "merged.vcf",
					imputedPanelsMergedMergedOut + "merged-comparedToSeqUnimputed.vcf",
					imputedPanelsMergedMergedOut + "merged-comparedToSeqUnimputed.txt", false, true);

			v.summarizeVCF(imputedPanelsMergedMergedOut + "merged.vcf", imputedPanelsMergedMergedOut + "merged-summarystats.txt");

			// sort vcfs
//			t.sortVCF(imputedPanelsOut + "seqpanel-1kgimputed.vcf.gz", imputedPanelsOut + "seqpanel-1kgimputed-sorted.vcf", imputedPanelsOut + "sort.sh");
//			t.sortVCF(imputedPanelsOut + "1kg-seqpanelimputed.vcf.gz", imputedPanelsOut + "1kg-seqpanelimputed-sorted.vcf", imputedPanelsOut + "sort.sh");

			// merge vcfs
//			String imputedPanelsMergedOut = outputPath + "8-PanelsImputedMerged/";
//			Gpio.createDir(imputedPanelsMergedOut);
//			// v.megeVCFSamples(imputedPanelsOut + "1kg-seqpanelimputed.vcf", imputedPanelsOut + "seqpanel-1kgimputed.vcf", imputedPanelsMergedOut+"merged.vcf");
//			pb = new ProcessBuilder(bcftools,
//					"merge",
//					"-m", "all",
//					"-o", imputedPanelsMergedOut + "merged.vcf",
//					"-O", "v",
//					imputedPanelsOut + "seqpanel-1kgimputed-sorted.vcf",
//					imputedPanelsOut + "1kg-seqpanelimputed-sorted.vcf"
//			);
			//t.run(pb);

			// prepare immunochip samples


		} catch (Exception e) {

			e.printStackTrace();
		}

	}

	public void cd28April21PlinkFiles() {

		try {

			String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/cd28region.bed";
			String[] datasets = new String[]{
					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_ES_QCgp",
					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_NL_QCgp",
					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_SE-E_QCgp",
					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_SE-U_QCgp",
					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_UK_QCgp",
					"/Data/ImmunoChip/RA/iChip_RACI_PhaseII_US_QCgp",
			};

			String[] datasetOutput = new String[]{
					"/Data/ImmunoChip/RA/CD28/ES/",
					"/Data/ImmunoChip/RA/CD28/NL/",
					"/Data/ImmunoChip/RA/CD28/SEE/",
					"/Data/ImmunoChip/RA/CD28/SEU/",
					"/Data/ImmunoChip/RA/CD28/UK/",
					"/Data/ImmunoChip/RA/CD28/US/",
			};

			String[] refNames = new String[]{"1kg", "seq", "1kgseq"};
			String[] refSets = new String[]{
					"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/5-1KG-SequencedRegions/filtered.vcf",
					"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/6-PanelsMatched/seqpanel-phased-sorted.vcf",
					"/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-21-CD28/9-PanelsMerged/merged.vcfsorted.vcf"};

			String mergedOut = "/Data/ImmunoChip/RA/CD28/";

			prepare(refSets, refNames, datasets, datasetOutput, regionFile, mergedOut);


			datasets = new String[]{
					"/Data/ImmunoChip/T1D/eur",
					"/Data/ImmunoChip/T1D/uk"
			};

			datasetOutput = new String[]{
					"/Data/ImmunoChip/T1D/CD28/eur/",
					"/Data/ImmunoChip/T1D/CD28/uk/"
			};

			mergedOut = "/Data/ImmunoChip/T1D/CD28/";

			prepare(refSets, refNames, datasets, datasetOutput, regionFile, mergedOut);


		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public void prepare(String[] refSets, String[] refNames, String[] datasets, String[] datasetOutput, String regionFile, String mergedOut) throws IOException {
		VCFFunctions f = new VCFFunctions();
		GenotypeTools t = new GenotypeTools();
		for (int j = 0; j < refSets.length; j++) {

			HashMap<Feature, Integer> overlap = new HashMap<Feature, Integer>();

			ArrayList<String> filesToMerge = new ArrayList<String>();
//
			for (int i = 0; i < datasets.length; i++) {
				String ds = datasets[i];

				String out = datasetOutput[i] + refNames[j] + "/";
				Gpio.createDir(out);
				String filename = preparePlinkDataset(ds, null, out, refSets[j], regionFile);
				ArrayList<Feature> variants = f.getVariantsFromVCF(filename);
				filesToMerge.add(filename);
				for (Feature feat : variants) {
					Integer ct = overlap.get(feat);
					if (ct == null) {
						ct = 0;
					}
					ct++;
					overlap.put(feat, ct);
				}
			}
//
			ArrayList<Feature> featuresToInclude = new ArrayList<Feature>();
			Set<Feature> set = overlap.keySet();
			for (Feature feat : set) {
				Integer ct = overlap.get(feat);
				if (ct == datasets.length) {
					featuresToInclude.add(feat);
				}
			}

			String[] mergeList = new String[datasets.length];
			for (int i = 0; i < datasets.length; i++) {
				String file = filesToMerge.get(i);
				String outdir = datasetOutput[i];
				TextFile tf = new TextFile(file, TextFile.R);
				TextFile tfout = new TextFile(file + "-filtered.vcf", TextFile.W);

				String ln = tf.readLine();
				while (ln != null) {

					if (ln.startsWith("#")) {
						tfout.writeln(ln);
					} else {
						String[] elems = ln.split("\t");
						Feature f2 = new Feature();
						f2.setChromosome(Chromosome.parseChr(elems[0]));
						f2.setStart(Integer.parseInt(elems[1]));
						f2.setStop(Integer.parseInt(elems[1]));
						if (featuresToInclude.contains(f2)) {
							tfout.writeln(ln);
						}
					}

					ln = tf.readLine();
				}

				tf.close();
				tfout.close();

				t.bgzipAndIndex(file + "-filtered.vcf",
						outdir + "bgzipindex.sh");

				mergeList[i] = file + "-filtered.vcf.gz";


			}
			String fileMergeStr = Strings.concat(mergeList, Pattern.compile(" "));
			t.mergeVCF(fileMergeStr, mergedOut + "merged.vcf", mergedOut + "merge.sh");

			// impute
			String beagle = "/Data/Tools/beagle/beagle.r1399.jar";
			ProcessBuilder pb = new ProcessBuilder("java",
					"-jar", beagle,
					"ref=" + refSets[j],
					//"ped=/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19variantsinsequencedregions/merge.fam",
					"gt=" + mergedOut + "merged.vcf",
					"out=" + mergedOut + "merged-" + refNames[j] + "-imputed",
					"nthreads=4"
			);
			t.run(pb);
		}
	}

	public void plotRegion(String[] vcfFiles, String[] vcfNames, String regionFile, String sequencedRegions, String output) throws IOException {


		try {
			VariantPlot plot = new VariantPlot(output, 1200, 500);

			plot.setMargin(100);
			// String[] variantFiles, String[] variantFileNames, String sequencedRegionFile, String regionFile


			TextFile tf2 = new TextFile(regionFile, TextFile.R);
			String[] elems = tf2.readLineElems(TextFile.tab);
			Feature region = new Feature();
			region.setChromosome(Chromosome.parseChr(elems[0]));
			region.setStart(Integer.parseInt(elems[1]));
			region.setStop(Integer.parseInt(elems[2]));
			tf2.close();

			plot.plot(vcfFiles, vcfNames, sequencedRegions, region);

		} catch (DocumentException e) {
			e.printStackTrace();
		}


	}

	public void old() {
		//		try {
//			GenotypeTools t = new GenotypeTools();


////			t.rewriteMapToBed("/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC.map", "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC-positions.bed");
////			t.rewriteMapToBed("/Data/Projects/2014-FR-Reseq/ImmunoChip/EurT1D/eur.map", "/Data/Projects/2014-FR-Reseq/ImmunoChip/EurT1D/eur-positions.bed");
//
//			t.convertPostLiftOverMAP("/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC.map"
//					, "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC-hg19.map",
//					"/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC-positions-lifted.bed");
//
//			t.convertPostLiftOverMAP("/Data/Projects/2014-FR-Reseq/ImmunoChip/EurT1D/eur.map"
//					, "/Data/Projects/2014-FR-Reseq/ImmunoChip/EurT1D/eur-hg19.map",
//					"/Data/Projects/2014-FR-Reseq/ImmunoChip/EurT1D/eur-positions_lifted.bed");
//
//
//			String dbsnpvcf = "/Data/Projects/2014-FR-Reseq/GATKResources/dbsnp_138.b37.vcf.gz";
//			t.updateRSNames(dbsnpvcf, "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC-hg19.map", "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC-hg19-updRS.map");
//			t.updateRSNames(dbsnpvcf, "/Data/Projects/2014-FR-Reseq/ImmunoChip/EurT1D/eur-hg19.map", "/Data/Projects/2014-FR-Reseq/ImmunoChip/EurT1D/eur-hg19-updRS.map");


//
//
//
////
////			t.compareFamFileWithSampleList("/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC.fam",
////					"/Data/Projects/2014-FR-Reseq/RASequencingSamples-ReWrite-seqIdToImmunoChip.txt");
////
////			String mapfile = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/filtered/ra.map";
////			String mapfileout = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/filtered/ra_liftover.map";
////			String liftoverbed = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/filtered/ra_liftover.bed";
//////			t.convertPostLiftOverMAP(mapfile, mapfileout, liftoverbed);
//
//
//			String rsMergeFile = "/Data/dbSNP/b142/RsMergeArch.bcp";
////			String rsRemoveFile = "";
////			String mapfileIn = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/filtered/ra_liftover.map";
////			String mapfileOut = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/filtered/ra_liftover_newRSIds.map";
////			t.rewriteRSNamesUsingDBSNP();
//
//			String mapin = "/Data/Projects/2014-FR-Reseq/ImmunoChip/NARAC_reseq_immuno_March_2015/PLINK_250315_0909/NARAC_reseq_immuno_March_2015.map";
//			String bedout = "/Data/Projects/2014-FR-Reseq/ImmunoChip/NARAC_reseq_immuno_March_2015/PLINK_250315_0909/NARAC_reseq_immuno_March_2015.bed";
//
////			t.rewriteMapToBed(mapin, bedout);
//
//			String bedin = "/Data/Projects/2014-FR-Reseq/ImmunoChip/NARAC_reseq_immuno_March_2015/LiftOver/lifted.bed";
//			String mapout = "/Data/Projects/2014-FR-Reseq/ImmunoChip/NARAC_reseq_immuno_March_2015/PLINK_250315_0909/NARAC_reseq_immuno_March_2015_lifted.map";
//
//			t.convertPostLiftOverMAP(mapin, mapout, bedin);
//
//			String dbsnpvcf = "/Data/Projects/2014-FR-Reseq/GATKResources/dbsnp_138.b37.vcf.gz";
//			String mapoutrsupdate = "/Data/Projects/2014-FR-Reseq/ImmunoChip/NARAC_reseq_immuno_March_2015/PLINK_250315_0909/NARAC_reseq_immuno_March_2015_lifted_updatedRS.map";
//			t.updateRSNames(dbsnpvcf, mapout, mapoutrsupdate);
//
//			String mapoutrsmerge = "/Data/Projects/2014-FR-Reseq/ImmunoChip/NARAC_reseq_immuno_March_2015/PLINK_250315_0909/NARAC_reseq_immuno_March_2015_lifted_updatedRS_dupsMerged.map";
//			t.rewriteRSNamesUsingDBSNP(rsMergeFile, mapoutrsupdate, mapoutrsmerge);
////			String vcf = args[0];
////			String plink = args[1];
////			String out = args[2];
////			t.summarizeVCF(vcf, out + "Summary.txt");
////			t.determineVCFSummaryStatistics(vcf, out + "SampleCallrate.txt");
////			t.compareVCFGenotypesToPedAndMap(vcf, plink, out);
//
////			String seqIdToImmunoChip = "/Data/Projects/2014-FR-Reseq/RASequencingSamples-ReWrite-seqIdToImmunoChip.txt";
////			String famfileIn = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/filtered/ra.ped";
////			String famfileOut = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/filtered/ra_rewrite.ped";
////			t.rewriteFamFileSampleNames(seqIdToImmunoChip, famfileIn, famfileOut);
//
//
//		} catch (IOException e) {
//			e.printStackTrace();
//
//		}

		// 2015-03-30 finalrun analysis
		try {
			GenotypeTools t = new GenotypeTools();
			VCFFunctions vcfFunctions = new VCFFunctions();
			PedAndMapFunctions pedAndMapFunctions = new PedAndMapFunctions();

//			String mapfile1 = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC-RsIDsFromT1DStudy.map";
//			String mapfile2 = "/Data/Projects/2014-FR-Reseq/ImmunoChip/NARAC_reseq_immuno_March_2015/PLINK_250315_0909/NARAC_reseq_immuno_March_2015-RACIChromosomeIds.map";
//			String outmapfile = "/Data/Projects/2014-FR-Reseq/ImmunoChip/NARAC_reseq_immuno_March_2015/PLINK_250315_0909/NARAC_reseq_immuno_March_2015-RsIDsFromT1DStudy.map";

////			pedAndMapFunctions.updateMapFileRsIdsUsingMapFile(mapfile1, mapfile2, outmapfile);
//
//

//			String file1 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/narac/narac.map";
//			String file2 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/narac/narac2.map";
//			pedAndMapFunctions.deduplicateMAP(file1,file2);

//			String freqFile1 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/eurt1d/plink.frq";
//			String freqFile2 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/raci/plink.frq";
//			String freqFile3 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/narac/plink.frq";
//			String freqFile4 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/eurt1dandracimerged/plink.frq";
//
////			String out = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/eurvsracifreq.txt";
////			String out2 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/eurvsnaracfreq.txt";
//			String out3 = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/mergedvsnarac.txt";
////			pedAndMapFunctions.comparePlinkAlleleFrequencies(freqFile1, freqFile2, out);
////			pedAndMapFunctions.comparePlinkAlleleFrequencies(freqFile1, freqFile3, out2);
//			pedAndMapFunctions.comparePlinkAlleleFrequencies(freqFile4, freqFile3, out3);

//			String[] pedFiles = new String[]{"/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/eurt1d/eur.ped",
//					"/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/raci/racifiltered_filteredlowcallratesnps.ped"};
//			String[] mapFiles = new String[]{"/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/eurt1d/eur.map",
//					"/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/raci/racifiltered_filteredlowcallratesnps.map"};
//			String includeTheseVariantsFilter = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/variantsToMergeT1DandRACI.txt";
//			String outpedFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/eurt1dandracimerged/merge";
//			t.mergeICPedFiles(pedFiles, mapFiles, includeTheseVariantsFilter, outpedFile);

//			String sampleList = "/Data/Projects/2014-FR-Reseq/2015-finalRun/AllSequencedImmunoChipIDs.txt";
//			String fam1 = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC.fam";
//			String out = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/allSequencedSamplesRACI.txt";
//			pedAndMapFunctions.filterFAM(sampleList, fam1, out);


//			t.rewriteMapFileChromosomeNames("/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC.map",
//					"/Data/Projects/2014-FR-Reseq/ImmunoChip/NARAC_reseq_immuno_March_2015/PLINK_250315_0909/NARAC_reseq_immuno_March_2015.map",
//					"/Data/Projects/2014-FR-Reseq/ImmunoChip/NARAC_reseq_immuno_March_2015/PLINK_250315_0909/NARAC_reseq_immuno_March_2015-RACIChromosomeIds.map");

//			String vcfStart = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-03-30-CalledGenotypes/merged.vcf";
//			String seqIdToImmunoChipID = "/Data/Projects/2014-FR-Reseq/RASequencingSamples-ReWrite-seqIdToImmunoChip.txt";
//			String vcfImmunoChipIDs = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/merged-ICIds.vcf";
//
//			// replace the header
////			t.replaceHeaderWithImmunoChipIDs(vcfStart, seqIdToImmunoChipID, vcfImmunoChipIDs);


//			t.mergeICPedFiles(immunoChipPEDFiles, immunoChipMAPFiles, vcfImmunoChipIDs, mergedImmunoChipOutput);

//			String[] pedFiles = new String[]{"/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/eurt1d/eur.ped",
//					"/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/raci/racifiltered_filteredlowcallratesnps.ped"};
//			String[] mapFiles = new String[]{"/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/eurt1d/eur.map",
//					"/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/raci/racifiltered_filteredlowcallratesnps.map"};
//			String includeTheseVariantsFilter = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/variantsToMergeT1DandRACI.txt";
//			String outpedFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/eurt1dandracimerged/merge";
//			t.mergeICPedFiles(pedFiles, mapFiles, includeTheseVariantsFilter, outpedFile);
//
//			String[] pedFiles = new String[]{"/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/eurt1dandracimerged/merge.ped",
//					"/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/narac/narac-overlapwithmerged-lowcallratevariantsremoved.ped"};
//			String[] mapFiles = new String[]{"/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/eurt1dandracimerged/merge.map",
//					"/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/narac/narac-overlapwithmerged-lowcallratevariantsremoved.map"};
//			String includeTheseVariantsFilter = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/variantsToMergeT1DandRACI.txt";
//			String exclusionList = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/variantsToExcludeWhileMergingWithNARAC.txt";
//			String newVariantList = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/variantsToMergeT1DandRACIAndNarac.txt";
//			t.removeVariantsFromList(includeTheseVariantsFilter, exclusionList, newVariantList);
//			String outpedFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmerged/merge";
//			t.mergeICPedFiles(pedFiles, mapFiles, newVariantList, outpedFile);
//


//			String mergedMapfile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmerged/merge.hg18map";
//			String mergedBedfile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmerged/merge.bed";
//			pedAndMapFunctions.rewriteMapToBed(mergedMapfile, mergedBedfile);

//			String liftedmergedBedfile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmerged/merge-lifted.bed";
//			String liftedmergedMapfile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmerged/merge.map";

//			pedAndMapFunctions.convertPostLiftOverMAP(mergedMapfile
//					, liftedmergedMapfile,
//					liftedmergedBedfile);

//			String liftedmergedMapFileDedupped = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmerged/merge-hg19-dedup.map";
//			pedAndMapFunctions.identifyDuplicatesInMap(liftedmergedMapfile, liftedmergedMapFileDedupped);
//
//			String dbsnpvcf = "/Data/Projects/2014-FR-Reseq/GATKResources/dbsnp_138.b37.vcf.gz";

//			String liftedmergedMapfile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19/merge.map";
//			String outmap = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19/merge-updRS.map";

//			pedAndMapFunctions.updateRSNames(dbsnpvcf,
//					liftedmergedMapfile,
//					outmap);
//			pedAndMapFunctions.identifyDuplicatesInMap(outmap, outmap+"dups");

			// compare VCF to genotyped individuals
//			String sequencingVCF = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/merged-ICIds.vcf";
//			VCFFunctions vcfFunctions = new VCFFunctions();
//			String sequencingVCFSummary = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/merged-ICIds-summary.txt";
//			String VCFPedSampleComparison = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/merged-ICIds-summary.txt";
			//vcfFunctions.summarizeVCF(sequencingVCF, sequencingVCFSummary);
////
//			String pedfile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19/merge";
//			boolean onlySharedSamples = true;
			//t.compareVCFGenotypesToPedAndMap(sequencingVCF, pedfile, VCFPedSampleComparison, onlySharedSamples);

//			String vcfGenotypes = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/merged-ICIds-summary.txtvcfGenotypes.txt";
//			t.performPCAOnVCFSampleCorrelationMatrix(vcfGenotypes);

//			String sampletoSampleList = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/SampleMixupsCorrected.txt";

//			VCFPedSampleComparison = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/MixupsFixed/merged-ICIds-MixupsFixed.txt";
//			sequencingVCFSummary = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/merged-ICIds-MixupsFixed-summary.txt";
//			t.replaceSampleNames(sequencingVCF, sampletoSampleList, vcfOut);
//			t.compareVCFGenotypesToPedAndMap(vcfOut, pedfile, VCFPedSampleComparison, onlySharedSamples);
//			vcfFunctions.summarizeVCF(vcfOut, sequencingVCFSummary);

//			String vcfOut = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/MixupsFixed/merged-ICIds-MixupsFixed.vcf";
//			String lowFreqVariantOut = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/MixupsFixed/stats/";
//			vcfFunctions.filterLowFrequencyVariants(vcfOut, lowFreqVariantOut, true);
//
//
//			String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/loci.txt";
//			String mapFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19/merge.map";
//			String variantsToKeep = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19/variantsOverlappingSequencedRegions.txt";
//			pedAndMapFunctions.filterMap(mapFile, regionFile, variantsToKeep);
//			String samplesToKeep = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/MixupsFixed/merged-ICIds-MixupsFixed.txtsharedSamples.txt";
//			String famFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19/merge.fam";
//			String out = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19/samplesToKeep.txt";
//			pedAndMapFunctions.filterFAM(samplesToKeep, famFile, out);
//
//			String famFileIn = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19variantsinsequencedregions/merge.fam";
//			String samplesToExclude = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/MixupsFixed/stats/incompleteTrios.txt";
//			pedAndMapFunctions.filterOutIncompleteTrios(famFileIn, samplesToExclude);

//			String vcfIn = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/MixupsFixed/stats/filtered.vcf";
//			String vcfVariantsToFilter = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/MixupsFixed/stats/violations.vcf";
//			String vcfOut = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/PoorVariantsFiltered/filtered.vcf";
//			vcfFunctions.filterVCF(vcfIn, vcfVariantsToFilter, vcfOut);
//			String outdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/PoorVariantsFiltered/filtered";
//			vcfFunctions.summarizeVCF(vcfOut, outdir);

//			String mapIn = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19variantsinsequencedregions/merge.map";
//			String vcfIn = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/PoorVariantsFiltered/filtered.vcf";
//			String listOfVariantsToExclude = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19variantsinsequencedregions/variantsToExclude.txt";
//			pedAndMapFunctions.filterMapForVCFVariants(mapIn, vcfIn, listOfVariantsToExclude);


//			String ped = "/Data/Projects/2014-FR-Reseq/2015-finalRun/ImmunoChipDataFiltered/allmergedhg19variantsinsequencedregions/merge-variantsUniqueToIC-dupsremoved";
//			String vcfIn = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/PoorVariantsFiltered/filtered.vcf";
//			String vcfOut = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/MergedWithIC/merged-lowfreqvariantsRemoved.vcf";
//
//			vcfFunctions.mergeWithPed(ped, vcfIn, vcfOut);
////
//			for (int i = 22; i < 23; i++) {
//				String reference = "/Data/Ref/BeagleRef/chr"+i+".1kg.phase3.v5.vcf.gz";
//				String referenceOut = "/Data/Ref/BeagleRef/chr"+i+".1kg.phase3.v5.vcf-filtered.gz";
//				String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/loci.txt";
//				vcfFunctions.filterBeagleReference(reference, referenceOut, regionFile);
//			}


			// filter IC map
//			String map = "/Data/ImmunoChip/US/raci_us.hg18map";
//			String refmap = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC.original.map";
//			String chrupdate = "/Data/ImmunoChip/US/raci_us-raciChrNames.map";
//			pedAndMapFunctions.rewriteMapFileChromosomeNames(refmap, map, chrupdate);
//
//			String rsnameref = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC-RsIDsFromT1DStudy.map";
//			String rsupdate = "/Data/ImmunoChip/US/raci_us-raciChrNames.map";
//			pedAndMapFunctions.updateMapFileRsIdsUsingMapFile(rsnameref, map, rsupdate);
//			String rsupdatededup = "/Data/ImmunoChip/US/raci_us-raciChrNames-dedup.map";
//			pedAndMapFunctions.deduplicateMAP(rsupdate, rsupdatededup);
//			String bedout = "/Data/ImmunoChip/US/raci_us-raciChrNames-dedup.bed";
//			pedAndMapFunctions.rewriteMapToBed(rsupdatededup, bedout);
//			String hg19map = "/Data/ImmunoChip/US/raci_us-raciChrNames-hg19.map";
//			String liftbed = "/Data/ImmunoChip/US/raci_us-raciChrNames-dedup-lifted.bed";
//			pedAndMapFunctions.convertPostLiftOverMAP(rsupdatededup, hg19map, liftbed);
//			String hg19mapdedup = "/Data/ImmunoChip/US/raci_us-raciChrNames-hg19-dedup.map";
//			pedAndMapFunctions.identifyDuplicatesInMap(hg19map, hg19mapdedup);
//
//			String variantSelect = "/Data/ImmunoChip/US/raci_us-raciChrNames-hg19-dedup-selectVariants.txt";
//			String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/loci.txt";
//			pedAndMapFunctions.filterMap(hg19map, regionFile, variantSelect);
//			String ped = "/Data/ImmunoChip/US/raci_us_filtered";
//			String vcf = "/Data/ImmunoChip/US/raci_us_filtered.vcf";
//			vcfFunctions.convertPEDToVCF(ped, vcf);

//			String origmap = "/Data/ImmunoChip/US/raci_us";
//			String removeFile = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/RAMatchingSample.txt";
//			String outdir = "/Data/ImmunoChip/US/";
//			t.preparePlinkDataset(origmap, removeFile, outdir);

//			String famin = "/Data/ImmunoChip/US/imputationresult/unimputed_unfiltered.fam";
//			String refFam = "/Data/ImmunoChip/iChip_RACI_PhaseII_US_QCgp.fam";
//			String famout = "/Data/ImmunoChip/US/imputationresult/unimputed_unfiltered-wPheno.fam";
////			pedAndMapFunctions.replaceFamPhenotypes(famin, refFam, famout);
//
//			String assoc = "/Data/ImmunoChip/US/imputationresult/plink.assoc";
//			String assocout = "/Data/ImmunoChip/US/imputationresult/plink-formanhattan.assoc";
//			t.rewritePlinkResults(assoc, assocout);
////
//			String vcf = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/MixupsFixed/merged-ICIds-MixupsFixed.vcf";
//			String imputedOnly = "/Data/ImmunoChip/US/imputationresult/plink-formanhattan-imputedonly.assoc.txt";
//			String filter = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/MixupsFixed/stats/variantsToKeep.txt";
//			t.getAssocForImputedVariants(assocout, vcf, filter, imputedOnly);


//			String vcfIn = "";
//			String vcfOut = "";
//		}
		} catch (Exception e) {

		}
	}

	public String preparePlinkDataset(String origmap, String removeFile, String outdir, String referenceVCF, String regionFile) throws IOException {
		PedAndMapFunctions pedAndMapFunctions = new PedAndMapFunctions();
		VCFFunctions vcfFunctions = new VCFFunctions();
		ProcessBuilder pb = null;
		GenotypeTools t = new GenotypeTools();
//		// String referenceVCF = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/MergedWithIC/merged-lowfreqvariantsremoved-sorted-phased-sorted.vcf";
////
//		String map = origmap + ".map";
//		String refmap = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC.original.map";
//		String chrupdate = "/Data/ImmunoChip/RA/US/raci_us-raciChrNames.map";
//		pedAndMapFunctions.rewriteMapFileChromosomeNames(refmap, map, chrupdate);
////
//		String rsnameref = "/Data/Projects/2014-FR-Reseq/ImmunoChip/RA_US/immunochip_us_preqc/Immunochip_RACI_PhaseII_US_PreQC-RsIDsFromT1DStudy.map";
//		String rsupdate = origmap + "-raciChrNames.map";
//		pedAndMapFunctions.updateMapFileRsIdsUsingMapFile(rsnameref, map, rsupdate);
//
//
//		String rsupdatededup = origmap + "-raciChrNames-dedup.map";
//		pedAndMapFunctions.deduplicateMAP(rsupdate, rsupdatededup);
//
//		String bedout = origmap + "-raciChrNames-dedup.bed";
//		pedAndMapFunctions.rewriteMapToBed(rsupdatededup, bedout);
//
//		// liftover
//
//		String lifted = origmap + "-lifted.bed";
//		String unlifted = origmap + "-unlifted.bed";
//
//		pb = new ProcessBuilder("/Data/Projects/2014-FR-Reseq/ImmunoChip/liftOver", bedout,
//				"/Data/Projects/2014-FR-Reseq/ImmunoChip/hg18ToHg19.over.chain.gz", lifted, unlifted);
//		t.run(pb);
//
//		String hg19map = origmap + "-raciChrNames-hg19.map";
//
//		pedAndMapFunctions.convertPostLiftOverMAP(rsupdatededup, hg19map, lifted);
//		String hg19mapdedup = origmap + "-raciChrNames-hg19-dedup.map";
//		pedAndMapFunctions.identifyDuplicatesInMap(hg19map, hg19mapdedup);
//
//		String hg19mapupd = origmap + "-raciChrNames-hg19-updRS.map";
//		String dbsnpvcf = "/Data/Projects/2014-FR-Reseq/GATKResources/dbsnp_138.b37.vcf.gz";
//		pedAndMapFunctions.updateRSNames(dbsnpvcf, hg19map, hg19mapupd);
//
////
//		String variantSelect = origmap + "-raciChrNames-hg19-dedup-selectVariants.txt";
////		String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/loci.txt";
//		pedAndMapFunctions.filterMap(hg19mapupd, regionFile, variantSelect);
//
//		Gpio.copyFile(origmap + ".map", origmap + ".hg18map");
//		Gpio.copyFile(hg19mapupd, origmap + ".map");
//
//		// select variants from plink
//		if (removeFile != null) {
//			pb = new ProcessBuilder("/Data/Tools/plink-1.07-mac-intel/plink1.9",
//					"--extract", variantSelect,
//					"--file", origmap,
//					"--recode", "--out", outdir + "-filtered",
//					"--remove", removeFile);
//		} else {
//			pb = new ProcessBuilder("/Data/Tools/plink-1.07-mac-intel/plink1.9",
//					"--extract", variantSelect,
//					"--file", origmap,
//					"--recode", "--out", outdir + "-filtered");
//		}
//		t.run(pb);
//
//		Gpio.copyFile(origmap + ".hg18map", origmap + ".map");

//		String ped = outdir + "-filtered";
//		String vcf = outdir + "-filtered.vcf";
//		vcfFunctions.convertPEDToVCF(ped, vcf);
//
		String sortedvcf = origmap + "-filtered-sorted.vcf";
//		String bashCommand = "cat " + vcf + " | /Data/Tools/vcftools/bin/vcf-sort > " + sortedvcf;
//		String bashfilename = origmap + "-filter-sort.sh";
//		TextFile tf = new TextFile(bashfilename, TextFile.W);
//		tf.writeln("#!/bin/bash\n" + bashCommand);
//		tf.close();
//		pb = new ProcessBuilder("bash", bashfilename);
//		t.run(pb);
//
//		// next up: beagle compare
//
		vcfFunctions.compareAndCorrectVCFVariants(referenceVCF,
				outdir + "refadj.vcf",
				sortedvcf,
				outdir + "vcf-sorted-matched.vcf",
				outdir + "matchlog.txt", false, false);
		return outdir + "vcf-sorted-matched.vcf";
		// impute filtered
//		pb = new ProcessBuilder("java",
//				"-jar",
//				"/Data/Tools/beagle/beagle.r1399.jar",
//				"ref=" + referenceVCF,
//				"gt=" + mergedsorted,
//				"out=" + tmpout + "phasedAndImputed",
//				"nthreads=4");
//		run(pb);
//
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

	}
}
