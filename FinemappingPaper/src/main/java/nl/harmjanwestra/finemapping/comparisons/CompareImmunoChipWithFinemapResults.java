package nl.harmjanwestra.finemapping.comparisons;

import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFTabix;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by Harm-Jan on 06/15/16.
 */
public class CompareImmunoChipWithFinemapResults {

	public static void main(String[] args) {
		try {
			CompareImmunoChipWithFinemapResults m = new CompareImmunoChipWithFinemapResults();
			String outfile = "";
			String fmAssocFile = "";
			String icTabFile = "";
			String locusFile = "";
			String tabixprefix = "/Data/Ref/beagle_1kg/1kg.phase3.v5a.chrCHR.vcf.gz";
			String samplefilter = "/Data/Ref/1kg-europeanpopulations.txt.gz";


			String indir = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/NormalHWEP1e4/";
			String outdir = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/NormalHWEP1e4/2017-02-27-ComparisonsWithICAnalysis/";
			Gpio.createDir(outdir);

//			/// RA
//			icTabFile = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/ImmunoBase/hg19_gwas_ra_okada_4_19_1.tab.gz";
//			fmAssocFile = indir + "RA-assoc0.3-COSMO-merged-posterior.txt.gz";
//			locusFile = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseRALoci.bed";
//			outfile = outdir + "2016-06-19-RA-LocusComparisonWithOkada-OkadaSignificant.txt";
//			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, samplefilter, outfile);

			icTabFile = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/ImmunoBase/hg19_gwas_ra_okada_4_19_1.tab.gz";
			fmAssocFile = indir + "RA-assoc0.3-COSMO-merged-posterior.txt.gz";
			locusFile = indir + "RA-significantloci-75e7.bed";
			outfile = outdir + "2016-06-19-RA-LocusComparisonWithOkada-SignificantLoci.txt";
			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, samplefilter, outfile);

//			// T1D
//			icTabFile = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/ImmunoBase/hg19_gwas_ic_t1d_onengut_cc_4_19_1.tab.gz";
//			fmAssocFile = indir + "/T1D-assoc0.3-COSMO-merged-posterior.txt.gz";
//			locusFile = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci-woMHC.txt";
//			outfile = outdir + "2016-06-19-T1D-LocusComparisonWithOnengutCC.txt";
//			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, samplefilter, outfile);

			icTabFile = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/ImmunoBase/hg19_gwas_ic_t1d_onengut_cc_4_19_1.tab.gz";
			fmAssocFile = indir + "T1D-assoc0.3-COSMO-merged-posterior.txt.gz";
			locusFile = indir + "T1D-significantloci-75e7.bed";
			outfile = outdir + "2016-06-19-T1D-LocusComparisonWithOnengutCC-SignificantLoci.txt";
			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, samplefilter, outfile);

			// ComparisonsToMeta
			icTabFile = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/NormalHWEP1e4/META-assoc0.3-COSMO-merged-posterior.txt.gz";
			fmAssocFile = indir + "RA-assoc0.3-COSMO-merged-posterior.txt.gz";
			locusFile = indir + "META-significantloci-75e7.bed";
			outfile = outdir + "2016-06-19-RA-LocusComparisonWithMeta-SignificantLoci.txt";
			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, samplefilter, outfile);

			icTabFile = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/NormalHWEP1e4/META-assoc0.3-COSMO-merged-posterior.txt.gz";
			fmAssocFile = indir + "T1D-assoc0.3-COSMO-merged-posterior.txt.gz";
			locusFile = indir + "META-significantloci-75e7.bed";
			outfile = outdir + "2016-06-19-T1D-LocusComparisonWithMeta-SignificantLoci.txt";
			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, samplefilter, outfile);

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	enum ALLELEASSESSED {
		MINOR,
		ALTERNATIVE
	}

	;

	public void make(String locusfile,
					 String icTabFile,
					 String fmAssocFile,
					 String tabixprefix,
					 String tabixSampleFile,
//					 ALLELEASSESSED alleleassessedIC,
//					 ALLELEASSESSED alleleassessedFM,
					 String outfile) throws IOException {

		String[] header3 = new String[]{
				"Locus",
				"TopVariantFinemap",
				"AltAFControls", "AltAFCases", "ImpQual",
				"Alleles", "Minor", "OR", "Beta", "Pval",
				"TopVariantGWAS",
				"Alleles", "Minor", "OR", "Beta", "Pval",
				"Distance",
				"LD(dprime)",
				"LD(rsquared)"
		};

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(locusfile);
		AssociationFile assocfilereader = new AssociationFile();

		TextFile out = new TextFile(outfile, TextFile.W);
//		TextFile outflip = new TextFile(outfile + "-flipped.txt", TextFile.W);
		out.writeln(Strings.concat(header3, Strings.tab));
//		outflip.writeln(Strings.concat(header1, Strings.tab));
//		outflip.writeln(Strings.concat(header2, Strings.tab));
//		outflip.writeln(Strings.concat(header3, Strings.tab));
		for (int r = 0; r < regions.size(); r++) {
			System.out.println(r + "/" + regions.size());
			Feature region = regions.get(r);

			ArrayList<AssociationResult> icResults = assocfilereader.read(icTabFile, region);
			ArrayList<AssociationResult> fmResults = assocfilereader.read(fmAssocFile, region);

			AssociationResult bestICResult = getBestResult(icResults);
			AssociationResult bestFMResult = getBestResult(fmResults);

			if (bestICResult != null && bestFMResult != null) {

				// get the variants in Cosmo
				Pair<VCFVariant, VCFVariant> variants = getVariants(tabixprefix, tabixSampleFile, region, bestICResult.getSnp(), bestFMResult.getSnp());
				DetermineLD ldcalc = new DetermineLD();
				Pair<Double, Double> ld = ldcalc.getLD(variants.getLeft(), variants.getRight());

//				AssociationResult bestICResultInFM = findResult(bestICResult, fmResults);
//				AssociationResult bestFMResultInIC = findResult(bestFMResult, icResults);


				boolean flip2 = false;
				out.writeln(
						region.toString()
								+ "\t" + bestFMResult.getSnp().toString()
								+ "\t" + (1 - bestFMResult.getSnp().getAFControls())
								+ "\t" + (1 - bestFMResult.getSnp().getAFCases())
								+ "\t" + bestFMResult.getSnp().getImputationQualityScore()
								+ "\t" + summaryStr(bestFMResult, false)
								+ "\t" + bestICResult.getSnp().toString()
								+ "\t" + summaryStr(bestICResult, false)
								+ "\t" + (bestFMResult.getSnp().getStart() - bestICResult.getSnp().getStart())
								+ "\t" + ld.getLeft()
								+ "\t" + ld.getRight()
				);

				System.out.println();

			} else {
				System.err.println("No associations for region: " + region.toString());
			}
		}
		out.close();
//		outflip.close();

	}

	private Pair<VCFVariant, VCFVariant> getVariants(String tabixrefprefix, String samplesToInclude, Feature region, Feature snp1, Feature snp2) throws IOException {


		String tabixfile = tabixrefprefix.replaceAll("CHR", "" + region.getChromosome().getNumber());

		boolean[] samplesToIncludeArr = null;

		VCFTabix tabix = new VCFTabix(tabixfile);
		if (samplesToInclude != null) {
			samplesToIncludeArr = tabix.getSampleFilter(samplesToInclude);
		}

		TabixReader.Iterator window = tabix.query(region);
		String next = window.next();
		VCFVariant variant1 = null;
		VCFVariant variant2 = null;
		int nr = 0;
		while (next != null) {
			VCFVariant variant = new VCFVariant(next, VCFVariant.PARSE.HEADER);
			if (variant.asFeature().overlaps(snp1)) {
//				if (variant.getId().equals(snp1.getName())) {

				variant1 = new VCFVariant(next, VCFVariant.PARSE.ALL, samplesToIncludeArr);
//				}
			}
			if (variant.asFeature().overlaps(snp2)) {
//				if (variant.getId().equals(snp2.getName())) {
				variant2 = new VCFVariant(next, VCFVariant.PARSE.ALL, samplesToIncludeArr);
//				}
			}
			next = window.next();
			nr++;
		}

		tabix.close();
		System.out.println(nr + " variants total for region in vcf: " + tabixfile);
		System.out.println((variant1 != null) + "\t" + (variant2 != null));
		return new Pair<>(variant1, variant2);
	}

	private String summaryStr(AssociationResult result, boolean flip) {
		// , "Alleles", "Minor", "OR", "Beta", "Pval"

		if (result == null) {
			return "NA\tNA\tNA\tNA\tNA";
		}

		double[] beta = result.getBeta();
		double[] ors = result.getORs();
		String[] alleles = result.getSnp().getAlleles();
		if (beta == null) {
			beta = new double[]{Math.log(result.getORs()[0])};
		}

		if (flip) {

			String[] tmpAlleles = new String[alleles.length];
			for (int a = 0; a < alleles.length; a++) {
				tmpAlleles[a] = alleles[alleles.length - 1 - a];
			}
			alleles = tmpAlleles;
			for (int i = 0; i < beta.length; i++) {
				beta[i] = 1 / beta[i];
				ors[i] = 1 / ors[i];

			}
		}

		String output = Strings.concat(alleles, Strings.forwardslash)
				+ "\t" + result.getSnp().getMinorAllele()
				+ "\t" + Strings.concat(ors, Strings.semicolon)
				+ "\t" + Strings.concat(beta, Strings.semicolon)
				+ "\t" + result.getLog10Pval();

		return output;
	}

	private AssociationResult findResult(AssociationResult bestICResult, ArrayList<AssociationResult> results) {
		HashMap<Integer, AssociationResult> r = new HashMap<>();
		for (AssociationResult a : results) {
			int pos = a.getSnp().getStart();
			if (r.containsKey(pos)) {

			} else {
				r.put(pos, a);
			}
		}
		AssociationResult output = r.get(bestICResult.getSnp().getStart());
		return output;
	}

	private AssociationResult getBestResult(ArrayList<AssociationResult> results) {
		Double maxp = 0d;
		AssociationResult output = null;
		for (AssociationResult r : results) {
			if (r.getLog10Pval() > maxp) {
				maxp = r.getLog10Pval();
				output = r;
			}
		}
		return output;
	}

}
