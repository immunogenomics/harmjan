package nl.harmjanwestra.finemapping;

import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by Harm-Jan on 06/15/16.
 */
public class MakeTableWithICResults {

	public static void main(String[] args) {
		try {
			MakeTableWithICResults m = new MakeTableWithICResults();
			String outfile = "";
			String fmAssocFile = "";
			String icTabFile = "";
			String locusFile = "";
			String tabixprefix = "/Data/Ref/beagle_1kg/1kg.phase3.v5a.chr";
			String tabixSampleFile = "";

			/// RA
//			icTabFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ImmunoBase/hg19_gwas_ic_ra_eyre_4_19_1.tab.gz";
//			fmAssocFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/RA-assoc0.3-COSMO-merged-posterior.txt.gz";
//			locusFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci-woMHC.txt";
//			outfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/ComparisonsToIC/2016-06-19-RA-LocusComparisonWithEyre.txt";
//			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, outfile);
//
//			icTabFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ImmunoBase/hg19_gwas_ic_ra_eyre_4_19_1.tab.gz";
//			fmAssocFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/RA-assoc0.3-COSMO-merged-posterior.txt.gz";
//			locusFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/RA-assoc0.3-COSMO-merged-posterior-significantloci.txt.gz";
//			outfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/ComparisonsToIC/2016-06-19-RA-LocusComparisonWithEyre-SignificantLoci.txt";
//			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, outfile);

			icTabFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ImmunoBase/hg19_gwas_ra_okada_4_19_1.tab.gz";
			fmAssocFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/RA-assoc0.3-COSMO-merged-posterior.txt.gz";
			locusFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci-woMHC.txt";
			outfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/ComparisonsToIC/2016-06-19-RA-LocusComparisonWithOkada.txt";
			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, tabixSampleFile, outfile);

			icTabFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ImmunoBase/hg19_gwas_ra_okada_4_19_1.tab.gz";
			fmAssocFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/RA-assoc0.3-COSMO-merged-posterior.txt.gz";
			locusFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/RA-assoc0.3-COSMO-significantregions-75e7.bed";
			outfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/ComparisonsToIC/2016-06-19-RA-LocusComparisonWithOkada-SignificantLoci.txt";
			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, tabixSampleFile, outfile);
//
//			// T1D
			icTabFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ImmunoBase/hg19_gwas_ic_t1d_onengut_cc_4_19_1.tab.gz";
			fmAssocFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/T1D-assoc0.3-COSMO-merged-posterior.txt.gz";
			locusFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci-woMHC.txt";
			outfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/ComparisonsToIC/2016-06-19-T1D-LocusComparisonWithOnengutCC.txt";
			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, tabixSampleFile, outfile);

			icTabFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ImmunoBase/hg19_gwas_ic_t1d_onengut_cc_4_19_1.tab.gz";
			fmAssocFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/T1D-assoc0.3-COSMO-merged-posterior.txt.gz";
			locusFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/T1D-assoc0.3-COSMO-significantregions-75e7.bed";
			outfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/ComparisonsToIC/2016-06-19-T1D-LocusComparisonWithOnengutCC-SignificantLoci.txt";
			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, tabixSampleFile, outfile);

			// ComparisonsToMeta
			icTabFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/META-assoc0.3-COSMO-merged-posterior.txt.gz";
			fmAssocFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/RA-assoc0.3-COSMO-merged-posterior.txt.gz";
			locusFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/META-assoc0.3-COSMO-significantregions-75e7.bed";
			outfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/ComparisonsToMeta/2016-06-19-RA-LocusComparisonWithOkada.txt";
			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, tabixSampleFile, outfile);

			icTabFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/META-assoc0.3-COSMO-merged-posterior.txt.gz";
			fmAssocFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/T1D-assoc0.3-COSMO-merged-posterior.txt.gz";
			locusFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/META-assoc0.3-COSMO-significantregions-75e7.bed";
			outfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/ComparisonsToMeta/2016-06-19-T1D-LocusComparisonWithMeta-SignificantLoci.txt";
			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, tabixSampleFile, outfile);

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void make(String locusfile, String icTabFile, String fmAssocFile, String tabixprefix, String tabixSampleFile, String outfile) throws IOException {

		String[] header1 = new String[]{"",
				"",
				"",
				"",
				"",
				""
				, "GWAS", "", "", "", ""
				, "", "", "", "", ""
				, "Finemap", "", "", "", ""
				, "", "", "", "", ""
		};
		String[] header2 = new String[]{"",
				"",
				"",
				"",
				"",
				""
				, "TopVariantGWAS", "", "", "", ""
				, "TopVariantFinemap", "", "", "", ""
				, "TopVariantFinemap", "", "", "", ""
				, "TopVariantGWAS", "", "", "", ""
		};
		String[] header3 = new String[]{"Locus",
				"TopVariantGWAS",
				"TopVariantFinemap",
				"Distance",
				"LD(dprime)",
				"LD(rsquared)"
				, "Alleles", "Minor", "OR", "Beta", "Pval"
				, "Alleles", "Minor", "OR", "Beta", "Pval"
				, "Alleles", "Minor", "OR", "Beta", "Pval"
				, "Alleles", "Minor", "OR", "Beta", "Pval"
		};

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(locusfile);
		AssociationFile assocfilereader = new AssociationFile();

		TextFile out = new TextFile(outfile, TextFile.W);
		out.writeln(Strings.concat(header1, Strings.tab));
		out.writeln(Strings.concat(header2, Strings.tab));
		out.writeln(Strings.concat(header3, Strings.tab));
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

				AssociationResult bestICResultInFM = findResult(bestICResult, fmResults);
				AssociationResult bestFMResultInIC = findResult(bestFMResult, icResults);

				out.writeln(
						region.toString()
								+ "\t" + bestICResult.getSnp().toString()
								+ "\t" + bestFMResult.getSnp().toString()
								+ "\t" + (bestICResult.getSnp().getStart() - bestFMResult.getSnp().getStart())
								+ "\t" + ld.getLeft()
								+ "\t" + ld.getRight()
								+ "\t" + summaryStr(bestICResult)
								+ "\t" + summaryStr(bestFMResultInIC)
								+ "\t" + summaryStr(bestFMResult)
								+ "\t" + summaryStr(bestICResultInFM)
				);
				System.out.println();

			} else {
				System.err.println("No associations for region: " + region.toString());
			}
		}
		out.close();


	}

	private Pair<VCFVariant, VCFVariant> getVariants(String tabixrefprefix, String samplesToInclude, Feature region, Feature snp1, Feature snp2) throws IOException {

		String tabixfile = tabixrefprefix + region.getChromosome().getNumber() + ".vcf.gz";

		boolean[] samplesToIncludeArr = null;

		if (samplesToInclude != null) {
			HashSet<String> sampleIncludeHash = new HashSet<String>();
			TextFile tf = new TextFile(samplesToInclude, TextFile.R);
			sampleIncludeHash.addAll(tf.readAsArrayList());
			tf.close();

			VCFGenotypeData g = new VCFGenotypeData(tabixfile);
			ArrayList<String> samplesInVCF = g.getSamples();
			samplesToIncludeArr = new boolean[samplesInVCF.size()];
			for (int i = 0; i < samplesInVCF.size(); i++) {
				if (sampleIncludeHash.contains(samplesInVCF.get(i))) {
					samplesToIncludeArr[i] = true;
				}
			}

		}

		TabixReader reader = new TabixReader(tabixfile);
		TabixReader.Iterator window = reader.query(region.getChromosome().getNumber() + ":" + (region.getStart() - 10) + "-" + (region.getStop() + 10));
		String next = window.next();
		VCFVariant variant1 = null;
		VCFVariant variant2 = null;

		int nr = 0;
		while (next != null) {
			VCFVariant variant = new VCFVariant(next, VCFVariant.PARSE.HEADER);
			if (variant.asFeature().overlaps(snp1)) {
				if (variant.getId().equals(snp1.getName())) {

					variant1 = new VCFVariant(next, VCFVariant.PARSE.ALL, samplesToIncludeArr);
				}
			}
			if (variant.asFeature().overlaps(snp2)) {
				if (variant.getId().equals(snp2.getName())) {
					variant2 = new VCFVariant(next, VCFVariant.PARSE.ALL, samplesToIncludeArr);
				}
			}
			next = window.next();
			nr++;
		}
		reader.close();
		System.out.println(nr + " variants total for region in vcf: " + tabixfile);
		System.out.println((variant1 != null) + "\t" + (variant2 != null));
		return new Pair<>(variant1, variant2);
	}

	private String summaryStr(AssociationResult result) {
		// , "Alleles", "Minor", "OR", "Beta", "Pval"

		if (result == null) {
			return "NA\tNA\tNA\tNA\tNA";
		}

		double[] beta = result.getBeta();
		if (beta == null) {
			beta = new double[]{Math.log(result.getORs()[0])};
		}

		String output = Strings.concat(result.getSnp().getAlleles(), Strings.forwardslash)
				+ "\t" + result.getSnp().getMinorAllele()
				+ "\t" + Strings.concat(result.getORs(), Strings.semicolon)
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
