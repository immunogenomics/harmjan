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
import umcg.genetica.io.trityper.util.BaseAnnot;
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
			String tabixSampleFile = "/Data/Ref/1kg-europeanpopulations-meh.txt";

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
			locusFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseRALoci.bed";
			outfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/ComparisonsToIC/2016-06-19-RA-LocusComparisonWithOkada-OkadaSignificant.txt";
			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, tabixSampleFile, outfile);

//			icTabFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ImmunoBase/hg19_gwas_ra_okada_4_19_1.tab.gz";
//			fmAssocFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/RA-assoc0.3-COSMO-merged-posterior.txt.gz";
//			locusFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci-woMHC.txt";
//			outfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/ComparisonsToIC/2016-06-19-RA-LocusComparisonWithOkada.txt";
//			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, tabixSampleFile, outfile);

			icTabFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ImmunoBase/hg19_gwas_ra_okada_4_19_1.tab.gz";
			fmAssocFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/RA-assoc0.3-COSMO-merged-posterior.txt.gz";
			locusFile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/RA-assoc0.3-COSMO-significantregions-75e7.bed";
			outfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/ComparisonsToIC/2016-06-19-RA-LocusComparisonWithOkada-SignificantLoci.txt";
			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, tabixSampleFile, outfile);

			System.exit(-1);
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
				"", "",
				"",
				"", "",
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
				"", "",
				"",
				"", "",
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
				"TopVariantGWASPresentInRef", "TopVariantGWASMultiAllele",
				"TopVariantFinemap",
				"TopVariantFinemapPresentInRef", "TopVariantFinemapMultiAllele",
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
		TextFile outflip = new TextFile(outfile + "-flipped.txt", TextFile.W);
		out.writeln(Strings.concat(header1, Strings.tab));
		out.writeln(Strings.concat(header2, Strings.tab));
		out.writeln(Strings.concat(header3, Strings.tab));
		outflip.writeln(Strings.concat(header1, Strings.tab));
		outflip.writeln(Strings.concat(header2, Strings.tab));
		outflip.writeln(Strings.concat(header3, Strings.tab));
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

				boolean var1PresentInRef = true;
				boolean var1MultiAllele = false;
				boolean var2PresentInRef = true;
				boolean var2MultiAllele = false;
				if (variants.getLeft() == null) {
					var1PresentInRef = false;
				}
				if (variants.getLeft() != null && variants.getLeft().getNrAlleles() > 2) {
					var1MultiAllele = true;
				}
				if (variants.getRight() == null) {
					var2PresentInRef = false;
				}

				if (variants.getRight() != null && variants.getRight().getNrAlleles() > 2) {
					var2MultiAllele = true;
				}


				out.writeln(
						region.toString()
								+ "\t" + bestICResult.getSnp().toString()
								+ "\t" + var1PresentInRef
								+ "\t" + var1MultiAllele

								+ "\t" + bestFMResult.getSnp().toString()
								+ "\t" + var2PresentInRef
								+ "\t" + var2MultiAllele

								+ "\t" + (bestICResult.getSnp().getStart() - bestFMResult.getSnp().getStart())
								+ "\t" + ld.getLeft()
								+ "\t" + ld.getRight()
								+ "\t" + summaryStr(bestICResult, false)
								+ "\t" + summaryStr(bestFMResultInIC, false)
								+ "\t" + summaryStr(bestFMResult, false)
								+ "\t" + summaryStr(bestICResultInFM, false)
				);

				// check whether the alleles are flipped for whatever reason...


				boolean flip1 = false;
				boolean flip2 = false;
				if (variants.getLeft() != null) {
					String[] alleles = bestICResult.getSnp().getAlleles();
					String major = alleles[0];
					String minor = alleles[1];
					String snp1 = variants.getLeft().getAlleles()[0] + "/" + variants.getLeft().getAlleles()[1];
					String assessed1 = variants.getLeft().getAlleles()[0];
					String snp2 = bestICResult.getSnp().getAlleles()[0] + "/" + bestICResult.getSnp().getAlleles()[1];
					String assessed2 = minor;

					Boolean flipQ = BaseAnnot.flipalleles(snp1, assessed1, snp2, assessed2);


					if (flipQ != null) {
						flip1 = flipQ;
					}


				}
				if (variants.getRight() != null) {
					String[] alleles = bestFMResult.getSnp().getAlleles();
					String major = alleles[0];
					String minor = alleles[1];
					String snp1 = variants.getRight().getAlleles()[0] + "/" + variants.getRight().getAlleles()[1];
					String assessed1 = variants.getRight().getAlleles()[0];
					String snp2 = bestFMResult.getSnp().getAlleles()[0] + "/" + bestFMResult.getSnp().getAlleles()[1];
					String assessed2 = minor;

					Boolean flipQ = BaseAnnot.flipalleles(snp1, assessed1, snp2, assessed2);


					if (flipQ != null) {
						flip2 = flipQ;
					}

				}


				outflip.writeln(
						region.toString()
								+ "\t" + bestICResult.getSnp().toString()
								+ "\t" + bestFMResult.getSnp().toString()
								+ "\t" + (bestICResult.getSnp().getStart() - bestFMResult.getSnp().getStart())
								+ "\t" + ld.getLeft()
								+ "\t" + ld.getRight()
								+ "\t" + summaryStr(bestICResult, flip1)
								+ "\t" + summaryStr(bestFMResultInIC, flip2)
								+ "\t" + summaryStr(bestFMResult, false)
								+ "\t" + summaryStr(bestICResultInFM, false)
				);


				System.out.println();

			} else {
				System.err.println("No associations for region: " + region.toString());
			}
		}
		out.close();
		outflip.close();

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
			int n = 0;
			for (int i = 0; i < samplesInVCF.size(); i++) {
				if (sampleIncludeHash.contains(samplesInVCF.get(i))) {
					samplesToIncludeArr[i] = true;
					n++;
				}
			}

			System.out.println(n + " samples will be included.");
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
