package nl.harmjanwestra.finemapping;

import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by Harm-Jan on 06/15/16.
 */
public class MakeTableWithICResults {

	public static void main(String[] args) {
		try {
			MakeTableWithICResults m = new MakeTableWithICResults();
			String locusFile = "D:\\tmp\\2016-06-15\\AllICLoci.bed";
			String icTabFile = "D:\\Data\\ImmunoBase\\hg19_gwas_ra_okada_4_19_1.tab";
			String fmAssocFile = "D:\\tmp\\2016-06-02\\RA-assoc0.3-COSMO-merged.txt";
			String outfile = "D:\\tmp\\2016-06-15\\testout.txt";
			String tabixprefix = "d:\\Data\\Ref\\1kg\\cosmo.1kg.phase3.v5.chr";
			m.make(locusFile, icTabFile, fmAssocFile, tabixprefix, outfile);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void make(String locusfile, String icTabFile, String fmAssocFile, String tabixprefix, String outfile) throws IOException {

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

			// get the variants in Cosmo
			Pair<VCFVariant, VCFVariant> variants = getVariants(tabixprefix, region, bestICResult.getSnp(), bestFMResult.getSnp());
			DetermineLD ldcalc = new DetermineLD();
			Pair<Double,Double> ld = ldcalc.getLD(variants.getLeft(), variants.getRight());

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
		}
		out.close();


	}

	private Pair<VCFVariant, VCFVariant> getVariants(String tabixrefprefix, Feature region, Feature snp1, Feature snp2) throws IOException {

		String tabixfile = tabixrefprefix + region.getChromosome().getNumber() + ".vcf.gz";
		TabixReader reader = new TabixReader(tabixfile);
		TabixReader.Iterator window = reader.query(region.getChromosome().getNumber() + ":" + (region.getStart() - 1000) + "-" + (region.getStart() + 1000));
		String next = window.next();
		VCFVariant variant1 = null;
		VCFVariant variant2 = null;

		while (next != null) {
			VCFVariant variant = new VCFVariant(next, VCFVariant.PARSE.HEADER);
			if (variant.asFeature().overlaps(snp1)) {
				if (variant.getId().equals(snp1.getName())) {
					variant1 = new VCFVariant(next, VCFVariant.PARSE.ALL);
				}
			}
			if (variant.asFeature().overlaps(snp2)) {
				if (variant.getId().equals(snp2.getName())) {
					variant2 = new VCFVariant(next, VCFVariant.PARSE.ALL);
				}
			}
			next = window.next();
		}
		reader.close();
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
