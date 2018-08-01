package nl.harmjanwestra.finemapping.comparisons;

import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by Harm-Jan on 08/13/16.
 */
public class CompareJointAnalysisToSingleDatasetAnalyses {

	public static void main(String[] args) {
		String[] files = new String[]{};
		String[] names = new String[]{"Joint", "RA", "T1D"};
		String jointBedFile = "";
		String tabixPrefix = "";
		double ldthreshold = 0.8;
	}

	public void run(String[] assocfiles, String[] names, String jointbedfile, String tabixprefix, String out) throws IOException {

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> bedregions = reader.readAsList(jointbedfile);

		TextFile tfout = new TextFile(out, TextFile.W);
		String header = "Region";
		for (int r = 0; r < names.length; r++) {
			header += "\tSNP-" + names[r];
			header += "\tAlleles-" + names[r];
			// header += "\tAlleleFrequencies-" + names[r];
			header += "\tLog10P-" + names[r];
			header += "\tBeta-" + names[r];
			header += "\tOR-" + names[r];
			header += "\tORInterval-" + names[r];
			if (r > 0) {
				header += "\tDistance-" + names[r];
				header += "\tLD(RSquared)-" + names[r];
			}
		}
		tfout.writeln(header);

		DetermineLD ld = new DetermineLD();

		for (int f = 0; f < bedregions.size(); f++) {
			Feature region = bedregions.get(f);

			AssociationFile af = new AssociationFile();
			AssociationResult[][] results = new AssociationResult[assocfiles.length][];
			AssociationResult[] top = new AssociationResult[assocfiles.length];

			for (int i = 0; i < assocfiles.length; i++) {
				ArrayList<AssociationResult> r = af.readRegion(assocfiles[i], region);
				results[i] = r.toArray(new AssociationResult[0]);
				double maxP = 0;
				for (AssociationResult a : r) {
					if (a.getLog10Pval() > maxP) {
						top[i] = a;
						maxP = a.getLog10Pval();
					}
				}
			}

			String line = region.toString();
			for (int r = 0; r < names.length; r++) {
				line += "\t" + top[r].getSnp().toString();
				line += "\t" + Strings.concat(top[r].getSnp().getAlleles(), Strings.comma);
				// line += "\t" + top[r].getSnp().get;
				line += "\t" + top[r].getLog10Pval();
				
				// TODO: does not work for multinomial analysis
				line += "\t" + Strings.concat(top[r].getBeta()[0], Strings.semicolon);
				line += "\t" + Strings.concat(top[r].getORs()[0], Strings.semicolon);
				line += "\t" + Strings.concat(top[r].getConf(false)[0], Strings.semicolon) + "," + Strings.concat(top[r].getConf(true)[0], Strings.comma);

				if (r > 0) {
					int distance = top[0].getSnp().getStart() - top[r].getSnp().getStart();
					line += "\t" + distance;
					VCFVariant variant1 = getVariant(tabixprefix, top[0].getSnp());
					VCFVariant variant2 = getVariant(tabixprefix, top[r].getSnp());
					if (variant1 != null && variant2 != null) {
						Pair<Double, Double> ldvals = ld.getLD(variant1, variant2);
						line += "\t" + ldvals.getRight();
					} else {
						line += "\tNaN";
					}

				}
			}
			tfout.writeln(line);

		}

		tfout.close();

	}

	private VCFVariant getVariant(String tabixrefprefix, SNPFeature snp1) throws IOException {

		String tabixfile = tabixrefprefix + snp1.getChromosome().getNumber() + ".vcf.gz";
		TabixReader reader = new TabixReader(tabixfile);
		TabixReader.Iterator window = reader.query(snp1.getChromosome().getNumber() + ":" + (snp1.getStart() - 10) + "-" + (snp1.getStop() + 10));
		String next = window.next();
		VCFVariant variant1 = null;

		while (next != null) {
			VCFVariant variant = new VCFVariant(next, VCFVariant.PARSE.HEADER);
			if (variant.asFeature().overlaps(snp1)) {
				if (variant.getId().equals(snp1.getName())) {
					variant1 = new VCFVariant(next, VCFVariant.PARSE.ALL);
				}
			}
			next = window.next();
		}
		reader.close();

		if (variant1 == null) {
			System.out.println(snp1.getName() + " not found in reference");
		}
		return variant1;


	}

}
