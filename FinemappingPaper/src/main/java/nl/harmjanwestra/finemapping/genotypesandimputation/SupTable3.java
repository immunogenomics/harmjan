package nl.harmjanwestra.finemapping.genotypesandimputation;

import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import nl.harmjanwestra.utilities.vcf.VCFVariant;

import java.io.IOException;
import java.util.ArrayList;

public class SupTable3 {
	
	
	public static void main(String[] args) {
		String disk = "c:";
		String in = disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\panels\\hapcaller-maf0005-cr0950-rd10-gq30.vcf.gz-samplenamefix-mixupfix-nonmatchingremoved.vcf.gz";
		String out = disk + "\\Sync\\Dropbox\\FineMap\\2018-01-Rebuttal\\StagingArea\\tables\\Supplementary Table 3 - SequencedVariants.txt";
		String bedregions = disk + "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		
		SupTable3 z = new SupTable3();
		try {
			z.run(in,bedregions,out);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String in, String regionsfile, String out) throws IOException {
		BedFileReader b = new BedFileReader();
		ArrayList<Feature> regions = b.readAsList(regionsfile);
		
		TextFile tf = new TextFile(in, TextFile.R);
		TextFile tfout = new TextFile(out, TextFile.W);
		tfout.writeln("Chromosome" +
				"\tPosition" +
				"\tRS Id" +
				"\tReference Allele" +
				"\tAlternate Allele" +
				"\tAllele Frequencies" +
				"\tCallRate" +
				"\tHardy Weinberg P-value");
		String ln = tf.readLine();
		while (ln != null) {
			if (!ln.startsWith("#")) {
				VCFVariant v = new VCFVariant(ln, VCFVariant.PARSE.ALL);
				if (v.asFeature().overlaps(regions) && v.getMAF() > 0.01) {
					tfout.writeln(
							v.getChr()
									+ "\t" + v.getPos()
									+ "\t" + v.getId()
									+ "\t" + v.getAlleles()[0]
									+ "\t" + Strings.concat(v.getAlleles(), Strings.comma, 1, v.getAlleles().length)
									+ "\t" + Strings.concat(v.getAlleleFrequencies(), Strings.comma)
									+ "\t" + v.getCallrate()
									+ "\t" + v.getHwep()
					);
				}
			}
			ln = tf.readLine();
		}
		tfout.close();
		tf.close();
	}
}
