package nl.harmjanwestra.goscializer;

import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 3/28/16.
 */
public class Main {

	public static void main(String[] args) {

		try {

			String snplist = "/Users/hwestra/Downloads/SNPs20160328.txt";
			String dbsnp = "/Data/Projects/2014-FR-Reseq/GATKResources/dbsnp_138.b37.vcf.gz";
			String out = "/Users/hwestra/Downloads/SNPs20160328-snpmap.txt";

			Main m = new Main();
			m.getAnnotationForSNPList(snplist, dbsnp, out);


		} catch (IOException e) {
			e.printStackTrace();

		}

	}

	public void getAnnotationForSNPList(String snplist, String dbsnpVCF, String out) throws IOException {
		TextFile tf = new TextFile(snplist, TextFile.R);
		ArrayList<String> rsIds = new ArrayList<String>();
		tf.close();

		HashSet<String> rsIdHash = new HashSet<String>();
		rsIdHash.addAll(rsIds);

		HashMap<String, String> rsIdToAnnotation = new HashMap<String, String>();

		VCFGenotypeData data = new VCFGenotypeData(dbsnpVCF);
		while (data.hasNext()) {
			VCFVariant variant = data.next();
			if (rsIdHash.contains(variant.getId())) {
				rsIdToAnnotation.put(variant.getId(), variant.getId() + "\t" + variant.getChr().toLowerCase() + "\t" + variant.getPos());
			}
		}

		TextFile tfout = new TextFile(out, TextFile.W);
		tfout.writeln("SNP\tChrom\tBP");
		for (String s : rsIds) {
			String output = rsIdToAnnotation.get(s);
			if (output != null) {
				tfout.writeln(output);
			}
		}
		tfout.close();

	}


}
