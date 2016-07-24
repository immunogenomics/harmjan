package nl.harmjanwestra.vcfutils;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by Harm-Jan on 07/24/16.
 */
public class VCFVariantFilter {

	public void run(String vcf, String out, String listfile, boolean exclude) throws IOException {

		TextFile tf = new TextFile(listfile, TextFile.R);

		ArrayList<String> variants = tf.readAsArrayList();
		HashSet<String> allowedSet = new HashSet<String>();
		allowedSet.addAll(variants);

		tf.close();

		TextFile vcfin = new TextFile(vcf, TextFile.R);
		TextFile vcfout = new TextFile(out, TextFile.W);

		String ln = vcfin.readLine();
		while (ln != null) {
			if (ln.startsWith("#")) {
				vcfout.writeln(ln);
			} else {
				String substr = ln.substring(0, 200);
				String[] elems = substr.split("\t");
				if( (exclude && !allowedSet.contains(elems[2])) || (!exclude && allowedSet.contains(elems[2]))){
					vcfout.writeln(ln);
				}
			}
			ln = vcfin.readLine();
		}
		vcfin.close();
		vcfout.close();
	}
}
