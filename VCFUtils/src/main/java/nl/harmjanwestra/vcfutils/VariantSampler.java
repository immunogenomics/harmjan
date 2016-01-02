package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by hwestra on 11/24/15.
 */
public class VariantSampler {

	public void sample(String vcf1, String referenceVCF, double percentage, String out) throws IOException {


		System.out.println("About to sample " + percentage + " from " + vcf1);
		// if ref available
		HashSet<String> refVariants = null;
		// create a list of reference variants (chr+start+name)
		if (referenceVCF != null) {
			System.out.println("Loading reference VCF variants.");
			refVariants = new HashSet<String>();
			refVariants.addAll(getListOfVariants(referenceVCF));
			System.out.println(refVariants.size() + " reference VCF variants found.");
		}

		// create a list of variants
		ArrayList<String> variants = getListOfVariants(vcf1);
		System.out.println(variants.size() + " variants found in " + vcf1);

		// if reference list loaded, limit selection to that list
		if (refVariants != null) {
			variants = intersect(variants, refVariants);
			System.out.println(variants.size() + " variants intersect with reference.");
		}

		// select variants from list
		HashSet<String> sample = sample(variants, percentage);

		// write list of selected variants to file
		System.out.println("List of selected variants is written to: " + out + "-variantList.txt");
		writeList(sample, out + "-variantList.txt");

		// write new vcf file and remove selected variants
		DecimalFormat format = new DecimalFormat("#.##");
		System.out.println("Output VCF: " + out + ".vcf.gz");
		filter(vcf1, sample, out + ".vcf.gz");

	}

	private void writeList(HashSet<String> sample, String out) throws IOException {
		TextFile outf = new TextFile(out, TextFile.W);
		for (String s : sample) {
			outf.writeln(s);
		}
		outf.close();
	}


	private void filter(String vcf, HashSet<String> sample, String out) throws IOException {
		TextFile tf = new TextFile(vcf, TextFile.R);
		TextFile outf = new TextFile(out, TextFile.W);
		String ln = tf.readLine();
		while (ln != null) {
			if (ln.startsWith("#")) {
				outf.writeln(ln);
			} else {
				String[] elems = Strings.tab.split(ln);
				String variant = elems[0] + "-" + elems[1] + "-" + elems[2];
				if (!sample.contains(variant)) {
					outf.writeln(ln);
				}
			}
			ln = tf.readLine();
		}
		outf.close();
		tf.close();
	}

	public ArrayList<String> getListOfVariants(String vcf) throws IOException {

		ArrayList<String> variants = new ArrayList<String>();
		VCFGenotypeData data = new VCFGenotypeData(vcf);
		System.out.println("Sampling variants with maf > 0 and CR > 0");
		while (data.hasNext()) {
			VCFVariant variant = data.next();
			if (variant.getCallrate() > 0 && variant.getMAF() > 0) {
				variants.add(variant.toString());
			}
		}
		return variants;
	}

	public ArrayList<String> intersect(ArrayList<String> input, HashSet<String> ref) {
		ArrayList<String> output = new ArrayList<String>();
		for (String s : input) {
			if (ref.contains(s)) {
				output.add(s);
			}
		}
		return output;
	}

	public HashSet<String> sample(ArrayList<String> variants, double perc) {
		int nr = (int) Math.ceil((double) variants.size() * perc);
		HashSet<String> selected = new HashSet<String>();
		System.out.println(nr + " variants to select out of " + variants.size() + " (" + ((double) nr / variants.size()) + ")");
		while (selected.size() < nr) {
			int idx = (int) Math.floor((Math.random() * variants.size()));
			if (idx < variants.size()) {
				selected.add(variants.get(idx));
			}
		}
		return selected;
	}
}
