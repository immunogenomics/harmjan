package nl.harmjanwestra.utilities.oxford;

import nl.harmjanwestra.utilities.vcf.VCFFunctions;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;

/**
 * Created by hwestra on 12/1/15.
 */
public class HapSample {

	public String covertHapsSampleToVCF(String haps, String vcf, boolean linux, String vcfsort) throws IOException {

		TextFile tfout = new TextFile(vcf + ".vcf.gz", TextFile.W);


		// write header
		// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT

		TextFile sampleIn = new TextFile(haps + ".sample", TextFile.R);

		sampleIn.readLine();
		sampleIn.readLine();

		String[] sampleElems = sampleIn.readLineElems(Strings.whitespace);
		String sampleStr = "";
		while (sampleElems != null) {
			sampleStr += "\t" + sampleElems[1];
			sampleElems = sampleIn.readLineElems(Strings.whitespace);
		}
		sampleIn.close();

		tfout.writeln("##fileformat=VCFv4.1");
		tfout.writeln("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" + sampleStr);

		TextFile tfIn = new TextFile(haps + ".haps", TextFile.R);

		String[] elems = tfIn.readLineElems(Strings.whitespace);
		while (elems != null) {
			// 1 rs7512482 2146966 C T 0 1
			// String ln = elems[0] + "\t" + elems[2] + "\t" + elems[1] + "\t" + elems[3] + "\t" + elems[4] + "\t.\t.\t.\tGT";
			StringBuilder ln = new StringBuilder();
			ln.append(elems[0]);
			ln.append("\t");
			ln.append(elems[2]);
			ln.append("\t");
			ln.append(elems[1]);
			ln.append("\t");
			ln.append(elems[3]);
			ln.append("\t");
			ln.append(elems[4]);
			ln.append("\t");
			ln.append(".");
			ln.append("\t");
			ln.append(".");
			ln.append("\t");
			ln.append(".");
			ln.append("\t");
			ln.append("GT");

			for (int i = 5; i < elems.length; i += 2) {
				ln.append("\t").append(elems[i]).append("|").append(elems[i + 1]);
			}

			tfout.writeln(ln.toString());
			elems = tfIn.readLineElems(Strings.whitespace);
		}

		tfIn.close();
		tfout.close();


		VCFFunctions t = new VCFFunctions();
		t.sortVCF(linux, vcfsort, vcf + ".vcf.gz", vcf + "-sorted.vcf.gz", vcf + "-sort.sh");

		return vcf + "-sorted.vcf.gz";
	}
}
