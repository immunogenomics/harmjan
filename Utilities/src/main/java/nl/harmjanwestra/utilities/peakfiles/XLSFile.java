package nl.harmjanwestra.utilities.peakfiles;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.PeakFeature;
import nl.harmjanwestra.utilities.features.Strand;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 8/31/15.
 */
public class XLSFile {

	public ArrayList<PeakFeature> readAllPeaks(String xlsfile, boolean filterfdr, double fdrthreshold) throws IOException {

		TextFile tf = new TextFile(xlsfile, TextFile.R);
		ArrayList<PeakFeature> features = new ArrayList<PeakFeature>();
		String ln = tf.readLine();

		double log10pq = -Math.log10(fdrthreshold);
		while (ln != null) {
			PeakFeature f = parseLine(ln, filterfdr, log10pq);
			if (f != null) {
				features.add(f);
			}
			ln = tf.readLine();
		}

		tf.close();
		return features;
	}

	private PeakFeature parseLine(String ln, boolean filterfdr, double fdrthreshold) throws IOException {
		if (ln.startsWith("chr\tstart\tend") || ln.trim().length() == 0) {
// useless header line
			return null;
		} else if (ln.startsWith("#")) {
			parseHeaderLine(ln);
			return null;
		} else {
			// chr	start	end	length	abs_summit	pileup	-log10(pvalue)	fold_enrichment	-log10(qvalue)	name
			String[] elems = Strings.whitespace.split(ln);
			if (elems.length < 9) {
				throw new IOException("Error parsing line: " + ln);
			}
			Chromosome chr = Chromosome.parseChr(elems[0]);
			int start = Integer.parseInt(elems[1]);
			int stop = Integer.parseInt(elems[2]);
			int summit = Integer.parseInt(elems[4]);
			double pileup = Double.parseDouble(elems[5]);
			double log10p = Double.parseDouble(elems[6]);
			double foldenrich = Double.parseDouble(elems[7]);
			double qval = Double.parseDouble(elems[8]);
			if (filterfdr) {
				if (qval > fdrthreshold) {
					PeakFeature f = new PeakFeature(chr, start, stop, summit, pileup, log10p, foldenrich, qval);
					f.setStrand(Strand.NA);
					return f;
				} else {
					return null;
				}
			} else {
				PeakFeature f = new PeakFeature(chr, start, stop, summit, pileup, log10p, foldenrich, qval);
				f.setStrand(Strand.NA);
				return f;
			}


		}


	}

	int fragmentsize = 0;
	int nrfragmentstreatment = 0;
	int nrfragmentstreatmentpostfilter = 0;
	double redundantratetreatment = 0;

	private void parseHeaderLine(String ln) {

		//
//		if (ln.startsWith("# fragment size is determined as ")) {
//			String str = ln.replaceAll("\\# fragment size is determined as ", "");
//			str = str.replaceAll(" bps", "");
//			fragmentsize = Integer.parseInt(str);
//		} else if (ln.startsWith("# total fragments in treatment: ")) {
//			String str = ln.replaceAll("\\# total fragments in treatment: ", "");
//			nrfragmentstreatment = Integer.parseInt(str);
//		} else if (ln.startsWith("# fragments after filtering in treatment: ")) {
//			String str = ln.replaceAll("\\# fragments after filtering in treatment: ", "");
//			nrfragmentstreatmentpostfilter = Integer.parseInt(str);
//		} else if (ln.startsWith("# Redundant rate in treatment: ")) {
//			String str = ln.replaceAll("\\# fragments after filtering in treatment: ", "");
//			redundantratetreatment = Double.parseDouble(str);
//		}


	}

}
