package nl.harmjanwestra.utilities.peakfiles;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.PeakFeature;
import nl.harmjanwestra.utilities.enums.Strand;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

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

	int chrcol = -1;
	int startcol = -1;
	int endcol = -1;
	int summitcol = -1;
	int pileupcol = -1;
	int log10pcol = -1;
	int foldenrichcol = -1;
	int qvalcol = -1;

	private PeakFeature parseLine(String ln, boolean filterfdr, double fdrthreshold) throws IOException {
		if (ln.startsWith("chr\tstart\tend") || ln.trim().length() == 0) {
// useless header line

			String[] elems = ln.split("\t");
			for (int i = 0; i < elems.length; i++) {
				if (elems[i].equals("chr")) {
					chrcol = i;
				} else if (elems[i].equals("start")) {
					startcol = i;
				} else if (elems[i].equals("end")) {
					endcol = i;
				} else if (elems[i].equals("summit")) {
					summitcol = i;
				} else if (elems[i].equals("pileup")) {
					pileupcol = i;
				} else if (elems[i].equals("-log10(pvalue)")) {
					log10pcol = i;
				} else if (elems[i].equals("fold_enrichment")) {
					foldenrichcol = i;
				} else if (elems[i].equals("-log10(qvalue)")) {
					qvalcol = i;
				}
			}

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


			Chromosome chr = Chromosome.NA;
			if (chrcol != -1) {
				chr = Chromosome.parseChr(elems[chrcol]);
			}

			int start = -1;
			if (startcol != -1) {
				start = Integer.parseInt(elems[startcol]);
			}
			int stop = -1;
			if (endcol != -1) {
				stop = Integer.parseInt(elems[endcol]);
			}
			int summit = -1;
			if (summitcol != -1) {
				Integer.parseInt(elems[summitcol]);
			}

			double pileup = 0;
			if (pileupcol != -1) {
				pileup = Double.parseDouble(elems[pileupcol]);
			}

			double log10p = 0;
			if (log10pcol != -1) {
				log10p = Double.parseDouble(elems[log10pcol]);
			}
			double foldenrich = 1;
			if (foldenrichcol != -1) {
				Double.parseDouble(elems[foldenrichcol]);
			}
			double qval = 1;
			if (qvalcol != -1) {
				qval = Double.parseDouble(elems[qvalcol]);
			}
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
