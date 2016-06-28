package nl.harmjanwestra.utilities.peakfiles;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.PeakFeature;
import nl.harmjanwestra.utilities.enums.Strand;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 9/13/15.
 */
public class ZinbaPeaksFile {


	public ArrayList<PeakFeature> readAllPeaks(String xlsfile, boolean filterfdr, double fdrthreshold) throws IOException {

		TextFile tf = new TextFile(xlsfile, TextFile.R);
		ArrayList<PeakFeature> features = new ArrayList<PeakFeature>();
		String ln = tf.readLine();
		while (ln != null) {
			PeakFeature f = parseLine(ln, filterfdr, fdrthreshold);
			if (f != null) {
				features.add(f);
			}
			ln = tf.readLine();
		}

		tf.close();
		return features;
	}

	private PeakFeature parseLine(String ln, boolean filterfdr, double fdrthreshold) throws IOException {
		if (ln.startsWith("PEAKID\tChrom\tStart") || ln.trim().length() == 0) {
// useless header line
			/*
			0 PEAKID
			1 Chrom
			2 Start
			3 Stop
			4 Strand
			5 Sig
			6 Maxloc
			7 Max
			8 pStart
			9 pStop
			10 Median
			1 1qValue

			 */
			return null;
		} else {    // chr	start	end	length	abs_summit	pileup	-log10(pvalue)	fold_enrichment	-log10(qvalue)	name
			String[] elems = Strings.whitespace.split(ln);
			if (elems.length < 12) {
				throw new IOException("Error parsing line: " + ln);
			}

			Chromosome chr = Chromosome.parseChr(elems[1]);
			int start = Integer.parseInt(elems[8]);
			int stop = Integer.parseInt(elems[9]);
			int summit = Integer.parseInt(elems[6]); //?


			double pileup = Double.parseDouble(elems[7]);
			double p = 1 - Double.parseDouble(elems[5]);
			double log10p = -Math.log10(p);
			if (p == 0) {
				log10p = Double.MAX_VALUE;
			}

			double foldenrich = 0;
			double qval = Double.parseDouble(elems[11]);
			if (filterfdr) {
				if (qval < fdrthreshold) {
					if (qval == 0) {
						qval = Double.MAX_VALUE;
					} else {
						qval = -Math.log10(qval);
					}
					PeakFeature f = new PeakFeature(chr, start, stop, summit, pileup, log10p, foldenrich, qval);
					f.setStrand(Strand.NA);
					return f;
				} else {
					return null;
				}
			} else {
				if (qval == 0) {
					qval = Double.MAX_VALUE;
				} else {
					qval = -Math.log10(qval);
				}

				PeakFeature f = new PeakFeature(chr, start, stop, summit, pileup, log10p, foldenrich, qval);
				f.setStrand(Strand.NA);
				return f;
			}

		}


	}


}
