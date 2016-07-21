package nl.harmjanwestra;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashMap;

/**
 * Created by hwestra on 2/9/16.
 */
public class HomIndivFile {
	double[] avgkb;
	double[] kb;
	int[] nseg;
	int[] phe;
	String[] samples;

	public double[] getAvgkb() {
		return avgkb;
	}

	public double[] getKb() {
		return kb;
	}

	public int[] getNseg() {
		return nseg;
	}

	public int[] getPhe() {
		return phe;
	}

	public String[] getSamples() {
		return samples;
	}

	public HomIndivFile(String file, boolean removezeros, HashMap<String, String> sampleLink) throws IOException {
		System.out.println("Loading " + file);
		TextFile tf = new TextFile(file, TextFile.R);

		int nrLines = tf.countLines() - 1;

		System.out.println(nrLines + " samples in path..");

		samples = new String[nrLines];
		avgkb = new double[nrLines];
		kb = new double[nrLines];
		nseg = new int[nrLines];
		phe = new int[nrLines];
		tf.close();
		tf.open();
		tf.readLine();
		String ln = tf.readLine();
		int ctr = 0;
		while (ln != null) {
			while (ln.contains("  ")) {
				ln = ln.replaceAll("  ", " ");
			}

			if (ln.startsWith(" ")) {
				ln = ln.substring(1);
			}
			if (ln.endsWith(" ")) {
				ln = ln.substring(0, ln.length() - 1);
			}
			String[] elems = ln.split(" ");

			if (elems.length == 6) {
				String sample = elems[1];
				String strPhe = elems[2];
				String strNseg = elems[3];
				String strKb = elems[4];
				String strAvgKb = elems[5];

				samples[ctr] = sample;
				avgkb[ctr] = Double.parseDouble(strAvgKb);
				kb[ctr] = Double.parseDouble(strKb);
				nseg[ctr] = Integer.parseInt(strNseg);
				phe[ctr] = Integer.parseInt(strPhe);

				ctr++;

			}
			ln = tf.readLine();

		}

		System.out.println(ctr + " samples loaded");

		tf.close();

		if (sampleLink != null) {
			filter(sampleLink);
			System.out.println(samples.length + " samples remain after sample link");
		}

		if (removezeros) {
			filterzeros();
			System.out.println(samples.length + " samples remain after removing zeros");
		}

	}

	private void filterzeros() {

		int nrShared = 0;
		for (int i = 0; i < samples.length; i++) {
			if (avgkb[i] > 0) {
				nrShared++;
			}
		}

		double[] tmpavgkb = new double[nrShared];
		double[] tmpkb = new double[nrShared];
		int[] tmpnseg = new int[nrShared];
		int[] tmpphe = new int[nrShared];
		String[] tmpsamples = new String[nrShared];

		int idx = 0;

		for (int i = 0; i < samples.length; i++) {
			String sample = samples[i];
			if (avgkb[i] > 0) {
				tmpavgkb[idx] = avgkb[i];
				tmpkb[idx] = kb[idx];
				tmpsamples[idx] = sample;
				tmpnseg[idx] = nseg[i];
				tmpphe[idx] = phe[i];
				idx++;
			}
		}


		avgkb = tmpavgkb;
		kb = tmpkb;
		nseg = tmpnseg;
		samples = tmpsamples;
		phe = tmpphe;

	}

	private void filter(HashMap<String, String> sampleLink) {
		int nrShared = 0;
		for (int i = 0; i < samples.length; i++) {
			if (sampleLink.containsKey(samples[i])) {
				nrShared++;
			}
		}
		double[] tmpavgkb = new double[nrShared];
		double[] tmpkb = new double[nrShared];
		int[] tmpnseg = new int[nrShared];
		int[] tmpphe = new int[nrShared];
		String[] tmpsamples = new String[nrShared];

		int idx = 0;

		for (int i = 0; i < samples.length; i++) {
			String sample = samples[i];
			if (sampleLink.containsKey(sample)) {
				tmpavgkb[idx] = avgkb[i];
				tmpkb[idx] = kb[idx];
				tmpsamples[idx] = sampleLink.get(sample);
				tmpnseg[idx] = nseg[i];
				tmpphe[idx] = phe[i];
				idx++;
			}
		}


		avgkb = tmpavgkb;
		kb = tmpkb;
		nseg = tmpnseg;
		samples = tmpsamples;
		phe = tmpphe;
	}

	public void reOrder(HashMap<String, Integer> newIndex) {

		double[] tmpavgkb = new double[newIndex.size()];
		double[] tmpkb = new double[newIndex.size()];
		int[] tmpnseg = new int[newIndex.size()];
		int[] tmpphe = new int[newIndex.size()];
		String[] tmpsamples = new String[newIndex.size()];

		for (int i = 0; i < samples.length; i++) {
			String sample = samples[i];
			Integer idx = newIndex.get(sample);
			if (idx != null) {
				tmpavgkb[idx] = avgkb[i];
				tmpkb[idx] = kb[idx];
				tmpsamples[idx] = sample;
				tmpnseg[idx] = nseg[i];
				tmpphe[idx] = phe[i];
			}
		}


		avgkb = tmpavgkb;
		kb = tmpkb;
		nseg = tmpnseg;
		samples = tmpsamples;
		phe = tmpphe;

	}
}
