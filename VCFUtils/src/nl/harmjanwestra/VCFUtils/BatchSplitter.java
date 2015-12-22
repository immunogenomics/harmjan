package nl.harmjanwestra.VCFUtils;

import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

/**
 * Created by hwestra on 11/24/15.
 */
public class BatchSplitter {


	public void splitVCFOverRandomBatches(String vcfIn, String ped, String outprefix, int nrSamplesPerBatch, long seed) throws IOException {

		System.out.println("splitting: " + vcfIn);
		System.out.println("ped: " + ped);
		System.out.println("out: " + outprefix);
		System.out.println("spb: " + nrSamplesPerBatch);
		HashMap<String, Triple<String, String, String>> trios = new HashMap<String, Triple<String, String, String>>();

		if (ped != null) {
			TextFile tf = new TextFile(ped, TextFile.R);
			String[] elems = tf.readLineElems(Strings.whitespace);
			while (elems != null) {

				if (elems.length > 3) {
					String sample = elems[1];
					String f1 = elems[2];
					String f2 = elems[3];

					boolean father = false;
					if (!f1.equals("0")) {
						father = true;
					}
					boolean mother = false;
					if (!f2.equals("0")) {
						mother = true;
					}

					if (father || mother) {
						Triple<String, String, String> trio = new Triple<String, String, String>(sample, f1, f2);
						trios.put(sample, trio);
					}
				}

				elems = tf.readLineElems(Strings.whitespace);
			}

			tf.close();

			System.out.println(trios.size() + " possible trios from ped: " + ped);
		}

		TextFile tf = new TextFile(vcfIn, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);

		ArrayList<String> samples = new ArrayList<String>();
		HashMap<String, Integer> sampleIndex = new HashMap<String, Integer>();
		while (elems != null) {
			if (elems[0].startsWith("##")) {

			} else if (elems[0].startsWith("#CHROM")) {

				for (int i = 9; i < elems.length; i++) {
					samples.add(elems[i]);
					sampleIndex.put(elems[i], i - 9);
				}

			} else {
				break;
			}
			elems = tf.readLineElems(TextFile.tab);

		}
		tf.close();

		// reassign samples
		int remainder = samples.size() % nrSamplesPerBatch; // don't give the last batch less samples
		int nrBatches = (int) Math.ceil((double) (samples.size() - remainder) / nrSamplesPerBatch);

		int remainderDistributed = (int) Math.ceil((double) remainder / nrBatches);
		nrSamplesPerBatch += remainderDistributed; // update nr samples per batch to accomodate remainder

		System.out.println("nrBatches: " + nrBatches);

		// determine samples that are in trios first, if any
		HashMap<String, Integer> sampleToBatch = new HashMap<String, Integer>();

		HashSet<String> visitedSamples = new HashSet<String>();
		int batchctr = 0;
		TextFile batchIndexOut = new TextFile(outprefix + "-batches.txt", TextFile.W);
		System.out.println("batches will be written to: " + outprefix + "-batches.txt");
		for (String sample : samples) {
			if (!visitedSamples.contains(sample) && sampleIndex.containsKey(sample)) {
				Triple<String, String, String> trio = trios.get(sample);

				if (trio != null) {

					String father = trio.getMiddle();
					String mother = trio.getRight();

					boolean fatherpresent = false;
					boolean motherpresent = false;
					if (!father.equals("0")) {

						if (sampleIndex.containsKey(father)) {

							fatherpresent = true;
						}
					}
					if (!mother.equals("0")) {
						if (sampleIndex.containsKey(mother)) {

							motherpresent = true;
						}

					}

					if (fatherpresent || motherpresent) {
						// assign to batch

						int batch = batchctr % nrBatches;
						if (motherpresent) {
							sampleToBatch.put(mother, batch);
							visitedSamples.add(mother);
							batchIndexOut.writeln(mother + "\t" + batch);
						}
						if (fatherpresent) {
							sampleToBatch.put(father, batch);
							visitedSamples.add(father);
							batchIndexOut.writeln(father + "\t" + batch);
						}
						sampleToBatch.put(sample, batch);
						batchIndexOut.writeln(sample + "\t" + batch);
						visitedSamples.add(sample);
						batchctr++;
					}
				}
			}
		}

		// get a list of remaining samples
		ArrayList<String> remainingSamples = new ArrayList<String>();
		for (String sample : samples) {
			if (!visitedSamples.contains(sample)) {
				remainingSamples.add(sample);
			}
		}

		batchctr = 0;
		Random r = new Random(seed);
		// randomly assign unrelated samples to a batch number...
		while (!remainingSamples.isEmpty()) {
			String sample = remainingSamples.remove((int) (r.nextDouble() * remainingSamples.size()));
			int batch = batchctr % nrBatches;
			sampleToBatch.put(sample, batch);
			batchIndexOut.writeln(sample + "\t" + batch);
			batchctr++;
		}
		batchIndexOut.close();

		// determine final number of samples per batch and also index the samples
		int[] sampleToBatchIndex = new int[samples.size()];
		int[] sampleToNewPosition = new int[samples.size()];

		int[] batchSampleCounters = new int[nrBatches];
		for (String sample : samples) {
			Integer batch = sampleToBatch.get(sample);
			if (batch == null) {
				System.err.println("Unassigned sample: " + sample);
			} else {
				int nrSamplesInBatch = batchSampleCounters[batch];

				int previousIndex = sampleIndex.get(sample);
				sampleToBatchIndex[previousIndex] = batch;
				sampleToNewPosition[previousIndex] = nrSamplesInBatch;

				batchSampleCounters[batch]++;
			}
		}

		System.out.println("Nr samples per batch: ");
		TextFile[] batchout = new TextFile[nrBatches];
		for (int batch = 0; batch < nrBatches; batch++) {
			System.out.println(batch + "\t" + batchSampleCounters[batch]);
			batchout[batch] = new TextFile(outprefix + "-batch-" + batch + ".vcf.gz", TextFile.W);
		}


		// iterate the VCF file
		tf.open();
		elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems[0].startsWith("##")) {
				for (int batch = 0; batch < nrBatches; batch++) {
					batchout[batch].writeln(Strings.concat(elems, Strings.tab));
				}
			} else {

				StringBuilder[] builders = new StringBuilder[nrBatches];
				String[][] batchElems = new String[nrBatches][0];
				for (int batch = 0; batch < nrBatches; batch++) {
					if (builders[batch] == null) {
						builders[batch] = new StringBuilder(100000);
					}
					builders[batch].append(Strings.concat(elems, Strings.tab, 0, 9));
					batchElems[batch] = new String[batchSampleCounters[batch]];
				}

				for (int i = 9; i < elems.length; i++) {
					int batchIndex = sampleToBatchIndex[i - 9];
					int batchSamplePos = sampleToNewPosition[i - 9];
					batchElems[batchIndex][batchSamplePos] = elems[i];
				}

				for (int batch = 0; batch < nrBatches; batch++) {
					builders[batch].append("\t").append(Strings.concat(batchElems[batch], Strings.tab));
					batchout[batch].writeln(builders[batch].toString());
				}
			}
			elems = tf.readLineElems(TextFile.tab);
		}


		tf.close();

		for (int batch = 0; batch < nrBatches; batch++) {

			batchout[batch].close();
		}

	}
}
