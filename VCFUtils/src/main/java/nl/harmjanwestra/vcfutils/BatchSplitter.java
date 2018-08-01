package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.legacy.genetica.containers.Triple;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

/**
 * Created by hwestra on 11/24/15.
 */
public class BatchSplitter {


	public void splitVCFOverRandomBatches(String vcfIn, String ped, String outprefix, String limitSamplefile, int nrBatches, long seed) throws IOException {

		System.out.println("splitting: " + vcfIn);
		System.out.println("ped: " + ped);
		System.out.println("out: " + outprefix);
		System.out.println("#ba: " + nrBatches);
		System.out.println("samplelimit: " + limitSamplefile);
		HashMap<String, Triple<String, String, String>> trios = new HashMap<String, Triple<String, String, String>>();

		HashSet<String> samplesToInclude = null;
		if (limitSamplefile != null) {
			samplesToInclude = new HashSet<String>();
			TextFile tf = new TextFile(limitSamplefile, TextFile.R);
			String ln = tf.readLine();
			while (ln != null) {
				samplesToInclude.add(ln.trim());
				ln = tf.readLine();
			}
			tf.close();
			System.out.println("Read: " + samplesToInclude.size() + " unique samples to select from: " + limitSamplefile);
		}

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

		ArrayList<String> selectedSamples = new ArrayList<String>();
		ArrayList<String> originalSamples = new ArrayList<>();
		HashMap<String, Integer> sampleIndex = new HashMap<String, Integer>();
		while (elems != null) {
			if (elems[0].startsWith("##")) {

			} else if (elems[0].startsWith("#CHROM")) {

				for (int i = 9; i < elems.length; i++) {
					String sample = elems[i];
					originalSamples.add(sample);
					if (samplesToInclude == null || samplesToInclude.contains(sample)) {
						selectedSamples.add(elems[i]);
					}
					sampleIndex.put(elems[i], i - 9);
				}

			} else {
				break;
			}
			elems = tf.readLineElems(TextFile.tab);

		}
		tf.close();

		System.out.println(selectedSamples.size() + " will be divided over imputation batches..");

		// reassign samples
		int remainder = selectedSamples.size() % nrBatches; // don't give the last batch less samples
		int nrSamplesPerBatch = (selectedSamples.size() - remainder) / nrBatches;
//		int nrBatches = (int) Math.ceil((double) (samples.size() - remainder) / nrSamplesPerBatch);

		int remainderDistributed = (int) Math.ceil((double) remainder / nrBatches);
		nrSamplesPerBatch += remainderDistributed; // update nr samples per batch to accomodate remainder

		System.out.println("nrBatches: " + nrBatches);

		// determine samples that are in trios first, if any
		HashMap<String, Integer> sampleToBatch = new HashMap<String, Integer>();

		HashSet<String> visitedSamples = new HashSet<String>();
		int batchctr = 0;
		TextFile batchIndexOut = new TextFile(outprefix + "-batches.txt", TextFile.W);
		System.out.println("batches will be written to: " + outprefix + "-batches.txt");

		// divide trios over batches
		for (String sample : selectedSamples) {
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
		for (String sample : selectedSamples) {
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
		int[] sampleToBatchIndex = new int[originalSamples.size()];
		int[] sampleToNewPosition = new int[originalSamples.size()];

		// initialize with -1
		for (int i = 0; i < originalSamples.size(); i++) {
			sampleToBatchIndex[i] = -1;
			sampleToNewPosition[i] = -1;
		}

		int[] batchSampleCounters = new int[nrBatches];
		for (String sample : selectedSamples) {
			Integer batch = sampleToBatch.get(sample);
			if (batch == null) {
				System.err.println("Unassigned sample: " + sample);
			} else {
				int nrSamplesInBatch = batchSampleCounters[batch];

				Integer originalSampleIndex = sampleIndex.get(sample);
				sampleToBatchIndex[originalSampleIndex] = batch;
				sampleToNewPosition[originalSampleIndex] = nrSamplesInBatch;
				batchSampleCounters[batch]++;

			}
		}

		System.out.println("Nr samples per batch: ");
		TextFile[] batchout = new TextFile[nrBatches];
		for (int batch = 0; batch < nrBatches; batch++) {
			System.out.println(batch + "\t" + batchSampleCounters[batch]);
			batchout[batch] = new TextFile(outprefix + "-batch-" + batch + ".vcf.gz", TextFile.W);
		}

		// iterate the VCF path
		tf.open();
		elems = tf.readLineElems(TextFile.tab);
		int lnctr = 0;
		while (elems != null) {
			if (elems[0].startsWith("##")) {
				// VCF headers
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
					// append variant header for this variant to each of the batch outputs
					builders[batch].append(Strings.concat(elems, Strings.tab, 0, 9));
					batchElems[batch] = new String[batchSampleCounters[batch]];
				}

				for (int i = 9; i < elems.length; i++) {
					int batchIndex = sampleToBatchIndex[i - 9];
					int batchSamplePos = sampleToNewPosition[i - 9];
					if (batchIndex != -1) {
						batchElems[batchIndex][batchSamplePos] = elems[i];
					}
				}

				// build strings for each batch
				for (int batch = 0; batch < nrBatches; batch++) {
					builders[batch].append("\t").append(Strings.concat(batchElems[batch], Strings.tab));
					batchout[batch].writeln(builders[batch].toString());
					batchout[batch].flush();
				}
			}
			elems = tf.readLineElems(TextFile.tab);
			lnctr++;
			if (lnctr % 1000 == 0) {
				System.out.print(lnctr + " lines parsed..\r");
			}
		}
		System.out.println("Done parsing. " + lnctr + " lines finally parsed.");


		tf.close();

		for (int batch = 0; batch < nrBatches; batch++) {
			batchout[batch].close();
		}

	}
}
