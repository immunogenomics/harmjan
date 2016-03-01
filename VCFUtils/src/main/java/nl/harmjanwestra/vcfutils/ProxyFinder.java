package nl.harmjanwestra.vcfutils;

import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Correlation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.Callable;

/**
 * Created by hwestra on 2/29/16.
 */
public class ProxyFinder {

	public void find(String tabixrefprefix, int windowsize, double threshold, String snpfile, String output, int nrthreads) throws IOException {


		ArrayList<Feature> snps = new ArrayList<Feature>();
		TextFile tf = new TextFile(snpfile, TextFile.R);
		String ln = tf.readLine();
		while (ln != null) {

			Feature f = Feature.parseFeature(ln);

			snps.add(f);

			ln = tf.readLine();
		}
		tf.close();

		for (Feature f : snps) {


		}

	}

	public class ProxyFinderTask implements Callable<ArrayList<String>> {

		String tabixrefprefix;
		int windowsize;
		Feature f;
		double threshold;

		public ProxyFinderTask() {

		}

		@Override
		public ArrayList<String> call() throws Exception {

			ArrayList<String> output = new ArrayList<>();
			output.add(f.getName() + "\t" + f.getName() + "\t" + f.getStart() + "\t" + 0 + "\t" + 1);
			TabixReader reader = new TabixReader(tabixrefprefix + f.getChromosome().getNumber());


			TabixReader.Iterator inputSNPiter = reader.query(f.getChromosome().toString(), f.getStart() - 1, f.getStart() + 1);

			String snpStr = inputSNPiter.next();
			VCFVariant testSNPObj = null;
			System.out.println("Query: " + f.getName());
			while (snpStr != null) {
				VCFVariant variant = new VCFVariant(snpStr);
				if (variant.asFeature().overlaps(f)) {
					if (variant.getId().equals(f.getName())) {
						testSNPObj = variant;
					}
				}
				snpStr = inputSNPiter.next();
			}

			if (testSNPObj == null) {
				System.out.println("Query: " + f.getName() + " not in reference");
				return output;
			}

			VariantCorrelationMatrix correlator2 = new VariantCorrelationMatrix();
			double[] genotypes1 = correlator2.convertToDouble(testSNPObj);

			TabixReader.Iterator window = reader.query(f.getChromosome().toString(), f.getStart() - windowsize, f.getStart() + windowsize);
			String next = window.next();
			while (next != null) {
				// correlate
				VCFVariant variant = new VCFVariant(next);
				if (!variant.equals(testSNPObj)) {
					// correlate
					if (variant.getMAF() > 0.005 && variant.getAlleles().length == 2) {
						double[] genotypes2 = correlator2.convertToDouble(variant);
						double corr = Correlation.correlate(genotypes1, genotypes2);
						corr *= corr;
						if (corr >= threshold) {
							output.add(f.getName() + "\t" + variant.getId() + "\t" + variant.getPos() + "\t" + (Math.abs(f.getStart() - variant.getPos())) + "\t" + corr);
						}

					}
				}
				next = window.next();
			}
			return output;
		}
	}

}
