package nl.harmjanwestra.gwas;

import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.gwas.CLI.ProxyFinderOptions;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Correlation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.*;

/**
 * Created by hwestra on 2/29/16.
 */
public class ProxyFinder {

	public ProxyFinder(ProxyFinderOptions options) throws IOException {

		ArrayList<Feature> snps = new ArrayList<Feature>();
		TextFile tf = new TextFile(options.snpfile, TextFile.R);
		String ln = tf.readLine();
		while (ln != null) {
			Feature f = Feature.parseFeature(ln);
			snps.add(f);
			ln = tf.readLine();
		}
		tf.close();

		ExecutorService threadPool = Executors.newFixedThreadPool(options.nrthreads);
		CompletionService<ArrayList<String>> jobHandler = new ExecutorCompletionService<>(threadPool);

		int submit = 0;
		for (Feature snp : snps) {
			ProxyFinderTask task = new ProxyFinderTask(options.tabixrefprefix, options.windowsize, snp, options.threshold);
			jobHandler.submit(task);
			submit++;
		}

		int returned = 0;
		TextFile outFile = new TextFile(options.output, TextFile.W);
		while (returned < submit) {
			try {
				ArrayList<String> proxies = jobHandler.take().get();
				if (proxies != null) {
					for (String s : proxies) {
						outFile.writeln(s);
					}
				}
				returned++;
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}

		}

		outFile.close();
	}

	public class ProxyFinderTask implements Callable<ArrayList<String>> {

		String tabixrefprefix;
		int windowsize;
		Feature queryVariantFeature;
		double threshold;

		public ProxyFinderTask(String tabix, int windowsize, Feature snp, double threshold) {
			this.tabixrefprefix = tabix;
			this.windowsize = windowsize;
			this.queryVariantFeature = snp;
			this.threshold = threshold;
		}

		@Override
		public ArrayList<String> call() throws Exception {

			ArrayList<String> output = new ArrayList<>();
			output.add(queryVariantFeature.getName() + "\t" + queryVariantFeature.getName() + "\t" + queryVariantFeature.getStart() + "\t" + 0 + "\t" + 1);
			TabixReader reader = new TabixReader(tabixrefprefix + queryVariantFeature.getChromosome().getNumber());


			TabixReader.Iterator inputSNPiter = reader.query(queryVariantFeature.getChromosome().toString(), queryVariantFeature.getStart() - 1, queryVariantFeature.getStart() + 1);

			String snpStr = inputSNPiter.next();
			VCFVariant testSNPObj = null;
			System.out.println("Query: " + queryVariantFeature.getName());
			while (snpStr != null) {
				VCFVariant variant = new VCFVariant(snpStr);
				if (variant.asFeature().overlaps(queryVariantFeature)) {
					if (variant.getId().equals(queryVariantFeature.getName())) {
						testSNPObj = variant;
					}
				}
				snpStr = inputSNPiter.next();
			}

			if (testSNPObj == null) {
				System.out.println("Query: " + queryVariantFeature.getName() + " not in reference");
				return output;
			}

			double[] genotypes1 = convertToDouble(testSNPObj);
			TabixReader.Iterator window = reader.query(queryVariantFeature.getChromosome().toString(), queryVariantFeature.getStart() - windowsize, queryVariantFeature.getStart() + windowsize);
			String next = window.next();
			while (next != null) {
				// correlate
				VCFVariant variant = new VCFVariant(next);
				if (!variant.equals(testSNPObj)) {
					// correlate
					if (variant.getMAF() > 0.005 && variant.getAlleles().length == 2) {
						double[] genotypes2 = convertToDouble(variant);
						double corr = Correlation.correlate(genotypes1, genotypes2);
						corr *= corr;

						// TODO: should replace with actual linkage
						if (corr >= threshold) {
							output.add(queryVariantFeature.getChromosome().toString()
									+ "\t" + queryVariantFeature.getStart()
									+ "\t" + queryVariantFeature.getName()
									+ "\t" + variant.asFeature().getChromosome().toString()
									+ "\t" + variant.asFeature().getStart()
									+ "\t" + variant.asFeature().getName()
									+ "\t" + (Math.abs(queryVariantFeature.getStart() - variant.getPos()))
									+ "\t" + corr);
						}
					}
				}
				next = window.next();
			}
			return output;
		}
	}

	public double[] convertToDouble(VCFVariant vcfVariant) {
		byte[][] alleles = vcfVariant.getGenotypeAlleles();
		double[] output = new double[alleles[0].length];
		for (int i = 0; i < alleles[0].length; i++) {
			if (alleles[0][i] == -1) {
				output[i] = -1;
			} else {
				output[i] = (alleles[0][i] + alleles[1][i]);
			}
		}

		return output;
	}
}
