package nl.harmjanwestra.utilities.vcf;

import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.vcf.filter.variantfilters.VCFVariantFilters;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.*;

/**
 * Created by hwestra on 3/21/17.
 */
public class VCFVariantLoader {

	public ArrayList<VCFVariant> run(String vcf, ArrayList<Feature> regions) throws IOException {
		return run(vcf, regions, null);
	}

	public ArrayList<VCFVariant> run(String vcf, ArrayList<Feature> regions, VCFVariantFilters filters) throws IOException {
		System.out.println("VCF Loader.\n-------------------------\nBooting up threadpool for " + Runtime.getRuntime().availableProcessors() + " CPUs");
		System.out.println("Loading: " + vcf);
		if (regions != null) {
			System.out.println("Filtering for: " + regions.size() + " regions");
		}
		if (filters != null) {
			System.out.println(filters.toString());
		}

		ExecutorService exService = Executors.newWorkStealingPool(Runtime.getRuntime().availableProcessors());
		CompletionService<ArrayList<VCFVariant>> queuee = new ExecutorCompletionService<ArrayList<VCFVariant>>(exService);

		TextFile tf = new TextFile(vcf, TextFile.R);
		String ln = tf.readLine();
		int buffersize = 25;
		int maxLinesInQueue = buffersize * Runtime.getRuntime().availableProcessors();

		String[] buffer = new String[buffersize];
		int ctr = 0;
		ArrayList<VCFVariant> output = new ArrayList<>();
		int submitted = 0;

		while (ln != null) {
			buffer[ctr] = ln;
			ctr++;
			if (ctr == buffersize) {
				queuee.submit(new Parser(buffer, regions, filters));
				ctr = 0;
				submitted++;
				buffer = new String[buffersize];

				System.out.print(submitted + " jobs submitted (" + (submitted * buffersize) + "lines)\t" + output.size() + " variants in buffer so far\r");
				if (submitted * buffersize == maxLinesInQueue) {
					cleanup(output, queuee, submitted);
					submitted = 0;
				}
			}
			ln = tf.readLine();
		}
		tf.close();

		if (ctr > 0) {
			queuee.submit(new Parser(buffer, regions, filters));
			submitted++;
		}

		cleanup(output, queuee, submitted);
		submitted = 0;

		System.out.println();
		System.out.println("Done. " + output.size() + " variants loaded.");

		exService.shutdown();
		return output;
	}

	private void cleanup(ArrayList<VCFVariant> buffer, CompletionService<ArrayList<VCFVariant>> queuee, int toReceive) {

		int received = 0;
		while (received < toReceive) {
			try {
				ArrayList<VCFVariant> future = queuee.take().get();
				if (future != null) {
					buffer.addAll(future);
					received++;
				}
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}
	}

	public class Parser implements Callable<ArrayList<VCFVariant>> {

		String[] buffer;
		ArrayList<Feature> regions;
		VCFVariantFilters filters;

		public Parser(String[] buffer, ArrayList<Feature> regions, VCFVariantFilters filters) {
			this.buffer = buffer;
			this.filters = filters;
			this.regions = regions;
		}

		@Override
		public ArrayList<VCFVariant> call() throws Exception {
			ArrayList<VCFVariant> output = new ArrayList<>();

			for (String s : buffer) {
				if (s != null && !s.startsWith("#")) {
					String substr = s.substring(0, 200);

					boolean parseln = false;
					if (regions == null) {
						parseln = true;
					} else {
						VCFVariant v = new VCFVariant(substr, VCFVariant.PARSE.HEADER);
						parseln = v.asFeature().overlaps(regions);
					}

					if (parseln) {
						VCFVariant v = new VCFVariant(s);
						if (filters != null) {
							if (filters.passesFilters(v)) {
								output.add(v);
							}
						} else {
							output.add(v);
						}

					}
				}
			}
			return output;
		}
	}

}
