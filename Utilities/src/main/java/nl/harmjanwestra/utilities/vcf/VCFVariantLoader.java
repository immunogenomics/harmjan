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

	boolean[] genotypeSamplesWithCovariatesAndDiseaseStatus;
	SampleAnnotation sampleAnnotation;

	public VCFVariantLoader(boolean[] genotypeSamplesWithCovariatesAndDiseaseStatus, SampleAnnotation sampleAnnotation, ExecutorService exService) {
		this.genotypeSamplesWithCovariatesAndDiseaseStatus = genotypeSamplesWithCovariatesAndDiseaseStatus;
		this.sampleAnnotation = sampleAnnotation;
	}

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

		boolean sharedExService = false;
		ExecutorService exService = Executors.newWorkStealingPool(Runtime.getRuntime().availableProcessors());
		CompletionService<ArrayList<VCFVariant>> queuee = new ExecutorCompletionService<ArrayList<VCFVariant>>(exService);


		KillSwitchEngage killswitch = new KillSwitchEngage();
		Collector c = new Collector(queuee, killswitch);
		Future<ArrayList<VCFVariant>> outputFuture = exService.submit(c);

		String[] files = new String[1];
		files[0] = vcf;
		if (vcf.contains("CHR")) {
			files = new String[22];
			for (int i = 1; i < 23; i++) {
				files[i - 1] = vcf.replaceAll("CHR", "" + i);
			}
		}

		int buffersize = 50;
		String[] buffer = new String[buffersize];
		int ctr = 0;
		int submitted = 0;
		int maxLinesInQueue = buffersize * Runtime.getRuntime().availableProcessors();

		for (String file : files) {
			TextFile tf = new TextFile(file, TextFile.R);
			String ln = tf.readLine();
			while (ln != null) {
				buffer[ctr] = ln;
				ctr++;
				if (ctr == buffersize) {
					queuee.submit(new Parser(buffer, regions, filters));
					ctr = 0;
					submitted++;
					buffer = new String[buffersize];
					System.out.print(submitted + " jobs submitted (" + (submitted * buffersize) + "lines)\t" + killswitch.getNrVariantsReceived() + " variants in buffer so far\r");
//				if (submitted * buffersize == maxLinesInQueue) {
//					cleanup(output, queuee, submitted);
//					submitted = 0;
//				}
				}
				ln = tf.readLine();
			}
			tf.close();
		}

		System.out.println("Done reading file. Waiting for output to return.");
		if (ctr > 0) {
			System.out.println("Submitting last batch");
			queuee.submit(new Parser(buffer, regions, filters));
			submitted++;
		}


		System.out.println("Switching the killswitch");
		killswitch.setNrSubmitted(submitted);
		killswitch.setFinishedLoading(true);

//		cleanup(output, queuee, submitted);
		System.out.println("Waiting for output");
		ArrayList<VCFVariant> output = null;
		try {
			output = outputFuture.get();
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		}


		System.out.println();
		System.out.println("Done. " + output.size() + " variants loaded.");

		if (!sharedExService) {
			exService.shutdown();
		}
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

	public class KillSwitchEngage {
		private int nrSubmitted = 0;
		private boolean finishedLoading = false;

		private int nrVariantsReceived = 0;

		public synchronized void setNrSubmitted(int nr) {
			this.nrSubmitted = nr;
		}

		public synchronized void setNrVariantsReceived(int n) {
			this.nrVariantsReceived = n;
		}

		public synchronized int getNrVariantsReceived() {
			return this.nrVariantsReceived;
		}

		public synchronized void setFinishedLoading(boolean b) {
			this.finishedLoading = b;
		}

		public synchronized boolean getFinishedLoading() {
			return finishedLoading;
		}

		public synchronized int getNrSubmitted() {
			return nrSubmitted;
		}
	}

	public class Collector implements Callable<ArrayList<VCFVariant>> {

		private KillSwitchEngage killswitch;
		private CompletionService<ArrayList<VCFVariant>> queuee;

		public Collector(CompletionService<ArrayList<VCFVariant>> queuee, KillSwitchEngage g) {
			this.killswitch = g;
			this.queuee = queuee;
		}

		@Override
		public ArrayList<VCFVariant> call() throws Exception {
			int received = 0;
			ArrayList<VCFVariant> buffer = new ArrayList<>();
			while (!killswitch.getFinishedLoading()) {
				try {
					ArrayList<VCFVariant> future = queuee.take().get();
					if (future != null) {
						buffer.addAll(future);
						killswitch.setNrVariantsReceived(buffer.size());
						received++;
					}
				} catch (InterruptedException e) {
					e.printStackTrace();
				} catch (ExecutionException e) {
					e.printStackTrace();
				}
			}

			// check whether there are any packages left in the queue that we should wait for
			while (received != killswitch.getNrSubmitted()) {
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

			return buffer;
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
						VCFVariant v = new VCFVariant(s, VCFVariant.PARSE.ALL, genotypeSamplesWithCovariatesAndDiseaseStatus, sampleAnnotation);
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
