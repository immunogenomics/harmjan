package nl.harmjanwestra.utilities.vcf;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.vcf.filter.variantfilters.VCFVariantFilter;
import nl.harmjanwestra.utilities.vcf.filter.variantfilters.VCFVariantFilters;
import nl.harmjanwestra.utilities.vcf.filter.variantfilters.VCFVariantMAFFilter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.concurrent.*;

/**
 * Created by hwestra on 3/21/17.
 */
public class VCFVariantLoader {
	
	boolean[] samplesToInclude;
	SampleAnnotation sampleAnnotation;
	
	public VCFVariantLoader() {
	
	}
	
	public VCFVariantLoader(boolean[] samplesToInclude, SampleAnnotation sampleAnnotation) {
		this.samplesToInclude = samplesToInclude;
		this.sampleAnnotation = sampleAnnotation;
	}
	
	public ArrayList<VCFVariant> run(String vcf) throws IOException {
		return run(vcf, null);
	}
	
	public ArrayList<VCFVariant> run(String vcf, VCFVariantFilters filters) throws IOException {
		return run(vcf, filters, -1, null);
	}
	
	public ArrayList<VCFVariant> run(String vcf, VCFVariantFilters filters, int nrThreads) throws IOException {
		return run(vcf, filters, nrThreads, null);
	}
	
	public ArrayList<VCFVariant> run(String vcf, VCFVariantFilters filters, int nrThreads, String samplesToIncludeStr) throws IOException {
		
		if (samplesToIncludeStr != null) {
			VCFTabix t = new VCFTabix();
			samplesToInclude = t.getSampleFilter(samplesToIncludeStr, vcf);
		}
		boolean sharedExService = false;
		int threadsToUse = Runtime.getRuntime().availableProcessors();
		if (nrThreads > 0) {
			threadsToUse = nrThreads;
		}
		if (threadsToUse < 2) {
			threadsToUse = 2;
		}
		
		
		//ExecutorService exService = Executors.newWorkStealingPool(threadsToUse);
		int buffersize = 100;
		int maxLinesInMemory = 250 * threadsToUse;
		int workQueueSize = maxLinesInMemory / buffersize;
		LinkedBlockingQueue<Runnable> queue = new LinkedBlockingQueue<Runnable>(workQueueSize);
		ExecutorService exService = new ThreadPoolExecutor(threadsToUse, threadsToUse, 0L, TimeUnit.MILLISECONDS, queue);
		System.out.println("VCF Loader.\n-------------------------\n" +
				"Booting up threadpool for " + threadsToUse + " CPUs");
		System.out.println("Work queue size: " + workQueueSize + " jobs of " + buffersize + " lines");
		System.out.println("Loading: " + vcf);
		
		if (filters != null) {
			System.out.println(filters.toString());
		}
		
		
		CompletionService<LinkedList<VCFVariant>> queuee = new ExecutorCompletionService<LinkedList<VCFVariant>>(exService);
		
		
		KillSwitchEngage killswitch = new KillSwitchEngage();
		Collector collector = new Collector(queuee, killswitch);
		Future<LinkedList<VCFVariant>> outputFuture = exService.submit(collector);
		
		String[] files = new String[1];
		files[0] = vcf;
		if (vcf.contains("CHR")) {
			files = new String[22];
			for (int i = 1; i < 23; i++) {
				files[i - 1] = vcf.replaceAll("CHR", "" + i);
			}
		}
		
		
		String[] buffer = new String[buffersize];
		int ctr = 0;
		int submitted = 0;
		
		
		for (String file : files) {
			TextFile tf = new TextFile(file, TextFile.R);
			String ln = tf.readLine();
			while (ln != null) {
				buffer[ctr] = ln;
				ctr++;
				if (ctr == buffersize) {
					while (queue.size() >= workQueueSize) {
						try {
//						System.out.println(queue.size() + " still in buffer out of " + workQueueSize);
							Thread.sleep(1);
						} catch (InterruptedException e) {
							e.printStackTrace();
						}
						
					}
					queuee.submit(new Parser(buffer, filters));
					ctr = 0;
					submitted++;
					buffer = new String[buffersize];
					System.out.print(submitted + " jobs submitted (" + (submitted * buffersize) + "lines)\t" + killswitch.getNrParsedVariants() + " variants in buffer so far, " + queue.size() + " jobs pending.\r");
//				if (submitted * buffersize == maxLinesInQueue) {
//					cleanup(output, queuee, submitted);
//					submitted = 0;
//				}
				}
				ln = tf.readLine();
			}
			tf.close();
		}
		System.out.println();
		System.out.println("Done reading file. Waiting for output to return.");
		if (ctr > 0) {
			while (queue.size() >= workQueueSize) {
				try {
//						System.out.println(queue.size() + " still in buffer out of " + workQueueSize);
					Thread.sleep(1);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
				
			}
			System.out.println("Submitting last batch");
			queuee.submit(new Parser(buffer, filters));
			submitted++;
		}
		
		
		System.out.println("Switching the killswitch");
		killswitch.setNrSubmitted(submitted);
		killswitch.setFinishedLoading(true);

//		cleanup(output, queuee, submitted);
		System.out.println("Waiting for output");
		LinkedList<VCFVariant> output = null;
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
		ArrayList<VCFVariant> v = new ArrayList<>(output);
		return v;
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
		private int nrReceived;
		private int nrParsedVariants;
		
		public synchronized void setNrSubmitted(int nr) {
			this.nrSubmitted = nr;
		}
		
		public synchronized void setNrReceived(int n) {
			this.nrReceived = n;
		}
		
		public synchronized int getNrReceived() {
			return this.nrReceived;
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
		
		
		public synchronized int getNrParsedVariants() {
			return nrParsedVariants;
		}
		
		public synchronized void setNrParsedVariants(int nrParsedVariants) {
			this.nrParsedVariants = nrParsedVariants;
		}
		
	}
	
	public class Collector implements Callable<LinkedList<VCFVariant>> {
		
		private KillSwitchEngage killswitch;
		private CompletionService<LinkedList<VCFVariant>> queuee;
		
		
		public Collector(CompletionService<LinkedList<VCFVariant>> queuee, KillSwitchEngage g) {
			this.killswitch = g;
			this.queuee = queuee;
		}
		
		@Override
		public LinkedList<VCFVariant> call() throws Exception {
			int received = 0;
			LinkedList<VCFVariant> buffer = new LinkedList<>();
			
			// loop until we've finished reading the VCF file
			while (!killswitch.getFinishedLoading()) {
				try {
					LinkedList<VCFVariant> future = queuee.take().get();
					
					if (future != null) {
						buffer.addAll(future);
						killswitch.setNrParsedVariants(buffer.size());
						received++;
						killswitch.setNrReceived(received);
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
					LinkedList<VCFVariant> future = queuee.take().get();
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
	
	public class Parser implements Callable<LinkedList<VCFVariant>> {
		
		String[] buffer;
		VCFVariantFilters filters;
		
		public Parser(String[] buffer, VCFVariantFilters filters) {
			this.buffer = buffer;
			this.filters = filters;
			
		}
		
		Feature tyk2 = new Feature(Chromosome.NINETEEN, 10396336, 10628468);
		
		@Override
		public LinkedList<VCFVariant> call() throws Exception {
			LinkedList<VCFVariant> output = new LinkedList<>();
			
			for (int i = 0; i < buffer.length; i++) {
				String s = buffer[i];
				if (s != null && !s.startsWith("#")) {
					
					
					boolean parseln = false;
					if (!filters.hasRegionOrVariantSetFilter()) {
						parseln = true;
					} else {
						String substr = s.substring(0, 200);
						VCFVariant v = new VCFVariant(substr, VCFVariant.PARSE.HEADER);
						parseln = filters.passesRegionOrVariantFilter(v);
					}
					
					if (parseln) {
						VCFVariant v = new VCFVariant(s, VCFVariant.PARSE.ALL, samplesToInclude, sampleAnnotation);
						if (v.asFeature().overlaps(tyk2)) {
							if (filters != null) {
								VCFVariantFilters filters2 = new VCFVariantFilters();
								ArrayList<VCFVariantFilter> setFilters = filters.getFilters();
								boolean maffilterset = false;
								for (VCFVariantFilter f : setFilters) {
									if (f instanceof VCFVariantMAFFilter) {
										maffilterset = true;
									} else {
										filters2.add(f);
									}
								}
								if (maffilterset) {
									filters2.addFilter(new VCFVariantMAFFilter(0.005));
								}
								if (filters2.passesFilters(v)) {
									output.add(v);
								} else {
									v.clean();
									v = null;
								}
							}
						} else if (filters != null) {
							if (filters.passesFilters(v)) {
								output.add(v);
							} else {
								v.clean();
								v = null;
							}
						} else {
							output.add(v);
						}
					}
				}
				s = null;
				buffer[i] = null;
			}
			buffer = null;
			return output;
		}
	}
	
}
