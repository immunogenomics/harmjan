package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.*;

/**
 * Created by hwestra on 4/22/16.
 */
public class VCFBedFilter {


	public void filter(String vcfIn, String referenceOut, String regionFile, String chr, boolean list) throws IOException {
		// boot up threadpool


		if (list) {
//			TextFile tf = new TextFile(vcfIn, TextFile.R);
//			ArrayList<String[]> listArr = new ArrayList<>();
//			String[] elems = tf.readLineElems(TextFile.tab);
//			if (elems.length < 3) {
//				System.out.println("Expecting three columns in input");
//				System.exit(-1);
//			}
//			while (elems != null) {
//				listArr.add(elems);
//				elems = tf.readLineElems(TextFile.tab);
//			}
//			tf.close();
//
//
//			int cores = Runtime.getRuntime().availableProcessors();
//			System.out.println("Detected " + cores + " Processors ");
//			ExecutorService threadPool = Executors.newFixedThreadPool(cores);
//			CompletionService<String> jobHandler = new ExecutorCompletionService<>(threadPool);
//
//			int submitted = 0;
//			for (String[] f : listArr) {
//				if (f.length == 3) {
//					VCFBedFilterTask task = new VCFBedFilterTask(f[0], f[1], regionFile, f[2]);
//					jobHandler.submit(task);
//					submitted++;
//				}
//			}
//
//			System.out.println("Submitted " + submitted + " jobs");
//
//			int returned = 0;
//			while (returned < submitted) {
//
//
//				Future<String> future = null;
//				try {
//					future = jobHandler.take();
//					if (future != null) {
//						String out = future.get();
////						System.out.println(out);
//						returned++;
//						System.out.println(returned + " / " + submitted + " returned");
//					}
//				} catch (InterruptedException e) {
//					e.printStackTrace();
//				} catch (ExecutionException e) {
//					e.printStackTrace();
//				}
//			}
//
//			if (!threadPool.isShutdown()) {
//				threadPool.shutdown();
//			}

		} else {
			filterFile(vcfIn, referenceOut, regionFile, chr);
		}
	}

	public void filterFile(String vcfIn, String referenceOut, String regionFile, String chr) throws IOException {
		System.out.println("BED filter");
		System.out.println("in: " + vcfIn);
		System.out.println("out: " + referenceOut);
		System.out.println("bed: " + regionFile);

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> set = reader.readAsList(regionFile);
		if (chr != null) {
			ArrayList<Feature> tmp = new ArrayList<>();
			Chromosome chrom = Chromosome.parseChr(chr);
			for (Feature f : set) {
				if (f.getChromosome().equals(chrom)) {
					tmp.add(f);
				}
			}
			set = tmp;
			System.out.println(set.size() + " regions overlapping chr: " + chr);
		}


		TextFile vcfoutf = new TextFile(referenceOut, TextFile.W);
		TextFile vcfoutf2 = new TextFile(referenceOut + "-filterdOutVariants.txt.gz", TextFile.W);
		TextFile vcftf = new TextFile(vcfIn, TextFile.R);
		int nrHeaderElems = 9;
		String ln = vcftf.readLine();

		int parsed = 0;
		int saved = 0;

		int cores = Runtime.getRuntime().availableProcessors();
		System.out.println("Detected " + cores + " Processors ");
		ExecutorService threadPool = Executors.newFixedThreadPool(cores);
		CompletionService<String[]> jobHandler = new ExecutorCompletionService<>(threadPool);

		int submitted = 0;
		int returned = 0;
		int lnnum = 0;
		int written = 0;
		while (ln != null) {
			if (ln.startsWith("##")) {
				vcfoutf.writeln(ln);
			} else if (ln.startsWith("#CHROM")) {
				vcfoutf.writeln(ln);
			} else {

				VCFBedFilterTask t = new VCFBedFilterTask(ln, set);
				jobHandler.submit(t);
				submitted++;
				if (submitted % 10000 == 0) {
					while (returned != submitted) {
						try {
							Future<String[]> f = jobHandler.take();
							if (f != null) {
								String[] output = f.get();
								if (output[0] != null) {
									vcfoutf.writeln(output[0]);
									written++;
								} else {
									vcfoutf2.writeln(output[1]);
								}
							}
							returned++;
						} catch (InterruptedException e) {
							e.printStackTrace();
						} catch (ExecutionException e) {
							e.printStackTrace();
						}
					}
					System.out.print(lnnum + " lines parsed. " + submitted + " submitted. " + written + " written. " + returned + " returned \r");
				}
			}
			ln = vcftf.readLine();
			lnnum++;
		}

		while (returned != submitted) {
			try {
				Future<String[]> f = jobHandler.take();
				if (f != null) {
					String[] output = f.get();
					if (output[0] != null) {
						vcfoutf.writeln(output[0]);
						written++;
					} else {
						vcfoutf2.writeln(output[1]);
					}
				}
				returned++;
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}
		System.out.print(lnnum + " lines parsed. " + submitted + " submitted. " + written + " written. " + returned + " returned \r");

		vcfoutf.close();
		vcfoutf2.close();
		vcftf.close();
		threadPool.shutdown();

	}

	public class VCFBedFilterTask implements Callable<String[]> {

		private final String ln;
		private final ArrayList<Feature> set;

		public VCFBedFilterTask(String ln, ArrayList<Feature> set) {
			this.ln = ln;
			this.set = set;
		}

		@Override
		public String[] call() throws Exception {
			String header = ln.substring(0, 500);
			String[] headerElems = header.split("\t");

			String[] output = new String[2];
			Feature f = new Feature();
			f.setChromosome(Chromosome.parseChr(headerElems[0]));
			f.setStart(Integer.parseInt(headerElems[1]));
			f.setStop(Integer.parseInt(headerElems[1]));

			for (Feature region : set) {
				if (f.overlaps(region)) {
					output[0] = ln;
					return output;
				}
			}

			output[1] = Strings.concat(headerElems, Strings.tab, 0, 9);

			return output;
		}
	}


}

