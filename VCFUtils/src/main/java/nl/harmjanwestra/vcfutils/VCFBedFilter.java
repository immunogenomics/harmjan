package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import nl.harmjanwestra.utilities.vcf.VCFVariant;

import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.*;

/**
 * Created by hwestra on 4/22/16.
 */
public class VCFBedFilter {


	public void filter(String vcfIn, String referenceOut, String regionFile, String chr, boolean list) throws IOException {
		// boot up threadpool


		if (vcfIn.contains("CHR")) {
			list = true;
		}
		if (list) {


			ArrayList<String> filelist = new ArrayList<String>();
			ArrayList<String> fileoutlist = new ArrayList<String>();
			for (int i = 1; i < 23; i++) {
				String f = vcfIn.replaceAll("CHR", "" + i);
				String f2 = referenceOut.replaceAll("CHR", "" + i);
				if (Gpio.exists(f)) {
					filelist.add(f);
					fileoutlist.add(f2);
				}
			}

			System.out.println(filelist.size() + " files found.");
			ExecutorService threadPool = Executors.newFixedThreadPool(filelist.size());
			for (int i = 0; i < filelist.size(); i++) {
				FilterTask t = new FilterTask(i, filelist.get(i), fileoutlist.get(i), regionFile, i + 1);
				threadPool.submit(t);
			}

			threadPool.shutdown();
			while (!threadPool.isTerminated()) {
			}
			System.out.println("Finished all threads");

		} else {
			FilterTask t = new FilterTask(0, vcfIn, referenceOut, regionFile, 0);
			t.run();
		}
	}

	class FilterTask implements Runnable {

		private int t;
		String vcfIn;
		String referenceOut;
		String regionFile;
		int chr;

		public FilterTask(int t, String in, String out, String region, int chr) throws IOException {
			this.t = t;
			this.vcfIn = in;
			this.referenceOut = out;
			this.regionFile = region;
			this.chr = chr;
		}

		@Override
		public void run() {
			try {
				System.out.println("BED filter");
				System.out.println("in: " + vcfIn);
				System.out.println("out: " + referenceOut);
				System.out.println("bed: " + regionFile);

				BedFileReader reader = new BedFileReader();
				ArrayList<Feature> set = reader.readAsList(regionFile);

				if (chr > 0) {
					ArrayList<Feature> settmp = new ArrayList<>();
					for (Feature region : set) {
						if (region.getChromosome().getNumber() == chr) {
							settmp.add(region);
						}
					}
					set = settmp;
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

						String substr = ln.substring(0, 1000);
						VCFVariant v = new VCFVariant(substr, VCFVariant.PARSE.HEADER);
						if (v.asFeature().overlaps(set)) {
							vcfoutf.writeln(ln);
							written++;
						}
//				VCFBedFilterTask t = new VCFBedFilterTask(ln, set);
//				jobHandler.submit(t);
//				submitted++;
//				if (submitted % 10000 == 0) {
//					while (returned != submitted) {
//						try {
//							Future<String[]> f = jobHandler.take();
//							if (f != null) {
//								String[] output = f.get();
//								if (output[0] != null) {
//									vcfoutf.writeln(output[0]);
//									written++;
//								} else {
//									vcfoutf2.writeln(output[1]);
//								}
//							}
//							returned++;
//						} catch (InterruptedException e) {
//							e.printStackTrace();
//						} catch (ExecutionException e) {
//							e.printStackTrace();
//						}
//					}

//				}
					}
					ln = vcftf.readLine();

					lnnum++;
					if (lnnum % 1000 == 0) {
						System.out.println("Task:\t" + t + "\t" + lnnum + " lines parsed. " + submitted + " submitted. " + written + " written. ");
					}
				}
				System.out.print("Task complete:\t" + t + "\t" + lnnum + " lines parsed. " + submitted + " submitted. " + written + " written. ");

				vcfoutf.close();
				vcfoutf2.close();
				vcftf.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
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

