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
			TextFile tf = new TextFile(vcfIn, TextFile.R);
			ArrayList<String[]> listArr = new ArrayList<>();
			String[] elems = tf.readLineElems(TextFile.tab);
			if (elems.length < 3) {
				System.out.println("Expecting three columns in input");
				System.exit(-1);
			}
			while (elems != null) {
				listArr.add(elems);
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();


			int cores = Runtime.getRuntime().availableProcessors();
			System.out.println("Detected " + cores + " Processors ");
			ExecutorService threadPool = Executors.newFixedThreadPool(cores);
			CompletionService<String> jobHandler = new ExecutorCompletionService<>(threadPool);

			int submitted = 0;
			for (String[] f : listArr) {
				if(f.length == 3) {
					VCFBedFilterTask task = new VCFBedFilterTask(f[0], f[1], regionFile, f[2]);
					jobHandler.submit(task);
					submitted++;
				}
			}

			System.out.println("Submitted " + submitted + " jobs");

			int returned = 0;
			while (returned < submitted) {


				Future<String> future = null;
				try {
					future = jobHandler.take();
					if (future != null) {
						String out = future.get();
//						System.out.println(out);
						returned++;
						System.out.println(returned + " / " + submitted + " returned");
					}
				} catch (InterruptedException e) {
					e.printStackTrace();
				} catch (ExecutionException e) {
					e.printStackTrace();
				}
			}

			if (!threadPool.isShutdown()) {
				threadPool.shutdown();
			}

		} else {
			VCFBedFilterTask task = new VCFBedFilterTask(vcfIn, referenceOut, regionFile, chr);
			task.call();
		}

	}

	public class VCFBedFilterTask implements Callable<String> {
		String vcfIn;
		String referenceOut;
		String regionFile;
		String chr;

		public VCFBedFilterTask(String vcfIn, String referenceOut, String regionFile, String chr) throws IOException {
			this.vcfIn = vcfIn;
			this.referenceOut = referenceOut;
			this.regionFile = regionFile;
			this.chr = chr;
		}

		@Override
		public String call() throws IOException {

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


			while (ln != null) {
				if (ln.startsWith("##")) {
					vcfoutf.writeln(ln);
				} else if (ln.startsWith("#CHROM")) {
					vcfoutf.writeln(ln);
				} else {

					String header = ln.substring(0,500);
					String[] headerElems = header.split("\t");

					Feature f = new Feature();
					f.setChromosome(Chromosome.parseChr(headerElems[0]));
					f.setStart(Integer.parseInt(headerElems[1]));
					f.setStop(Integer.parseInt(headerElems[1]));
					boolean overlap = false;
					for (Feature region : set) {
						if (f.overlaps(region)) {
							overlap = true;
							break;
						}
					}

					if (overlap) {
						vcfoutf.writeln(ln);
						saved++;
					} else {
						vcfoutf2.writeln(Strings.concat(headerElems, Strings.tab, 0, 9));
					}
					parsed++;

					if (parsed % 50000 == 0) {
						System.out.println(saved + "/" + parsed + " written");
					}
				}

				ln = vcftf.readLine();
			}

			vcfoutf.close();
			vcfoutf2.close();
			vcftf.close();

			return (saved + "/" + parsed + " written");

		}
	}
}

