package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.*;

/**
 * Created by Harm-Jan on 07/10/16.
 */
public class VCFVariantCounter {


	public static void main(String[] args) {

	}

	public void calc(String file, String bedfile) throws IOException {
		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bedfile);
		TextFile tf = new TextFile(file, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		VCFParserTask t = new VCFParserTask(regions, "#");
		int nr = 0;
		while (elems != null) {

			if (!elems[0].startsWith("#")) {
				int pos = Integer.parseInt(elems[1]);
				Feature f = new Feature(Chromosome.parseChr(elems[0]), pos, pos + 1);
				if (t.overlapregion(f)) {
					nr++;
				}
			}

			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

	}

	public void run(String fileStr, String bedfile) throws IOException {

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bedfile);
		String[] files = fileStr.split(",");

		int nrAboveMAFThreshold = 0;
		int nrAboveInfoThreshold = 0;
		int nrAboveInfoAndMafThreshold = 0;
		int nrTotal = 0;

		double mafthreshold = 0.01;
		double infothreshold = 0.8;

		int submitted = 0;
		int returned = 0;


		int cores = Runtime.getRuntime().availableProcessors();
		System.out.println("Detected " + cores + " Processors ");
		ExecutorService threadPool = Executors.newFixedThreadPool(cores);
		CompletionService<Pair<Double, Double>> jobHandler = new ExecutorCompletionService<>(threadPool);

		for (String f : files) {

			TextFile tf = new TextFile(f, TextFile.R);

			String ln = tf.readLine();
			while (ln != null) {

				VCFParserTask t = new VCFParserTask(regions, ln);

				jobHandler.submit(t);
				submitted++;

				if (submitted % 10000 == 0) {

					System.out.println(submitted + " lines submitted..");
					while (returned < submitted) {
						try {
							Future<Pair<Double, Double>> future = jobHandler.take();
							Pair<Double, Double> content = future.get();
							if (content != null) {
								double maf = content.getLeft();
								double info = content.getRight();
								if (maf > mafthreshold) {
									nrAboveMAFThreshold++;
								}
								if (info > infothreshold) {
									nrAboveInfoThreshold++;
								}
								if (info > infothreshold && maf > mafthreshold) {
									nrAboveInfoAndMafThreshold++;
								}
								nrTotal++;
							}
							returned++;
						} catch (InterruptedException e) {
							e.printStackTrace();
						} catch (ExecutionException e) {
							e.printStackTrace();
						}
					}

				}


				ln = tf.readLine();

			}

			tf.close();

		}

		while (returned < submitted) {
			try {
				Future<Pair<Double, Double>> future = jobHandler.take();
				Pair<Double, Double> content = future.get();
				if (content != null) {
					double maf = content.getLeft();
					double info = content.getRight();
					if (maf > mafthreshold) {
						nrAboveMAFThreshold++;
					}
					if (info > infothreshold) {
						nrAboveInfoThreshold++;
					}
					if (info > infothreshold && maf > mafthreshold) {
						nrAboveInfoAndMafThreshold++;
					}
					nrTotal++;
				}
				returned++;
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}


		System.out.println("Total: " + nrTotal);
		System.out.println("MAF threshold: " + nrAboveMAFThreshold);
		System.out.println("INFO threshold: " + nrAboveInfoThreshold);
		System.out.println("MAF and INFO threshold: " + nrAboveInfoAndMafThreshold);

		threadPool.shutdown();

	}

	public class VCFParserTask implements Callable<Pair<Double, Double>> {

		String ln;
		ArrayList<Feature> regions;

		public VCFParserTask(ArrayList<Feature> regions, String ln) {
			this.ln = ln;
			this.regions = regions;
		}

		public VCFParserTask() {

		}

		@Override
		public Pair<Double, Double> call() throws Exception {

			if (ln.startsWith("#")) {
				return null;
			}

			String substr = ln.substring(0, 200);
			String[] elems = substr.split("\t");
			Chromosome chr = Chromosome.parseChr(elems[0]);
			Integer pos = Integer.parseInt(elems[1]);
			Feature f = new Feature(chr, pos, pos + 1);

			if (!overlapregion(f)) {
				return null;
			}

			VCFVariant var = new VCFVariant(ln);
			double maf = var.getMAF();
			Double info = var.getImputationQualityScore();
			if (info == null) {
				info = 0d;
			}


			return new Pair<Double, Double>(maf, info);
		}

		boolean overlapregion(Feature f) {
			for (Feature r : regions) {
				if (r.overlaps(f)) {
					return true;
				}
			}
			return false;
		}
	}

}
