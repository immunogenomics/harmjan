package nl.harmjanwestra.gwas;

import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.gwas.CLI.ProxyFinderOptions;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Correlation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.*;

/**
 * Created by hwestra on 2/29/16.
 */
public class ProxyFinder {

	ProxyFinderOptions options;

	public ProxyFinder(ProxyFinderOptions options) throws IOException {
		this.options = options;

		if (options.pairwise) {
			pairwiseLD();
		} else {
			findProxies();
		}

	}

	public void findProxies() throws IOException {
		ArrayList<Feature> snps = new ArrayList<Feature>();
		TextFile tf = new TextFile(options.snpfile, TextFile.R);
		String ln = tf.readLine();
		while (ln != null) {
			String[] elems = ln.split("\t");
			if (elems.length >= 4) {
				Feature snp = new Feature();
				snp.setStart(Integer.parseInt(elems[2]));
				snp.setStop(snp.getStart());
				snp.setName(elems[3]);
				snp.setChromosome(Chromosome.parseChr(elems[1]));
				snps.add(snp);
			}
			ln = tf.readLine();
		}
		tf.close();

		if (snps.isEmpty()) {
			System.err.println("Error: no snps in " + options.snpfile);
			System.exit(-1);
		} else {
			// check whether the reference VCFs are here..

			boolean allfilespresent = true;
			for (Feature queryVariantFeature : snps) {
				if (!Gpio.exists(options.tabixrefprefix + queryVariantFeature.getChromosome().getNumber() + ".vcf.gz")) {
					System.out.println("Could not find required file: " + options.tabixrefprefix + queryVariantFeature.getChromosome().getNumber() + ".vcf.gz");
					allfilespresent = false;
				}
			}

			if (!allfilespresent) {
				System.out.println("Some VCF files are missing");
				System.exit(-1);
			}
		}

		System.out.println(snps.size() + " snps in " + options.snpfile);

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
				System.out.println(returned + " / " + submit + " returned");
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}

		outFile.close();
		threadPool.shutdown();
	}

	public void pairwiseLD() throws IOException {
		ArrayList<Pair<String, String>> pairs = new ArrayList<>();
		TextFile tf = new TextFile(options.snpfile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length >= 2) {
				pairs.add(new Pair(elems[0], elems[1]));
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		// check if all VCF files are there
		if (pairs.isEmpty()) {
			System.err.println("Error: no snp pairs in " + options.snpfile);
			System.exit(-1);
		} else {
			// check whether the reference VCFs are here..

			boolean allfilespresent = true;
			for (Pair<String, String> p : pairs) {
				String snp1 = p.getLeft();
				String snp2 = p.getRight();
				String[] snp1elems = snp1.split("_");
				Chromosome chr1 = Chromosome.parseChr(snp1elems[0]);
				String[] snp2elems = snp2.split("_");
				Chromosome chr2 = Chromosome.parseChr(snp2elems[0]);

				if (!Gpio.exists(options.tabixrefprefix + chr1.getNumber() + ".vcf.gz")) {
					System.out.println("Could not find required file: " + options.tabixrefprefix + chr1.getNumber() + ".vcf.gz");
					allfilespresent = false;
				}

				if (!Gpio.exists(options.tabixrefprefix + chr2.getNumber() + ".vcf.gz")) {
					System.out.println("Could not find required file: " + options.tabixrefprefix + chr2.getNumber() + ".vcf.gz");
					allfilespresent = false;
				}
			}

			if (!allfilespresent) {
				System.out.println("Some VCF files are missing");
				System.exit(-1);
			}
		}

		TextFile out = new TextFile(options.output, TextFile.W);
		int nr = 0;
		for (Pair<String, String> p : pairs) {

			String snp1 = p.getLeft();
			String snp2 = p.getRight();
			String[] snp1elems = snp1.split("_");
			Integer snp1pos = Integer.parseInt(snp1elems[1]);
			Chromosome chr1 = Chromosome.parseChr(snp1elems[0]);
			String[] snp2elems = snp2.split("_");
			Chromosome chr2 = Chromosome.parseChr(snp2elems[0]);
			Integer snp2pos = Integer.parseInt(snp2elems[1]);

			Feature snpfeature1 = new Feature();
			snpfeature1.setChromosome(chr1);
			snpfeature1.setStart(snp1pos);
			snpfeature1.setStart(snp1pos);
			snpfeature1.setName(snp1elems[2]);

			Feature snpfeature2 = new Feature();
			snpfeature2.setChromosome(chr2);
			snpfeature2.setStart(snp2pos);
			snpfeature2.setStart(snp2pos);
			snpfeature2.setName(snp2elems[2]);

			VCFVariant variant1 = getSNP(snpfeature1);
			VCFVariant variant2 = getSNP(snpfeature2);

			if (variant1 != null && variant2 != null) {
				double[] genotypes1 = convertToDouble(variant1);
				double[] genotypes2 = convertToDouble(variant2);
				Pair<double[], double[]> filtered = stripmissing(genotypes1, genotypes2);
				double corr = Correlation.correlate(filtered.getLeft(), filtered.getRight());

				double rsq = (corr * corr);
				out.writeln(snp1 + "\t" + snp2 + "\t" + rsq);

				if (rsq > options.threshold) {
					nr++;
				}

			}
		}
		out.close();

		System.out.println(nr + " of " + pairs.size() + " have rsq>" + options.threshold);
	}

	private Pair<double[], double[]> stripmissing(double[] a, double[] b) {
		boolean[] missing = new boolean[a.length];
		int nrmissing = 0;
		for (int i = 0; i < a.length; i++) {
			if (a[i] == -1 || b[i] == -1) {
				missing[i] = true;
				nrmissing++;
			}
		}
		if (nrmissing == 0) {
			return new Pair<double[], double[]>(a, b);
		} else {
			double[] tmpa = new double[a.length - nrmissing];
			double[] tmpb = new double[a.length - nrmissing];
			int ctr = 0;
			for (int i = 0; i < a.length; i++) {
				if (!missing[i]) {
					tmpa[ctr] = a[i];
					tmpb[ctr] = b[i];
					ctr++;
				}
			}

			return new Pair<double[], double[]>(tmpa, tmpb);
		}
	}

	private synchronized VCFVariant getSNP(Feature snp) throws IOException {
		String tabixfile = options.tabixrefprefix + snp.getChromosome().getNumber() + ".vcf.gz";
		TabixReader reader = new TabixReader(tabixfile);
		TabixReader.Iterator inputSNPiter = reader.query(snp.getChromosome().getNumber() + ":" + (snp.getStart() - 1) + "-" + (snp.getStart() + 1));
		String snpStr = inputSNPiter.next();
		VCFVariant testSNPObj1 = null;
//		System.out.println("Query: " + snp.toString());
		while (snpStr != null) {
			VCFVariant variant = new VCFVariant(snpStr);
			if (variant.asFeature().overlaps(snp)) {
				if (variant.getId().equals(snp.getName())) {
					testSNPObj1 = variant;
				}
			}

			snpStr = inputSNPiter.next();
		}
		reader.close();
		return testSNPObj1;
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
			String selfStr = queryVariantFeature.getChromosome().toString()
					+ "\t" + queryVariantFeature.getStart()
					+ "\t" + queryVariantFeature.getName()
					+ "\t" + queryVariantFeature.getChromosome().toString()
					+ "\t" + queryVariantFeature.getStart()
					+ "\t" + queryVariantFeature.getName()
					+ "\t" + 0
					+ "\t" + 1;

			VCFVariant testSNPObj = getSNP(queryVariantFeature);

			if (testSNPObj == null) {
				System.out.println("Query: " + queryVariantFeature.getName() + " not in reference. " + queryVariantFeature.getChromosome().getNumber() + ":" + (queryVariantFeature.getStart() - 1) + "-" + (queryVariantFeature.getStart() + 1));
				if (output.isEmpty()) {
					output.add(selfStr);
				}
				return output;
			}

			String tabixfile = tabixrefprefix + queryVariantFeature.getChromosome().getNumber() + ".vcf.gz";
			TabixReader reader = new TabixReader(tabixfile);
			double[] genotypes1 = convertToDouble(testSNPObj);
			TabixReader.Iterator window = reader.query(queryVariantFeature.getChromosome().getNumber() + ":" + (queryVariantFeature.getStart() - windowsize) + "-" + (queryVariantFeature.getStart() + windowsize));
			String next = window.next();

			while (next != null) {
				// correlate
				VCFVariant variant = new VCFVariant(next);
				if (!variant.equals(testSNPObj)) {
					// correlate
					if (variant.getMAF() > 0.005 && variant.getAlleles().length == 2) {
						double[] genotypes2 = convertToDouble(variant);

						Pair<double[], double[]> filtered = stripmissing(genotypes1, genotypes2);
						double corr = Correlation.correlate(filtered.getLeft(), filtered.getRight());
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
			reader.close();
			if (output.isEmpty()) {
				output.add(selfStr);
			}
			return output;
		}
	}
}
