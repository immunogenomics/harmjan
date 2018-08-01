import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.utilities.vcf.GeneticSimilarity;
import nl.harmjanwestra.utilities.vcf.SampleAnnotation;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Triple;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Created by hwestra on 4/12/17.
 */
public class Main {


	public static void main(String[] args) {

		Main m = new Main();
		try {

			String dir = "/Data/GGD/HapMap2r24CEU/";


			if (exists(dir)) {
				System.out.println("Exists!");
			}

			m.run(dir);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static boolean exists(String dir) {
		return existsAndReadable(new File(dir));
	}

	public static boolean existsAndReadable(File file) {
		return file.exists() && file.canRead();
	}


	public void run(String dsloc) throws IOException {


		int chr = 22;
		TriTyperGenotypeData ds = new TriTyperGenotypeData(dsloc);
		ArrayList<SNP> snps = new ArrayList<SNP>();

		{
			String[] snpNames = ds.getSNPs();
			SNPLoader loader = ds.createSNPLoader();
			ArrayList<Integer> snpsToLoad = new ArrayList<Integer>();
			for (int i = 0; i < snpNames.length; i++) {
				int snpchr = ds.getChr(i);
				if (snpchr == chr) {
					snpsToLoad.add(i);
				}
			}
			System.out.println(snpsToLoad.size() + " snps to load...");
			ArrayList<Pair<Integer, SNP>> snpsort = new ArrayList<Pair<Integer, SNP>>();
			for (int i = 0; i < snpsToLoad.size(); i++) {
				int snpid = snpsToLoad.get(i);
				SNP snp = ds.getSNPObject(snpid);
				snpsort.add(new Pair<Integer, SNP>(snp.getChrPos(), snp, Pair.SORTBY.LEFT));
			}
			Collections.sort(snpsort);

			for (int i = 0; i < snpsort.size(); i++) {
				Pair<Integer, SNP> p = snpsort.get(i);
				SNP s = p.getRight();
				loader.loadGenotypes(s);
				double maf = s.getMAF();
				double cr = s.getCR();
				double hwep = s.getHWEP();
				if (maf > 0.01 && cr > 0.95 && hwep > 0.001) {
					snps.add(p.getRight());
				}
				if (i % 1000 == 0) {
					System.out.print("Loaded: " + i + " snps.\r");
				}
			}
			loader.close();
		}

		System.out.print("\n");
		System.out.println("Done loading..");
		// snps loaded and sorted..
		// now compare IBS in windows of ~1000 variants, with 250 overlap
		int windowstart = 0;
		int windowstop = 1000;
		int stepsize = 250;
		boolean finished = false;
		SNP[] window = new SNP[windowstop - windowstart];
		GeneticSimilarity sim = new GeneticSimilarity();

		ExecutorService exService = Executors.newWorkStealingPool(Runtime.getRuntime().availableProcessors());


		while (!finished) {


			int ctr = 0;
			for (int i = windowstart; i < windowstop && i < snps.size(); i++) {
				window[ctr] = snps.get(i);
				ctr++;
			}

			// null padding
			for (int i = ctr; i < window.length; i++) {
				window[ctr] = null;
				ctr++;
			}

			// do some calculations

			VCFVariant[] variants = convertToVCFVariants(window);
			Triple<DoubleMatrix2D, DoubleMatrix2D, DoubleMatrix2D> similarityData = sim.calculateWithinDataset(variants, exService);

			// geneticSimilarity, geneticSimilaritySameGenotypes, geneticSimilarityCorrelation


			windowstart += stepsize;
			windowstop += stepsize;
			if (windowstart > snps.size()) {
				finished = true;
			}
			System.out.print("Window: " + windowstart + "\t" + windowstop + "\t" + snps.size());
		}

		exService.shutdown();

	}

	private VCFVariant[] convertToVCFVariants(SNP[] window) {

		int nrNotNull = 0;
		for (int i = 0; i < window.length; i++) {
			if (window[i] != null) {
				nrNotNull++;
			}
		}
		VCFVariant[] output = new VCFVariant[nrNotNull];
		DoubleMatrix2D dosages = null;
		SampleAnnotation sampleAnnotation = null;

		int ctr = 0;
		for (int i = 0; i < window.length; i++) {
			if (window[i] != null) {
				// String chr, Integer pos, String id, String alleleStr, String info, DoubleMatrix2D alleles, DoubleMatrix2D dosages, SampleAnnotation annotation
				SNP snp = window[i];
				String chr = new String("" + snp.getChr()).intern();
				Integer pos = snp.getChrPos();
				String id = snp.getName();
				String alleleStr = new String(snp.getAlleles()[0] + "," + snp.getAlleles()[1]).intern();
				String info = ".";
				DoubleMatrix2D allelesMat = DoubleFactory2D.dense.make(snp.getAllele1().length, 2, -1);// new DenseDoubleMatrix2D(nrSamples, 2);
				byte[] genotypes = snp.getGenotypes();
				for (int j = 0; j < snp.getAllele1().length; j++) {
					if (genotypes[j] != -1) {
						if (genotypes[j] == 2) {
							allelesMat.setQuick(j, 0, 1);
							allelesMat.setQuick(j, 1, 1);
						} else if (genotypes[j] == 1) {
							allelesMat.setQuick(j, 0, 0);
							allelesMat.setQuick(j, 1, 1);
						} else {
							allelesMat.setQuick(j, 0, 0);
							allelesMat.setQuick(j, 1, 0);
						}
					}

				}
				VCFVariant variant = new VCFVariant(chr, pos, id, alleleStr, info, allelesMat, dosages, sampleAnnotation);
				output[ctr] = variant;
				ctr++;
			}
		}
		return new VCFVariant[0];
	}
}
