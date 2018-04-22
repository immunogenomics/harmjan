package nl.harmjanwestra.proxyfinder;


import cern.colt.matrix.tdouble.DoubleMatrix2D;
import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.console.ProgressBar;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.Correlation;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFTabix;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.vcf.VCFVariantComparator;
import nl.harmjanwestra.utilities.vcf.VCFVariantLoader;
import nl.harmjanwestra.utilities.vcf.filter.variantfilters.VCFVariantFilters;
import nl.harmjanwestra.utilities.vcf.filter.variantfilters.VCFVariantSetFilter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.*;

/**
 * Created by hwestra on 2/29/16.
 */
public class ProxyFinder {
	
	ProxyFinderOptions options;
	
	public ProxyFinder(ProxyFinderOptions options) {
		this.options = options;
	}
	
	// use tabix to find proxies
	public void findProxies() throws IOException {
		System.out.println("Trying to find proxies for snps specified in :" + options.snpfile);
		ArrayList<SNPFeature> snps = new ArrayList<SNPFeature>();
		
		if (options.snpfile == null) {
			System.out.println("Please set snp file");
			System.exit(-1);
		} else if (!Gpio.exists(options.snpfile)) {
			System.out.println("Could not find file: " + options.snpfile);
			System.exit(-1);
		}
		
		TextFile tf = new TextFile(options.snpfile, TextFile.R);
		String ln = tf.readLine();
		if (ln == null) {
			System.out.println("Error with snp file. First line is null");
			System.exit(-1);
		}
		int chrCol = 0;
		int posCol = 1;
		int rsCol = 2;
		if (ln.trim().length() > 0) {
			String[] firstlnelems = Strings.whitespace.split(ln);
			if (firstlnelems.length >= 3 &&
					firstlnelems[0].toLowerCase().equals("snp") &&
					firstlnelems[1].toLowerCase().equals("chrom") &&
					firstlnelems[2].toLowerCase().equals("bp")) {
				// goshifter input file format
				System.out.println("File is in GoShifter format...");
				rsCol = 0;
				chrCol = 1;
				posCol = 2;
				ln = tf.readLine();
			} else {
				System.out.println("Assuming three column format");
			}
		}
		
		// format: chr pos rsid or chr_pos_rsid
		while (ln != null) {
			if (ln.trim().length() > 0) {
				String[] elems = ln.split("\t");
				if (elems.length >= 3) {
					SNPFeature snp = new SNPFeature();
					snp.setStart(Integer.parseInt(elems[posCol]));
					snp.setStop(snp.getStart());
					snp.setName(elems[rsCol]);
					snp.setChromosome(Chromosome.parseChr(elems[chrCol]));
					snps.add(snp);
				} else if (elems.length == 2) {
					// can't quite remember which format this was...
					SNPFeature snp = SNPFeature.parseSNPFeature(elems[1]);
					snps.add(snp);
				}
				ln = tf.readLine();
			}
		}
		tf.close();
		
		
		if (snps.isEmpty()) {
			System.err.println("Error: no snps in " + options.snpfile);
			System.exit(-1);
		} else {
			
			// check whether the reference VCFs are here..
			boolean allfilespresent = true;
			for (Feature queryVariantFeature : snps) {
				String tabixfile = options.tabixrefprefix.replaceAll("CHR", "" + queryVariantFeature.getChromosome().getNumber());
				//System.out.println("Trying to locate: " + tabixfile);
				if (!Gpio.exists(tabixfile)) {
					System.out.println("Could not find required path: " + tabixfile);
					allfilespresent = false;
				}
				// check whether the tabix indexes are there as well..
				if (!Gpio.exists(tabixfile + ".tbi")) {
					System.out.println("Could not find TABIX index path: " + options.tabixrefprefix + queryVariantFeature.getChromosome().getNumber() + ".vcf.gz.tbi");
					allfilespresent = false;
				}
			}
			
			if (!allfilespresent) {
				System.out.println("Some VCF files are missing");
				System.exit(-1);
			}
		}
		
		System.out.println(snps.size() + " snps in " + options.snpfile);
		if (options.nrthreads == 1) {
			System.out.println("Using: " + options.nrthreads + " thread");
		} else {
			System.out.println("Using: " + options.nrthreads + " threads");
		}
		
		ExecutorService threadPool = Executors.newFixedThreadPool(options.nrthreads);
		CompletionService<ArrayList<String>> jobHandler = new ExecutorCompletionService<>(threadPool);
		
		int submit = 0;
		for (Feature snp : snps) {
			if (snp != null) {
				ProxyFinderTask task = new ProxyFinderTask(snp);
				jobHandler.submit(task);
				submit++;
			}
		}
		
		System.out.println(snps.size() + " input variants. ");
		
		int returned = 0;
		TextFile outFile = new TextFile(options.output, TextFile.W);
		outFile.writeln("ChromA\tPosA\tRsIdA\tChromB\tPosB\tRsIdB\tDistance\tRSquared\tDprime");
		ProgressBar pb = new ProgressBar(submit);
		while (returned < submit) {
			try {
				ArrayList<String> proxies = jobHandler.take().get();
				if (proxies != null) {
					for (String s : proxies) {
						outFile.writeln(s);
					}
				}
				returned++;
				if (returned % 10 == 0) {
					System.out.println();
//					System.out.println(returned + " SNPs done out of " + submit + "\t" + ((double) returned / submit));
					
					pb.set(returned);
					pb.print();
					
					System.out.println();
				}
				
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}
//		pb.close();
		outFile.close();
		threadPool.shutdown();
		System.out.println("Done. Have a nice day, and don't forget that your output is here: " + outFile.getFullPath());
	}
	
	// calculate ld for all variants within a locus
	public void locusLD() throws IOException {
		
		BedFileReader br = new BedFileReader();
		if (options.regionfile == null) {
			System.out.println("Please set region file");
			System.exit(-1);
		} else if (!Gpio.exists(options.regionfile)) {
			System.out.println("Could not find region file: " + options.regionfile);
			System.exit(-1);
		}
		
		ArrayList<Feature> regions = br.readAsList(options.regionfile);
		System.out.println(regions.size() + " regions in " + options.regionfile);
		if (regions.isEmpty()) {
			System.out.println("No regions defined.");
			System.exit(-1);
		}
		// check whether the reference VCFs are here..
		
		boolean allfilespresent = true;
		for (Feature f : regions) {
			
			String tabixfile = options.tabixrefprefix.replaceAll("CHR", "" + f.getChromosome().getNumber());
			if (!Gpio.exists(tabixfile)) {
				System.out.println("Could not find required path: " + tabixfile);
				allfilespresent = false;
			}
		}
		
		System.out.println("All files are here. ");
		if (allfilespresent) {
			for (Feature f : regions) {
				System.out.println();
				System.out.println("Processing region: " + f.toString());
				ArrayList<VCFVariant> variants = new ArrayList<>();
				String tabixfile = options.tabixrefprefix.replaceAll("CHR", "" + f.getChromosome().getNumber());
				VCFTabix reader = new VCFTabix(tabixfile);
				boolean[] sampleSelect = reader.getSampleFilter(options.samplefilter);
				TabixReader.Iterator window = reader.query(f);
				String next = window.next();
				int ctr = 0;
				while (next != null) {
					// correlate
					VCFVariant variant = new VCFVariant(next, VCFVariant.PARSE.ALL, sampleSelect);
					if (variant.getMAF() > options.mafthreshold && variant.getAlleles().length == 2) {
						variant.calculateHWEP();
//						System.out.println(variant.getId() + "\t" + variant.getHwep());
						
						variants.add(variant);
					}
					ctr++;
					if (ctr % 100 == 0) {
						System.out.println(ctr + " variants parsed. " + variants.size() + " in memory");
					}
					next = window.next();
				}
				reader.close();
				
				
				System.out.println(variants.size() + " variants in region " + f.toString());
				
				
				System.out.println("Creating output: " + options.output + "-" + f.toString() + ".txt.gz");
				TextFile out = new TextFile(options.output + "-" + f.toString() + ".txt.gz", TextFile.W);
				
				// ChromA\tPosA\tRsIdA\tChromB\tPosB\tRsIdB\tDistance\tRSquared\tDprime
				out.writeln("ChromA\tPosA\tRsIdA\tAllelesA\tMinorAlleleA\tMAFA\tHWEPA\t" +
						"ChromB\tPosB\tRsIdB\tAllelesB\tMinorAlleleB\tMAFB\tHWEPB\t" +
						"Distance\tR-squared\tDprime");
				
				ProgressBar pb = new ProgressBar(variants.size(), f.toString());
				for (int i = 0; i < variants.size(); i++) {
					
					VCFVariant variant1 = variants.get(i);
					if (variant1.getAlleles().length == 2 && variant1.getMAF() > options.mafthreshold && variant1.getHwep() > options.hwepthreshold) {
						String variant1Str = variant1.getChr().toString() + "\t" + variant1.getPos() + "\t" + variant1.getId()
								+ "\t" + Strings.concat(variant1.getAlleles(), Strings.comma) + "\t" + variant1.getMinorAllele() + "\t" + variant1.getMAF() + "\t" + variant1.getHwep();
						for (int j = i + 1; j < variants.size(); j++) {
							DetermineLD ldcalc = new DetermineLD();
							VCFVariant variant2 = variants.get(j);
							if (variant2.getAlleles().length == 2 && variant2.getMAF() > options.mafthreshold && variant2.getHwep() > options.hwepthreshold) {
								Pair<Double, Double> ld = ldcalc.getLD(variant1, variant2);
								if (ld.getRight() > options.threshold) {
									String variant2Str = "\t" + variant2.getChr().toString() + "\t" + variant2.getPos() + "\t" + variant2.getId()
											+ "\t" + Strings.concat(variant2.getAlleles(), Strings.comma) + "\t" + variant2.getMinorAllele() + "\t" + variant2.getMAF() + "\t" + variant2.getHwep();
									
									variant2Str += "\t" + (variant2.getPos() - variant1.getPos()) + "\t" + ld.getRight() + "\t" + ld.getLeft();
									out.writeln(variant1Str + variant2Str);
								}
							}
						}
					}
					pb.set(i);
					
				}
				pb.close();
				out.close();
			}
			
		}
	}
	
	// use tabix for pairwise LD calculations
	public void pairwiseLD() throws IOException {
		
		if (options.snpfile == null) {
			System.out.println("Please set snp file");
			System.exit(-1);
		} else if (!Gpio.exists(options.snpfile)) {
			System.out.println("Could not find file: " + options.snpfile);
			System.exit(-1);
		}
		
		ArrayList<Pair<SNPFeature, SNPFeature>> pairs = new ArrayList<>();
		TextFile tf = new TextFile(options.snpfile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		ArrayList<SNPFeature> toCombine = null;
		while (elems != null) {
			if (elems.length == 1) {
				if (toCombine == null) {
					toCombine = new ArrayList<>();
				}
				if (elems[0].trim().length() > 0) {
					SNPFeature f = SNPFeature.parseSNPFeature(elems[0]);
					toCombine.add(f);
					System.out.println(f.toString());
				}
			} else if (elems.length == 2) {
				// format
				pairs.add(new Pair(SNPFeature.parseSNPFeature(elems[0]), SNPFeature.parseSNPFeature(elems[1])));
			} else if (elems.length == 3) {
				if (toCombine == null) {
					toCombine = new ArrayList<>();
				}
				SNPFeature snp1 = new SNPFeature();
				snp1.setStart(Integer.parseInt(elems[1]));
				snp1.setStop(snp1.getStart());
				snp1.setName(elems[2]);
				snp1.setChromosome(Chromosome.parseChr(elems[0]));
				toCombine.add(snp1);
			} else if (elems.length == 6) {
				SNPFeature snp1 = new SNPFeature();
				snp1.setStart(Integer.parseInt(elems[1]));
				snp1.setStop(snp1.getStart());
				snp1.setName(elems[2]);
				snp1.setChromosome(Chromosome.parseChr(elems[0]));
				
				SNPFeature snp2 = new SNPFeature();
				snp2.setStart(Integer.parseInt(elems[4]));
				snp2.setStop(snp2.getStart());
				snp2.setName(elems[5]);
				snp2.setChromosome(Chromosome.parseChr(elems[3]));
				
				pairs.add(new Pair<>(snp1, snp2));
				
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		if (toCombine != null) {
			System.out.println(toCombine.size() + " snps to combine.");
			for (int i = 0; i < toCombine.size(); i++) {
				for (int j = i + 1; j < toCombine.size(); j++) {
					Pair<SNPFeature, SNPFeature> p = new Pair(toCombine.get(i), toCombine.get(j));
					pairs.add(p);
					System.out.println(p.toString());
				}
			}
			
		}
		
		System.out.println(pairs.size() + " pairs loaded.");
		
		// check if all VCF files are there
		if (pairs.isEmpty()) {
			System.err.println("Error: no snp pairs in " + options.snpfile);
			System.exit(-1);
		} else {
			// check whether the reference VCFs are here..
			System.out.println(pairs.size() + " variant pairs loaded.");
			boolean allfilespresent = true;
			if (options.tabixrefprefix != null) {
				for (Pair<SNPFeature, SNPFeature> p : pairs) {
					
					SNPFeature snp1 = p.getLeft();
					SNPFeature snp2 = p.getRight();
					
					
					String tabixfile1 = options.tabixrefprefix.replaceAll("CHR", "" + snp1.getChromosome().getNumber());
					String tabixfile2 = options.tabixrefprefix.replaceAll("CHR", "" + snp2.getChromosome().getNumber());
					if (!Gpio.exists(tabixfile1)) {
						System.out.println("Could not find required path: " + tabixfile1);
						allfilespresent = false;
					}
					if (!Gpio.exists(tabixfile2)) {
						System.out.println("Could not find required path: " + tabixfile2);
						allfilespresent = false;
					}
				}
				
				if (!allfilespresent) {
					System.out.println("Some VCF files are missing");
					System.exit(-1);
				}
			} else {
				if (options.vcf == null || !Gpio.exists(options.vcf)) {
					System.out.println("Could not find input vcf: " + options.vcf);
					System.exit(-1);
				}
			}
		}
		
		int nr = 0;
		ArrayList<VCFVariant> allVariants = null;
		if (options.tabixrefprefix == null && options.vcf != null) {
			int nrthreads = Runtime.getRuntime().availableProcessors();
			if (options.nrthreads < nrthreads && options.nrthreads > 0) {
				nrthreads = options.nrthreads;
			}
			
			VCFVariantLoader loader = new VCFVariantLoader();
			VCFVariantFilters filter = new VCFVariantFilters();
			ArrayList<SNPFeature> snps = new ArrayList<>();
			for (Pair<SNPFeature, SNPFeature> p : pairs) {
				snps.add(p.getLeft());
				snps.add(p.getRight());
			}
			filter.addFilter(new VCFVariantSetFilter(snps));
			allVariants = loader.run(options.vcf, filter, nrthreads, options.samplefilter);
			Collections.sort(allVariants, new VCFVariantComparator());
		}
		
		TextFile out = new TextFile(options.output, TextFile.W);
		out.writeln("ChromA\tPosA\tRsIdA\tChromB\tPosB\tRsIdB\tDistance\tRSquared\tDprime");
		ProgressBar pb = new ProgressBar(pairs.size(), "Calculatinng LD...");
		for (Pair<SNPFeature, SNPFeature> p : pairs) {
			
			SNPFeature snpfeature1 = p.getLeft();
			SNPFeature snpfeature2 = p.getRight();
			
			VCFVariant variant1 = null;
			VCFVariant variant2 = null;
			
			if (options.tabixrefprefix != null) {
				variant1 = getSNP(snpfeature1);
				variant2 = getSNP(snpfeature2);
			} else {
				variant1 = getSNP(allVariants, snpfeature1);
				variant2 = getSNP(allVariants, snpfeature2);
			}
			
			if (variant1 != null && variant2 != null) {
				double[] genotypes1 = convertToDouble(variant1);
				double[] genotypes2 = convertToDouble(variant2);
				Pair<double[], double[]> filtered = stripmissing(genotypes1, genotypes2);
				double corr = Correlation.correlate(filtered.getLeft(), filtered.getRight());
				
				double rsq = (corr * corr);
				DetermineLD ldcalc = new DetermineLD();
				Pair<Double, Double> ldvals = ldcalc.getLD(variant1, variant2);
				
				String outStr = variant1.getChrObj().toString() + "\t" + variant1.getPos() + "\t" + variant1.getId()
						+ "\t" + variant2.getChrObj().toString() + "\t" + variant2.getPos() + "\t" + variant2.getId()
						+ "\t" + ldvals.getRight() + "\t" + ldvals.getLeft();
				
				out.writeln(outStr);
				
				if (rsq > options.threshold) {
					nr++;
				}
			} else {
				System.out.println();
				if (variant1 == null) {
					System.out.println(snpfeature1.toString() + " not found in reference");
				}
				if (variant2 == null) {
					System.out.println(snpfeature2.toString() + " not found in reference");
				}
			}
			pb.iterate();
		}
		pb.close();
		out.close();
		
		System.out.println(nr + " of " + pairs.size() + " have rsq>" + options.threshold);
	}
	
	private VCFVariant getSNP(ArrayList<VCFVariant> allVariants, SNPFeature snpfeature1) {
		for (VCFVariant v : allVariants) {
			if (v.asSNPFeature().getChromosome().equals(snpfeature1.getChromosome())) {
				if (v.getPos() == snpfeature1.getStart()) {
					return v;
				}
			}
		}
		return null;
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
		
		String tabixfile = options.tabixrefprefix.replace("CHR", "" + snp.getChromosome().getNumber());
		VCFTabix reader = new VCFTabix(tabixfile);
		boolean[] snpSampleFilter = reader.getSampleFilter(options.samplefilter);

//		System.out.println(tabixfile);
		
		Feature snpF = new Feature(snp);
		snpF.setStart(snpF.getStart() - 1);
		snpF.setStop(snpF.getStart() + 1);
		
		TabixReader.Iterator inputSNPiter = reader.query(snpF);
		String snpStr = inputSNPiter.next();
		VCFVariant testSNPObj1 = null;
//		System.out.println("Looking for: " + snp.toString());
		while (snpStr != null) {
			VCFVariant variant = null;
			try {
				variant = new VCFVariant(snpStr, VCFVariant.PARSE.ALL, snpSampleFilter);
//				System.out.println(variant.asFeature().toString());
			} catch (NullPointerException e) {
				String str = snpStr;
				if (str.length() > 1000) {
					str = str.substring(0, 1000);
				}
				System.out.println("Error parsing line: " + str);
				System.exit(-1);
			}
			if (variant.asFeature().overlaps(snp) && variant.asFeature().getStart() == snp.getStart()) {
				if (!options.matchrsid) {
					testSNPObj1 = variant;
				} else if (variant.getId().equals(snp.getName())) {
					testSNPObj1 = variant;
				}
			}
			
			snpStr = inputSNPiter.next();
		}
		reader.close();
		return testSNPObj1;
	}
	
	public double[] convertToDouble(VCFVariant vcfVariant) {
		DoubleMatrix2D alleles = vcfVariant.getGenotypeAllelesAsMatrix2D();
		double[] output = new double[alleles.rows()];
		for (int i = 0; i < alleles.rows(); i++) {
			if (alleles.getQuick(i, 0) == -1) {
				output[i] = -1;
			} else {
				output[i] = (alleles.getQuick(i, 0) + alleles.getQuick(i, 1));
			}
		}
		
		return output;
	}
	
	public void defineCMWindows() throws IOException {
		TextFile tf = new TextFile(options.snpfile, TextFile.R);
		String ln = tf.readLine();
		if (ln == null) {
			System.out.println("Error with snp file. First line is null");
			System.exit(-1);
		}
		int chrCol = 0;
		int posCol = 1;
		int rsCol = 2;
		if (ln.trim().length() > 0) {
			String[] firstlnelems = Strings.whitespace.split(ln);
			if (firstlnelems.length >= 3 &&
					firstlnelems[0].toLowerCase().equals("snp") &&
					firstlnelems[1].toLowerCase().equals("chrom") &&
					firstlnelems[2].toLowerCase().equals("bp")) {
				// goshifter input file format
				System.out.println("File is in GoShifter format...");
				rsCol = 0;
				chrCol = 1;
				posCol = 2;
				ln = tf.readLine();
			} else {
				System.out.println("Assuming three column format");
			}
		}
		ArrayList<SNPFeature> snps = new ArrayList<>();
		// format: chr pos rsid or chr_pos_rsid
		while (ln != null) {
			if (ln.trim().length() > 0) {
				String[] elems = ln.split("\t");
				if (elems.length >= 3) {
					SNPFeature snp = new SNPFeature();
					snp.setStart(Integer.parseInt(elems[posCol]));
					snp.setStop(snp.getStart());
					snp.setName(elems[rsCol]);
					snp.setChromosome(Chromosome.parseChr(elems[chrCol]));
					snps.add(snp);
				} else if (elems.length == 2) {
					// can't quite remember which format this was...
					SNPFeature snp = SNPFeature.parseSNPFeature(elems[1]);
					snps.add(snp);
				}
				ln = tf.readLine();
			}
		}
		tf.close();
		
		// now load the hapmap
		
		
	}
	
	
	public class ProxyFinderTask implements Callable<ArrayList<String>> {
		
		Feature queryVariantFeature;
		
		public ProxyFinderTask(Feature snp) {
			
			this.queryVariantFeature = snp;
			
		}
		
		@Override
		public ArrayList<String> call() throws Exception {
			if (queryVariantFeature == null) {
				return new ArrayList<String>();
			}
			int windowsize = options.windowsize;
			double threshold = options.threshold;
			ArrayList<String> output = new ArrayList<>();
			
			String selfStr = queryVariantFeature.getChromosome().toString()
					+ "\t" + queryVariantFeature.getStart()
					+ "\t" + queryVariantFeature.getName()
					+ "\t" + queryVariantFeature.getChromosome().toString()
					+ "\t" + queryVariantFeature.getStart()
					+ "\t" + queryVariantFeature.getName()
					+ "\t" + 0
					+ "\t" + 1
					+ "\t" + 1;
			
			VCFVariant testSNPObj = getSNP(queryVariantFeature);
			
			
			if (testSNPObj == null) {
				System.out.println("Query: " + queryVariantFeature.getName() + " not in reference. " + queryVariantFeature.getChromosome().getNumber() + ":" + (queryVariantFeature.getStart() - 1) + "-" + (queryVariantFeature.getStart() + 1));
				if (options.addSNPasProxyToItself) {
					output.add(selfStr);
				}
				return output;
			}
			
			
			String tabixfile = options.tabixrefprefix.replace("CHR", "" + testSNPObj.getChrObj().getNumber());
			VCFTabix tabix = new VCFTabix(tabixfile);
			Feature SNPRegion = testSNPObj.asFeature();
			int start = SNPRegion.getStart() - windowsize;
			if (start < 1) {
				start = 1;
			}
			SNPRegion.setStart(start);
			int stop = SNPRegion.getStop();
			if (stop < 0) {
				stop = Integer.MAX_VALUE;
			}
			SNPRegion.setStop(stop);
			
			TabixReader.Iterator window = tabix.query(SNPRegion);
			String next = window.next();
			
			boolean[] snpSampleFilter = null;
			if (options.samplefilter != null) {
				snpSampleFilter = tabix.getSampleFilter(options.samplefilter);
			}
			
			DetermineLD ldcalc = new DetermineLD();
			int variantctr = 0;
			int variantstested = 0;
			int variantsabovethresh = 0;
			while (next != null) {
				// correlate
				VCFVariant variant = new VCFVariant(next, VCFVariant.PARSE.ALL, snpSampleFilter);
				if (variant.equals(testSNPObj) && options.addSNPasProxyToItself) {
					output.add(selfStr);
				} else {
					// correlate
					if (variant.getMAF() > options.mafthreshold && variant.getHwep() > options.hwepthreshold && variant.getAlleles().length == 2) {
						
						Pair<Double, Double> ld = ldcalc.getLD(testSNPObj, variant);
						double rsq = ld.getRight();
						double dpr = ld.getLeft();
						
						if (rsq >= threshold) {
							output.add(queryVariantFeature.getChromosome().toString()
									+ "\t" + queryVariantFeature.getStart()
									+ "\t" + queryVariantFeature.getName()
									+ "\t" + variant.asFeature().getChromosome().toString()
									+ "\t" + variant.asFeature().getStart()
									+ "\t" + variant.asFeature().getName()
									+ "\t" + (Math.abs(queryVariantFeature.getStart() - variant.getPos()))
									+ "\t" + rsq
									+ "\t" + dpr);
							variantsabovethresh++;
						}
						variantstested++;
					}
				}

//				}
				variantctr++;
				next = window.next();
			}
			tabix.close();
			
			if (options.addSNPasProxyToItself && output.isEmpty()) {
				output.add(selfStr);
			}
			
			System.out.println(queryVariantFeature.toString() + " has\t" + variantctr + "\tnearby variants.\t" + variantstested + " were tested.\t" + variantsabovethresh + " actual proxies.");
			return output;
		}
	}
}
