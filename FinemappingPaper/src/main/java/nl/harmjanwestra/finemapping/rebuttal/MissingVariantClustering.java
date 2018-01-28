package nl.harmjanwestra.finemapping.rebuttal;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.BedFileFeature;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.console.ProgressBar;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import nl.harmjanwestra.utilities.legacy.genetica.util.RunTimer;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFTabix;
import nl.harmjanwestra.utilities.vcf.VCFVariant;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.concurrent.Executor;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;

public class MissingVariantClustering {
	
	public static void main(String[] args) {
		
		
		if (args.length < 4) {
			System.out.println("Usage: refvcf samplefile out threads");
			System.exit(-1);
		}
		String ref = args[0]; //"d:\\Data\\Ref\\1kg\\ALL.chrCHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";
		String samplefile = args[1];//"d:\\Data\\Ref\\1kg-europeanpopulations.txt.gz";
		String out = args[2]; //"d:\\Data\\Ref\\1kg\\allvars.txt";
		int window = 250000;
		int threads = Integer.parseInt(args[3]);
		MissingVariantClustering v = new MissingVariantClustering();
		try {
			if (!Gpio.exists(samplefile)) {
				samplefile = null;
				System.out.println("Sample file: " + samplefile + " not found. Setting to null");
			}
			v.variantStats(ref, samplefile, window, out, threads);
			
		} catch (IOException e) {
			e.printStackTrace();
		}


//		String str = "1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20";
//		boolean[] include = new boolean[]{
//				true,
//				true,
//				true,
//				true,
//				true,
//				false,
//				true,
//				true,
//				true,
//				true
//		};
//		int offset = 9;
//		String[] output = Strings.subsplit(str, Strings.semicolon, offset, include);
//		for (String s : output) {
//			System.out.println(s);
//		}
//		System.out.println(output.length + " vars");
	}
	
	
	public class KgVariant implements Comparable<KgVariant> {
		SNPFeature f;
		int nproxies;
		double maf;
		double hwep;
		
		@Override
		public int compareTo(KgVariant o) {
			FeatureComparator c = new FeatureComparator();
			return c.compare(this.f, o.f);
		}
		
		@Override
		public boolean equals(Object o) {
			if (this == o) return true;
			if (o == null || getClass() != o.getClass()) return false;
			
			KgVariant kgVariant = (KgVariant) o;
			
			return f != null ? f.equals(kgVariant.f) : kgVariant.f == null;
		}
		
		@Override
		public int hashCode() {
			return f != null ? f.hashCode() : 0;
		}
	}
	
	
	public class KgVariantPair implements Comparable<KgVariantPair> {
		
		KgVariant v1;
		KgVariant v2;
		
		public KgVariantPair(KgVariant v1, KgVariant v2) {
			this.v1 = v1;
			this.v2 = v2;
		}
		
		@Override
		public boolean equals(Object o) {
			if (this == o) return true;
			if (o == null || getClass() != o.getClass()) return false;
			
			KgVariantPair that = (KgVariantPair) o;
			if (this.v1.equals(v1) && this.v2.equals(that.v2) ||
					this.v2.equals(v1) && this.v1.equals(that.v2)) {
				return true;
			} else {
				return false;
			}
		}
		
		@Override
		public int hashCode() {
			int result = v1 != null ? v1.hashCode() : 0;
			result = 31 * result + (v2 != null ? v2.hashCode() : 0);
			return result;
		}
		
		@Override
		public int compareTo(KgVariantPair o) {
			if (this.equals(o)) {
				return 0;
			} else {
				return this.v1.compareTo(o.v1);
			}
		}
	}
	
	public void run(String stat1kgfile, String[] imputedVCFFiles, int ldthreshold, double mafthreshold, String regionsFile) throws IOException {
		
		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(regionsFile);
		
		
		// determine which variants are included
		// TODO: decide whether to include all imputed variants or just the ones with maf>1% and/or made it into the assoc analysis
		ArrayList<ArrayList<VCFVariant>> imputedVariants = new ArrayList<>();
		for (int d = 0; d < imputedVCFFiles.length; d++) {
			ArrayList<VCFVariant> vars = new ArrayList<>();
			TextFile tf = new TextFile(imputedVCFFiles[d], TextFile.R);
			String ln = tf.readLine();
			while (ln != null) {
				if (ln.startsWith("#")) {
					VCFVariant v = new VCFVariant(ln, VCFVariant.PARSE.HEADER);
					String af = v.getInfo().get("AF");
					if (af != null) {
						Double maf = Double.parseDouble(af);
						if (maf > 0.5) {
							maf = 1 - maf;
						}
						if (maf > mafthreshold) {
							KgVariant vs = new KgVariant();
							vs.f = v.asSNPFeature();
							vs.maf = maf;
							String hweps = v.getInfo().get("HWEP");
							if (hweps != null) {
								vs.hwep = Double.parseDouble(hweps);
							}
							vars.add(v);
						}
					}
				}
				ln = tf.readLine();
			}
			tf.close();
			imputedVariants.add(vars);
		}
		
		// determine which variants are in the reference panel
		ArrayList<KgVariant> kgVariantsInRegions = new ArrayList<>();
		ArrayList<KgVariant> kgVariantsNotInRegions = new ArrayList<>();
		TextFile tf = new TextFile(stat1kgfile, TextFile.R);
		tf.readLine();
		String ln = tf.readLine();
		while (ln != null) {
			String[] elems = Strings.tab.split(ln);
			// Variant MAF     HWEP    0       1       2       3       4       5       6       7       8       9
			String variant = new String(elems[0]);
			SNPFeature feature = SNPFeature.parseSNPFeature(variant);
			Double maf = Double.parseDouble(elems[1]);
			if (maf > mafthreshold) {
				Double hwep = Double.parseDouble(elems[2]);
				int nproxies = 0;
				for (int b = 3 + ldthreshold; b < elems.length; b++) {
					nproxies += Integer.parseInt(elems[b]);
				}
				KgVariant v = new KgVariant();
				v.f = feature;
				v.nproxies = nproxies;
				v.maf = maf;
				v.hwep = hwep;
				
				if (feature.overlaps(regions)) {
					kgVariantsInRegions.add(v);
				} else {
					kgVariantsNotInRegions.add(v);
				}
			}
			
			ln = tf.readLine();
		}
		tf.close();
		
		// compare imputed variants with reference variantset
		for (int d = 0; d < imputedVariants.size(); d++) {
			ArrayList<VCFVariant> diseasevariants = imputedVariants.get(d);
			HashSet<String> variantIds = new HashSet<String>();
			for (VCFVariant v : diseasevariants) {
				variantIds.add(v.asFeature().getChromosome().toString() + "_" + v.getPos());
			}
			
			// 2. determine which variants are missing after imputation
			ArrayList<KgVariant> missed1kgVariants = new ArrayList<>();
			for (KgVariant v : kgVariantsInRegions) {
				String id = v.f.getChromosome().toString() + "_" + v.f.getStart();
				if (!variantIds.contains(id)) {
					missed1kgVariants.add(v);
				}
			}
			System.out.println(variantIds.size() + " in imputed set. " + kgVariantsInRegions.size() + " variants in ref. " + missed1kgVariants.size() + " missing.");
			
			// sort variants because sorting is awesome
			Collections.sort(missed1kgVariants);
			
			// 3. measure distance between variants in regions that are missing
			ArrayList<Integer> distances = new ArrayList<>();
			HashSet<KgVariantPair> contributingVariants = new HashSet<>();
			if (missed1kgVariants.size() > 1) {
				for (int v = 0; v < missed1kgVariants.size(); v++) {
					KgVariant var = missed1kgVariants.get(v);
					KgVariant neighbor = null;
					Integer distance = null;
					if (v == 0) {
						// nearest is the next
						neighbor = missed1kgVariants.get(v + 1);
						if (neighbor.f.getChromosome().equals(var.f.getChromosome())) {
							// measure distance
							distance = Math.abs(neighbor.f.getStart() - var.f.getStart());
						}
					} else if (v == missed1kgVariants.size() - 1) {
						// previous one is nearest
						neighbor = missed1kgVariants.get(v - 1);
						if (neighbor.f.getChromosome().equals(var.f.getChromosome())) {
							// measure distance
							distance = Math.abs(neighbor.f.getStart() - var.f.getStart());
						}
					} else {
						
						// either previous or next is nearest
						KgVariant neighbor1 = missed1kgVariants.get(v - 1);
						Integer d1 = null;
						Integer d2 = null;
						if (neighbor1.f.getChromosome().equals(var.f.getChromosome())) {
							// measure distance
							d1 = Math.abs(neighbor1.f.getStart() - var.f.getStart());
						}
						KgVariant neighbor2 = missed1kgVariants.get(v + 1);
						if (neighbor2.f.getChromosome().equals(var.f.getChromosome())) {
							// measure distance
							d2 = Math.abs(neighbor2.f.getStart() - var.f.getStart());
						}
						
						if (d1 != null && d2 != null) {
							distance = Math.min(d1, d2);
							if (distance.equals(d1)) {
								neighbor = neighbor1;
							} else {
								neighbor = neighbor2;
							}
						} else if (d1 != null) {
							distance = d1;
							neighbor = neighbor1;
						} else if (d2 != null) {
							distance = d2;
							neighbor = neighbor2;
						}
					}
					if (distance != null) {
						// got ourselves a missing variant with some friends
						KgVariantPair p = new KgVariantPair(var, neighbor);
						if (contributingVariants.contains(p)) {
							distances.add(distance);
							contributingVariants.add(p); // make sure we're only counting pairs once
						}
					}
				}
			}
			
			// select similar snps from around the genome, measure distances to nearest neighbor?
			
			
		}
	}
	
	public void variantStats(String ref, String samplefile, int window, String out, int threads) throws IOException {
		TextFile tfo = new TextFile(out, TextFile.W);
		String header = "Variant\tMAF\tHWEP";
		for (int i = 0; i < 10; i++) {
			header += "\t" + i;
		}
		tfo.writeln(header);
		
		ExecutorService ex = Executors.newFixedThreadPool(threads);
		
		int rctr = 0;
		
		int totalr = 0;
		
		int chrstart = 1;
		int chrend = 23;
		
		for (int c = chrstart; c < chrend; c++) {
			
			Chromosome chr = Chromosome.parseChr("" + c);
			
			for (int rstart = 0; rstart < chr.getLength(); rstart += window) {
				totalr++;
			}
			
		}
		System.out.println("Total nr of windows: " + totalr);
		Monitor m = new Monitor();
		for (int c = chrstart; c < chrend; c++) {
			String vcf = ref.replaceAll("CHR", "" + c);
			
			Chromosome chr = Chromosome.parseChr("" + c);
			System.out.println((chr.getLength() / window) + " total windows ");
			
			for (int rstart = 0; rstart < chr.getLength(); rstart += window) {
				int overlap = window / 10;
				int rend = rstart + window;
				int overlapStart = rstart - overlap;
				int overlapEnd = rend + overlap;
				
				if (overlapStart < 0) {
					rstart = 0;
				}
				if (overlapEnd > chr.getLength()) {
					overlapEnd = chr.getLength();
				}
				
				if (overlapStart < 0) {
					overlapStart = 0;
				}
				if (overlapEnd > chr.getLength()) {
					overlapEnd = chr.getLength();
				}
				
				Feature testregion = new Feature(chr, rstart, rend);
				Feature getregion = new Feature(chr, overlapStart, overlapEnd);
				
				CalcLDPartnerTask t = new CalcLDPartnerTask(chr, getregion, testregion, samplefile, vcf, tfo, rctr, totalr, m);
				ex.submit(t);
				rctr++;
				
			}
			
		}
		RunTimer timer = new RunTimer();
		while (m.n < totalr) {
			try {
				double deltat = timer.getTimeDiff();
				int iterations = m.n;
				long diff = timer.getTimeDiff() / 1000000000;
				double timePerIter = (double) diff / iterations;
				double timeLeft = timePerIter * (totalr - iterations);
				String strTimeLeft = timer.getTimeDesc(((long) timeLeft) * 1000000000);
				System.out.println(Thread.currentThread().getName() + "\t" + m.n + "/" + totalr + " results returned.. sleeping.\t" + timer.getTimeDesc() + "\tExpecting to be done in: " + strTimeLeft);
				Thread.sleep(10000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		ex.shutdown();
		tfo.close();
		
	}
	
	public class Monitor {
		
		public int n = 0;
		
		public synchronized void iterate() {
			this.n++;
		}
		
	}
	
	public class CalcLDPartnerTask implements Runnable {
		
		private final Feature getregion;
		private final Feature testregion;
		private final int rid;
		private final int totalr;
		private final Monitor monitor;
		Chromosome chr;
		int window;
		String samplefile;
		String vcf;
		TextFile output;
		
		public CalcLDPartnerTask(Chromosome chr, Feature getWindow, Feature testWindow, String samplefile, String vcf, TextFile output, int rid, int totalr, Monitor m) {
			this.chr = chr;
			this.getregion = getWindow;
			this.testregion = testWindow;
			this.samplefile = samplefile;
			this.vcf = vcf;
			this.output = output;
			this.rid = rid;
			this.totalr = totalr;
			this.monitor = m;
		}
		
		@Override
		public void run() {
			try {
				
				System.out.println(Thread.currentThread().getName() + "\tGetting region " + getregion.toString());
				VCFTabix t = new VCFTabix(vcf);
				boolean[] samplestoinclude = null;
				if (samplefile != null) {
					System.out.println(Thread.currentThread().getName() + "\tReading samplefile:  " + samplefile);
					samplestoinclude = t.getSampleFilter(samplefile, vcf);
				}
				
				TabixReader.Iterator it = t.query(getregion);
				String ln = it.next();
				ArrayList<VCFVariant> regionVariants = new ArrayList<VCFVariant>();
				while (ln != null) {
					VCFVariant v = new VCFVariant(ln, VCFVariant.PARSE.ALL, samplestoinclude);
					if (v.getMAF() > 0) {
						regionVariants.add(v);
					}
					ln = it.next();
				}
				t.close();
				
				
				DoubleMatrix2D d2d = DoubleFactory2D.dense.make(regionVariants.size(), regionVariants.size(), -2);
				DetermineLD ldc = new DetermineLD();
				
				int nrCalc = (regionVariants.size() * regionVariants.size()) / 2;
				System.out.println(Thread.currentThread().getName() + "\t" + rid + "/" + totalr + "\tRegion " + getregion.toString() + " has " + regionVariants.size() + " variants. " + nrCalc + " calculations to perform");
				int ctr = 0;
				for (int v = 0; v < regionVariants.size(); v++) {
					VCFVariant var1 = regionVariants.get(v);
					if (var1.asSNPFeature().overlaps(testregion)) {
						for (int v2 = v + 1; v2 < regionVariants.size(); v2++) {
							VCFVariant var2 = regionVariants.get(v2);
							Pair<Double, Double> output = ldc.getLD(var1, var2);
							if (output != null) {
								d2d.setQuick(v, v2, output.getRight());
								d2d.setQuick(v2, v, output.getRight());
							}
							ctr++;
							if (ctr % 1E6 == 0) {
								double perc = ((double) ctr / nrCalc) * 100;
								System.out.println(Thread.currentThread().getName() + "\t" + rid + "/" + totalr + "\t" + getregion.toString() + "\t" + perc + "% done");
							}
						}
					}
				}
				
				System.out.println(Thread.currentThread().getName() + "\t" + rid + "/" + totalr + "\t" + getregion.toString() + "\t" + 100 + "% done");
				int[] nrPartners = new int[10];
				for (int v = 0; v < regionVariants.size(); v++) {
					VCFVariant var1 = regionVariants.get(v);
					
					if (var1.asSNPFeature().overlaps(testregion)) {
						for (int q = 0; q < nrPartners.length; q++) {
							nrPartners[q] = 0;
						}
						for (int v2 = 0; v2 < regionVariants.size(); v2++) {
							if (v != v2) {
								double ld = d2d.getQuick(v, v2);
								if (ld != -2) {
									int bin = (int) Math.floor(ld * nrPartners.length);
									if (bin == nrPartners.length) {
										bin = nrPartners.length - 1;
									}
									nrPartners[bin]++;
								}
							}
						}
						
						output.writelnsynced(var1.asSNPFeature().toString() + "\t" + var1.getMAF() + "\t" + var1.getHwep() + "\t" + Strings.concat(nrPartners, Strings.tab));
					}
				}
				monitor.iterate();
				
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
	}
	
}
