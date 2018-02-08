package nl.harmjanwestra.finemapping.rebuttal;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
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
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class CalculateLDVariantStats {
	
	
	public static void mainOld(String[] args) {
		
		
		if (args.length < 4) {
			System.out.println("Usage: refvcf samplefile out threads");
			System.exit(-1);
		}
		String ref = args[0]; //"d:\\Data\\Ref\\1kg\\ALL.chrCHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";
		String samplefile = args[1];//"d:\\Data\\Ref\\1kg-europeanpopulations.txt.gz";
		String out = args[2]; //"d:\\Data\\Ref\\1kg\\allvars.txt";
		int window = 250000;
		int threads = Integer.parseInt(args[3]);
		CalculateLDVariantStats v = new CalculateLDVariantStats();
		try {
			if (!Gpio.exists(samplefile)) {
				samplefile = null;
				System.out.println("Sample file: " + samplefile + " not found. Setting to null");
			}
			v.calculateVariantStats(ref, samplefile, window, out, threads);
			
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
	
	
	public void calculateVariantStats(String ref, String samplefile, int window, String out, int threads) throws IOException {
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
