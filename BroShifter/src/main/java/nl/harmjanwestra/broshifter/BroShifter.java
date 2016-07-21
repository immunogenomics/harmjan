package nl.harmjanwestra.broshifter;


import nl.harmjanwestra.broshifter.CLI.BroShifterOptions;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.*;

/**
 * Created by hwestra on 11/11/15.
 */
public class BroShifter {

	boolean DEBUG = false;

	public static void main(String[] args) {

//		/*
//		tring regionsFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";
//			String t1dDir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/Conditional/T1D/";
//			String[] refs = new String[]{"1kg", "seq", "1kg-seq-merged"};
//			String outfilename = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/Posteriors/T1D/";
//		 */
//
//		// TODO: option parsing
//
//
////		String regionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";
////		String outfile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/Posteriors/";
////		String indir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/Posteriors/";
////		String[] listOfAnnotationLists = new String[]{
////				"/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/Peaks/allatacpeaks.txt"
////		};
////		String[] listOfAnnotationListNames = new String[]{
////				"ATAC",
////				"H3K4me3",
////				"DNASE"};
//
//		BroShifter b = new BroShifter();
//		int nrIterations = 10000;
//		String regionFile = "/medpop/srlab/hwestra/tmp/allLoci.bed";
//		String outfile = "/medpop/srlab/hwestra/tmp/Posteriors/";
//		String indir = "/medpop/srlab/hwestra/tmp/Posteriors/";
//		String[] listOfAnnotationLists = new String[]{
//				"/medpop/srlab/hwestra/tmp/annotations/atac.txt",
//				"/medpop/srlab/hwestra/tmp/annotations/pilotatac.txt",
//				"/medpop/srlab/hwestra/tmp/annotations/h3k4me3.txt",
//				"/medpop/srlab/hwestra/tmp/annotations/dnase.txt"
//		};
//		String[] listOfAnnotationListNames = new String[]{
//				"ATAC",
//				"ATACPILOT",
//				"H3K4me3",
//				"DNASE"};
//
//		boolean[] conditional = new boolean[]{
//				false,
//				true,
//				false,
//				false
//		};
//		// String[] refs = new String[]{"1kg", "seq", "1kg-seq-merged"};
//		String[] refs = new String[]{"seq"};
//
//
//		boolean usePeakCenter = true;
//		int bpToExtendAroundCenter = 150;
//		boolean trimRegions = true;
//		int extendRegionEdges = 100;
//		int nrThreads = 64;
//
//		String[] ds = new String[]{"RA", "T1D"};
//
//		for (int d = 0; d < ds.length; d++) {
//			for (int i = 0; i < listOfAnnotationLists.length; i++) {
//				String listOfAnnotations = listOfAnnotationLists[i];
//				String annotationName = listOfAnnotationListNames[i];
//				for (int refId = 0; refId < refs.length; refId++) {
//					String out = outfile + "/" + ds[d] + "-" + refs[refId] + "-" + annotationName;
//					if (conditional[i]) {
//						out += "-conditional";
//					}
//					String posteriorFile = indir + ds[d] + "/";
//					String assocdir = posteriorFile + "/" + refs[refId] + "/";
//					try {
//
//						b.run(regionFile, assocdir, listOfAnnotations, nrIterations, out, conditional[i], usePeakCenter, bpToExtendAroundCenter, trimRegions, extendRegionEdges, nrThreads);
//					} catch (IOException e) {
//						e.printStackTrace();
//					}
//				}
//			}
//		}
	}

	public BroShifter(BroShifterOptions options) throws IOException {
		this.options = options;
		run();
	}

	BroShifterOptions options;

	public void run() throws IOException {

		String regionFile = options.regionFile;
		String directoryWithPosteriors = options.posteriorFile;
		String listOfAnnotations = options.listOfAnnotations;
		int nrIterations = options.nrIterations;
		String outfile = options.outfile;
		boolean conditional = options.conditional;
		boolean usePeakCenter = options.usePeakCenter;
		int bpToExtendAnnotation = options.bpToExtendAnnotation;
		boolean trimRegions = options.trimRegions;
		int defaultRegionExtend = options.defaultRegionExtend;
		int nrThreads = options.nrThreads;

		// TODO: this input path needs some standard formatting or something
		TextFile tf = new TextFile(listOfAnnotations, TextFile.R);
		ArrayList<String> annotationFiles = tf.readAsArrayList();
		System.out.println(annotationFiles.size() + " annotation files in: " + listOfAnnotations);
		tf.close();

		annotationFiles = checkFiles(annotationFiles);
		System.out.println(annotationFiles.size() + " annotation files found on disk.");

		// start threadpool

		ExecutorService threadPool = Executors.newFixedThreadPool(nrThreads);
		CompletionService<Pair<String, ArrayList<String>>> jobHandler = new ExecutorCompletionService<Pair<String, ArrayList<String>>>(threadPool);

		System.out.println("Spinning up a threadpool of " + nrThreads);

		String headerOverall = "P"
				+ "\tP(z)"
				+ "\tz"
				+ "\tSigmaPosterior"
				+ "\tnullMeanSigmaPosterior"
				+ "\tnullSDSigmaPosterior"
				+ "\tGoShifter-P"
				+ "\tGoShifter-P(z)"
				+ "\tGoShifter-z"
				+ "\tGoShifter-OverlappingLoci"
				+ "\tGoShifter-nullMean"
				+ "\tGoShifter-nullSD"
				+ "\tnrSNPsOverlap"
				+ "\tnrSNPs"
				+ "\t%nrSNPsOverlap"
				+ "\tAnnotation"
				+ "\tNrAnnotationRegions";
		if (conditional) {
			headerOverall += "\tAnnotation2";
			headerOverall += "\tNrAnnotation2Regions";
		}

		String headerLocus = "P"
				+ "\tP(z)"
				+ "\tz"
				+ "\tSigmaPosterior"
				+ "\tnullMeanSigmaPosterior"
				+ "\tnullSDSigmaPosterior"
				+ "\tGoShifter-LocusScore"
				+ "\tTestedRegion"
				+ "\tOriginalRegion"
				+ "\tnrSNPsOverlap"
				+ "\tnrSNPs"
				+ "\t%nrSNPsOverlap"
				+ "\tAnnotation"
				+ "\tNrAnnotationRegions";

		if (conditional) {
			headerLocus += "\tAnnotation2";
			headerLocus += "\tNrAnnotation2Regions";
		}


		TextFile outOverall = new TextFile(outfile + "-Overall.txt", TextFile.W);
		TextFile outLocus = new TextFile(outfile + "-Locus.txt", TextFile.W);
		outOverall.writeln(headerOverall);
		outLocus.writeln(headerLocus);


		// for each annotatios
		submitted = 0;
		returned = 0;
		if (conditional) {
			// conditional task
			for (int a1 = 0; a1 < annotationFiles.size(); a1++) {
				String annotation1 = annotationFiles.get(a1);
				for (int a2 = 0; a2 < annotationFiles.size(); a2++) {
					if (a2 != a1 || DEBUG) {

						String annotation2 = annotationFiles.get(a2);
						BroShifterTask task = new BroShifterTask(
								submitted,
								annotation1,
								annotation2,
								options);
						submitted++;
						task.DEBUG = DEBUG;
						jobHandler.submit(task);
					}
				}
			}
		} else {
			for (String annotation1 : annotationFiles) {
				// normal enrichment task
				BroShifterTask task = new BroShifterTask(submitted,
						annotation1,
						null,
						options);
				submitted++;
				task.DEBUG = DEBUG;
				jobHandler.submit(task);
			}
		}

		System.out.println(submitted + " jobs submitted.");

		clearQueue(outOverall, outLocus, jobHandler);


		outOverall.close();
		outLocus.close();

		System.out.println("Main: done testing. Overall results are here: " + outOverall.getFileName());
		System.out.println("Main: done testing. Locus results are here: " + outLocus.getFileName());
		System.out.println();

		threadPool.shutdown();


	}

	int returned = 0;
	int submitted = 0;

	private void clearQueue(TextFile outOverall, TextFile outLocus, CompletionService<Pair<String, ArrayList<String>>> jobHandler) throws IOException {
		System.out.println(submitted + " results to process.");
		while (returned < submitted) {
			try {
				Pair<String, ArrayList<String>> output = jobHandler.take().get();
				if (output != null) {
					String overallStr = output.getLeft();
					outOverall.writeln(overallStr);
					ArrayList<String> locusSpecific = output.getRight();
					for (String s : locusSpecific) {
						outLocus.writeln(s);
					}
					outOverall.flush();
					outLocus.flush();
					if (returned % 10 == 0) {
						System.out.println("\nMain: " + returned + " out of " + submitted + " jobs completed\n");
					}
				}
				returned++;
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}
		returned = 0;
		submitted = 0;
	}


	private ArrayList<String> checkFiles(ArrayList<String> annotationFiles) {
		ArrayList<String> filesPresent = new ArrayList<String>();
		for (int i = 0; i < annotationFiles.size(); i++) {
			if (Gpio.exists(annotationFiles.get(i))) {
				filesPresent.add(annotationFiles.get(i));
			} else {
				System.err.println("WARNING: could not find path: " + annotationFiles.get(i));
			}
		}
		return filesPresent;
	}


}
