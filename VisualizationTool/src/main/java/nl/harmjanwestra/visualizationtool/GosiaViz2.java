/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.visualizationtool;

import com.lowagie.text.DocumentException;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.*;
import nl.harmjanwestra.utilities.gtf.GTFAnnotation;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.NavigableSet;

/**
 * @author hwestra
 */
public class GosiaViz2 {

	private Track[][] bedFileTracks;
	private Track[] peakTracks;
	private boolean mergeMarks = true;

	public void run(String locusScoreFile, String snpOverlapFile,
					String gtfFile, String bedfileDir, String peakFile, String snpAnnotationFile, String outputFileName, int extraWindowSize) throws IOException, DocumentException {

//        if (enrichFile == null || !Gpio.exists(enrichFile)) {
//            throw new IllegalArgumentException(enrichFile + " does not exist or is null");
//        }
		if (locusScoreFile == null || !Gpio.exists(locusScoreFile)) {
			throw new IllegalArgumentException(locusScoreFile + " does not exist or is null");
		}
		if (snpOverlapFile == null || !Gpio.exists(snpOverlapFile)) {
			throw new IllegalArgumentException(snpOverlapFile + " does not exist or is null");
		}
		if (snpAnnotationFile == null || !Gpio.exists(snpAnnotationFile)) {
			throw new IllegalArgumentException(snpAnnotationFile + " does not exist or is null");
		}
		if (gtfFile == null || !Gpio.exists(gtfFile)) {
			throw new IllegalArgumentException(gtfFile + " does not exist or is null");
		}
		if (bedfileDir == null || !Gpio.exists(bedfileDir)) {
			throw new IllegalArgumentException(bedfileDir + " does not exist or is null");
		}
		if (outputFileName == null) {
			throw new IllegalArgumentException("outputfilename == null");
		}
		String outputFileNameLC = outputFileName.toLowerCase();
		if (!outputFileNameLC.endsWith("jpg") || !outputFileNameLC.endsWith("pdf")) {
			outputFileName += ".pdf";
		}

        /*
		 // load enrichfile
         nperm   nSnpOverlap     allSnps enrichment
         0       29      68      0.42647
         1       22      68      0.32353
         2       22      68      0.32353
         3       19      68      0.27941
         */
		// load the locusscores
		loadLocusScores(locusScoreFile);

		// load the LD variants for each locus
		loadLDSNPs(snpOverlapFile);

		// load SNP annotation for loaded SNPs
		loadSNPAnnotationAndSetLocusBounds(snpAnnotationFile, extraWindowSize);

		annotateLociWithGenes(gtfFile);
////        System.exit(-1);
//		// load the chromatin peaks
//		loadChromatinPeaks(peakFile);
//
//		// load the tracks of annotations
//		loadChromatinMarks(bedfileDir);

		// plot loci, genes, annotations etc..
		Plot plot = new Plot();
		plot.plot(loci, bedFileTracks, peakTracks, outputFileName);

	}

	ArrayList<Locus> loci = new ArrayList<Locus>();

	private void loadLocusScores(String locusScoreFile) throws IOException {
		TextFile tf = new TextFile(locusScoreFile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);

		while (elems != null) {

			String snp = elems[0];
			String overlap = elems[1];
			String locusscore = elems[2];
			Locus l = new Locus();
			if (overlap.equals("1")) {
				l.setOverlap(true);
			}
			l.setScore(Double.parseDouble(locusscore));
			l.setName(snp);
			loci.add(l);
			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();

		System.out.println(loci.size() + " loci loaded from: +" + locusScoreFile);
	}

	private void loadLDSNPs(String snpOverlapFile) throws IOException {
		HashMap<String, Locus> strToLocus = new HashMap<String, Locus>();
		for (Locus l : loci) {
			strToLocus.put(l.getName(), l);
		}

		TextFile tf = new TextFile(snpOverlapFile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);

		while (elems != null) {

			String snp = elems[0];
			String snp2 = elems[1];
			String overlap = elems[2];

			boolean b = false;
			Locus l = strToLocus.get(snp);
			if (l == null) {
				//System.err.println("Locus not loaded: " + snp);
			} else {
				if (overlap.equals("1")) {
					b = true;
				}
				if (!snp2.equals(snp)) {
					System.out.println("Adding proxy: " + snp2 + " for locus: " + l.getName() + " - " + l.getProxies().size());
					l.addProxy(snp2, b);
				}
			}
			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();
		System.out.println("Done loading proxies: " + snpOverlapFile);
	}

	private void loadSNPAnnotationAndSetLocusBounds(String snpAnnotationFile, int extraWindowSize) throws IOException {

		HashMap<String, Feature> strToFeature = new HashMap<String, Feature>();

		for (Locus l : loci) {
			if (strToFeature.containsKey(l.getName())) {
				System.err.println("ERROR: duplicate SNP in locus definitions: " + l.toString());
			} else {
				strToFeature.put(l.getName(), l);
				for (Feature f : l.proxies) {
					if (strToFeature.containsKey(f.getName())) {
						System.err.println("ERROR: duplicate SNP in locus proxy definitions: " + f.toString());
					} else {
						strToFeature.put(f.getName(), f);
					}
				}
			}
		}

		System.out.println("Loading SNP annotation: " + snpAnnotationFile);
		TextFile tf = new TextFile(snpAnnotationFile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);

		while (elems != null) {

			String snp = elems[0];

			Feature f = strToFeature.get(snp);
			if (f != null) {
				Chromosome chr = Chromosome.parseChr(elems[1]);
				Integer position = Integer.parseInt(elems[2]);

				f.setChromosome(chr);
				f.setStart(position);
				f.setStop(position);
			}

			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();

		for (Locus l : loci) {

			int maxStart = l.getStart();
			int maxStop = l.getStop();
			for (Feature f : l.getProxies()) {
				if (f.getStart() < maxStart) {
					maxStart = f.getStart();
				}
				if (f.getStop() > maxStop) {
					maxStop = f.getStop();
				}
			}
			maxStart -= extraWindowSize;
			maxStop += extraWindowSize;

			l.setLeftBound(maxStart);
			l.setRightBound(maxStop);

			System.out.println(l.getName() + "\t" + l.getChromosome() + ":" + l.getLeftBound() + "-" + l.getRightBound() + "\t" + (l.getRightBound() - l.getLeftBound()));
		}

	}

//	private void loadChromatinMarks(String bedfileDir) throws IOException {
//
//		bedfileDir = Gpio.formatAsDirectory(bedfileDir);
//		String[] files = Gpio.getListOfFiles(bedfileDir);
//		int nrBedFiles = 0;
//		for (String bedfile : files) {
//			if (bedfile.toLowerCase().endsWith(".bed") || bedfile.toLowerCase().endsWith(".bed.gz")) {
//				nrBedFiles++;
//			}
//		}
//
//		if (mergeMarks) {
//			bedFileTracks = new Track[loci.size()][1];
//			int bedFileNr = 0;
//			BedFileReader bfr = new BedFileReader();
//
//			for (String bedfile : files) {
//				bedfile = bedfileDir + bedfile;
//
//				if (bedfile.toLowerCase().endsWith(".bed") || bedfile.toLowerCase().endsWith(".bed.gz")) {
//					System.out.println("Loading: " + bedfile);
//					Track t = bfr.read(bedfile, bedfile, null, 0, 0, false);
//					System.out.println(loci.size() + " loci will be annotated.");
//					int locusId = 0;
//					for (Locus l : loci) {
//						int maxStart = l.getLeftBound();
//						int maxStop = l.getRightBound();
//
//						if (bedFileTracks[locusId][0] == null) {
//							bedFileTracks[locusId][0] = t.getSubset(l.getChromosome(), maxStart, maxStop);
//						} else {
//							bedFileTracks[locusId][0].addFeatures(t.getSubset(l.getChromosome(), maxStart, maxStop));
//						}
//
//						locusId++;
//					}
//
//				}
//
//			}
//
//		} else {
//			bedFileTracks = new Track[loci.size()][nrBedFiles];
//			int bedFileNr = 0;
//			BedFileReader bfr = new BedFileReader();
//			for (String bedfile : files) {
//				bedfile = bedfileDir + bedfile;
//
//				if (bedfile.toLowerCase().endsWith(".bed") || bedfile.toLowerCase().endsWith(".bed.gz")) {
//
//					int locusId = 0;
//					System.out.println("Loading: " + bedfile);
//					Track t = bfr.read(bedfile, bedfile, null, 0, 0, false);
//					System.out.println(loci.size() + " loci will be annotated.");
//					for (Locus l : loci) {
//						int maxStart = l.getLeftBound();
//						int maxStop = l.getRightBound();
//						bedFileTracks[locusId][bedFileNr] = t.getSubset(l.getChromosome(), maxStart, maxStop);
//						System.out.println(locusId + "\t" + bedFileTracks[locusId][bedFileNr].getNrFeatures());
//						locusId++;
//					}
//				}
//				bedFileNr++;
//			}
//		}
//
//	}

//	private void loadChromatinPeaks(String peakfile) throws IOException {
//
//		System.out.println("Loading peaks: " + peakfile);
//		peakTracks = new Track[loci.size()];
//
//		BedFileReader bfr = new BedFileReader();
//		int locusId = 0;
//		Track t = bfr.read(peakfile, peakfile, null, 0, 0, false);
//
//		for (Locus l : loci) {
//			int maxStart = l.getLeftBound();
//			int maxStop = l.getRightBound();
//			peakTracks[locusId] = t.getSubset(l.getChromosome(), maxStart, maxStop); //bfr.read(peakfile, peakfile, l.getChromosome(), maxStart, maxStop, false);
//			locusId++;
//		}
//	}
	private void annotateLociWithGenes(String gtfFile) throws IOException {
		GTFAnnotation ann = new GTFAnnotation(gtfFile);

		Gene ellgene = null;
		Collection<Gene> genes = ann.getGenes();
		for (Gene g : genes) {
			if (g.getGeneId().equals("DIS3L2")) {
				System.out.println("Found it! " + g.toString());
				ellgene = g;
			}
		}

		for (Locus l : loci) {
			int locusStart = l.getLeftBound();
			int locusStop = l.getRightBound();
			NavigableSet<Gene> s = ann.getGeneTree().subSet(new Gene("", l.getChromosome(), Strand.POS, locusStart, locusStart), true,
					new Gene("", l.getChromosome(), Strand.NEG, locusStop, locusStop), true);

			System.out.println(s.size() + " genes found overlapping with locus: " + l.getName());
			for (Gene g : s) {
				System.out.println("Adding gene: " + g.toString() + " to locus: " + l.getName());
				l.addGene(g);
			}

			if (ellgene != null) {

				NavigableSet<Gene> s2 = ann.getGeneTree().subSet(
						new Gene("", ellgene.getChromosome(), ellgene.getStrand(), ellgene.getStart() + 1, ellgene.getStart() + 1),
						true,
						new Gene("", ellgene.getChromosome(), ellgene.getStrand(), ellgene.getStop() - 1, ellgene.getStop() - 1),
						true
				);

				System.out.println(s2.size() + "genes match query gene");

				if (l.overlaps(ellgene)) {
					System.out.println("Overlap!");
				} else {
					System.out.println("No overlap");
				}

				System.out.println(locusStart + "-" + locusStop);
				System.out.println(ellgene.getStart() + "-" + ellgene.getStop());
			}
		}

	}

}
