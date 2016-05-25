package nl.harmjanwestra.broshifter;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.broshifter.CLI.BroShifterOptions;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.bedfile.BedGraphReader;
import nl.harmjanwestra.utilities.features.*;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.*;
import nl.harmjanwestra.utilities.gtf.GTFAnnotation;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

/**
 * Created by hwestra on 2/24/16.
 */
public class AnnotationOverlapPlot {


	private final BroShifterOptions options;

	public AnnotationOverlapPlot(BroShifterOptions options) throws IOException, DocumentException {
		this.options = options;
		if (options.overlapmatrix) {
			this.overlapMatrix();
		} else if (options.continuous) {
			this.plotContinuousTrait();
		} else {
			this.plotBinaryTrait();
		}

	}

	public static void main(String[] args) {
		// debug

		args = new String[11];
		args[0] = "--plotoverlap";
		args[1] = "--gtf";
		args[2] = "/Data//Ref/Annotation/UCSC/genes.gtf";
		args[3] = "-a";
		args[4] = "/Data/tmp/2016-03-08/meh.txt";

		args[5] = "-p";
		args[6] = "/Data/tmp/2016-03-08/swisscheesebeagle41-posterior.txt";

		args[7] = "-r";
		args[8] = "/Data/tmp/2016-03-08/cd28.bed";

		args[9] = "-o";
		args[10] = "/Data/tmp/2016-03-08/tmp.pdf";


		BroShifterOptions options = new BroShifterOptions(args);

		try {
			AnnotationOverlapPlot plot = new AnnotationOverlapPlot(options);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}


	}

	private Track[] loadAnnotations(ArrayList<Feature> regions) throws IOException {
		TextFile tf = new TextFile(options.listOfAnnotations, TextFile.R);

		String[] elems = tf.readLineElems(TextFile.tab);
		ArrayList<String> annotationFiles = new ArrayList<String>();
		HashMap<String, String> fileToName = null;
		while (elems != null) {
			if (elems.length == 1) {
				annotationFiles.add(elems[0]);
			} else {
				if (fileToName == null) {
					fileToName = new HashMap<String, String>();
				}
				fileToName.put(elems[1], elems[0]);
				annotationFiles.add(elems[1]);
			}
			elems = tf.readLineElems(TextFile.tab);
		}


		Collections.sort(annotationFiles);
		System.out.println(annotationFiles.size() + " annotation files in: " + options.listOfAnnotations);
		tf.close();


		// load annotations
		Track[] allAnnotations = new Track[annotationFiles.size()];
		AnnotationLoader loader = new AnnotationLoader();
		for (int i = 0; i < annotationFiles.size(); i++) {
			allAnnotations[i] = loader.loadAnnotations(annotationFiles.get(i), options.usePeakCenter, options.bpToExtendAnnotation, true, regions);

			if (fileToName != null) {
				String name = fileToName.get(annotationFiles.get(i));
				allAnnotations[i].setName(name);
			}

		}
		return allAnnotations;
	}

	public void plotBinaryTrait() throws IOException, DocumentException {

		if (options.geneAnnotationFile == null) {
			System.out.println("Error: provide --gtf");
			System.exit(-1);
		}
		GTFAnnotation geneannotation = new GTFAnnotation(options.geneAnnotationFile);


		// load bed regions to testNormal
		BedFileReader bf = new BedFileReader();
		ArrayList<Feature> regions = bf.readAsList(options.regionFile);

		Track[] allAnnotations = loadAnnotations(regions);

		// load posteriors
		BroShifterTask bs = new BroShifterTask(options);


		System.out.println("All annotations loaded. Now iterating " + regions.size() + " loci");
		// iterate regions
		for (int r = 0; r < regions.size(); r++) {

			Feature region = regions.get(r);

			AssociationFile assocFile = new AssociationFile();
			ArrayList<AssociationResult> results = assocFile.read(options.posteriorFile, region);
			ApproximateBayesPosterior ab = new ApproximateBayesPosterior();
			ArrayList<AssociationResult> credibleSet = ab.createCredibleSet(results, options.credibleSetThreshold);
			ArrayList<SNPFeature> snps = new ArrayList<>(results.size());
			ArrayList<SNPFeature> credibleSetSNPFeatures = new ArrayList<>(credibleSet.size());
			HashSet<SNPFeature> credibleSetHash = new HashSet<SNPFeature>();
			for (AssociationResult result : credibleSet) {
				Feature f = result.getSnp();
				SNPFeature f2 = new SNPFeature();
				f2.setChromosome(f.getChromosome());
				f2.setStart(f.getStart());
				f2.setStop(f.getStop());
				f2.setName(f.getName());
				f2.setP(result.getPosterior());
				credibleSetSNPFeatures.add(f2);
				credibleSetHash.add(f2);
			}


			ArrayList<Pair<Integer, Double>> allPvalues = new ArrayList<>();
			boolean[] mark = new boolean[results.size()];
			int ctr = 0;
			for (AssociationResult result : results) {
				Feature f = result.getSnp();

				SNPFeature f2 = new SNPFeature();
				f2.setChromosome(f.getChromosome());
				f2.setStart(f.getStart());
				f2.setStop(f.getStop());
				f2.setName(f.getName());

				f2.setP(result.getPosterior());
				allPvalues.add(new Pair<Integer, Double>(f.getStart(), result.getLog10Pval()));
				snps.add(f2);
				if (credibleSetHash.contains(f2)) {
					mark[ctr] = true;
				}
				ctr++;

			}

			ArrayList<Triple<Track, Double, Integer>> annotationsUnsorted = new ArrayList<>();

			System.out.println(snps.size() + " variants loaded");

			for (int i = 0; i < allAnnotations.length; i++) {
				Track annot = allAnnotations[i].getSubset(region.getChromosome(), region.getStart(), region.getStop());
				Pair<Double, Integer> overlap = bs.getOverlap(annot, credibleSetSNPFeatures);
				if (overlap.getRight() > 0) {
					Triple<Track, Double, Integer> t = new Triple<>(annot, overlap.getLeft(), overlap.getRight());
					annotationsUnsorted.add(t);
				}
			}


			System.out.println(annotationsUnsorted.size() + " annotations overlap with credible sets");
			if (annotationsUnsorted.isEmpty()) {
				System.out.println("No annotations to plotBinaryTrait");
				System.exit(0);
			}

			// sort the annotations by overlap
			Collections.sort(annotationsUnsorted, new AssocTripleComparator());

			ArrayList<Track> annotationsSorted = new ArrayList<>();
			DecimalFormat format = new DecimalFormat("#.###");
			for (int q = 0; q < annotationsUnsorted.size(); q++) {
				Triple<Track, Double, Integer> t = annotationsUnsorted.get(q);
				t.getLeft().setName(format.format(t.getMiddle()) + " - " + t.getRight() + " - " + t.getLeft().getName());
				annotationsSorted.add(t.getLeft());
			}

			// determine size of plotBinaryTrait
			AnnotationPanel annotPanel = new AnnotationPanel(annotationsSorted.size(), 1);
			annotPanel.setData(region, annotationsSorted);
			annotPanel.setOverlappingFeatures(credibleSetSNPFeatures);

			int pixelspertrack = annotPanel.getTrackheight() + annotPanel.getMarginBetween();
			int genepanelrows = (int) Math.ceil(100 / pixelspertrack);
			int assocPanelRows = (int) Math.ceil(300 / pixelspertrack);

			System.out.println(1000 + "x" + pixelspertrack + " grid");
			System.out.println(annotationsSorted.size() + genepanelrows + assocPanelRows + " rows ");

			Grid grid = new Grid(1000, pixelspertrack, 10 + annotationsSorted.size() + genepanelrows + assocPanelRows, 1, 100, 0);


			TreeSet<Gene> genes = geneannotation.getGeneTree();
			Gene geneStart = new Gene("", region.getChromosome(), Strand.POS, region.getStart(), region.getStart());
			Gene geneStop = new Gene("", region.getChromosome(), Strand.POS, region.getStop(), region.getStop());
			SortedSet<Gene> overlappingGenes = genes.subSet(geneStart, true, geneStop, true);

			ArrayList<Gene> overlappingGenesList = new ArrayList<>();
			overlappingGenesList.addAll(overlappingGenes);

			GenePanel genePanel = new GenePanel(genepanelrows, 1);
			genePanel.setData(region, overlappingGenesList);

			AssociationPanel assocPanel = new AssociationPanel(assocPanelRows, 1);
			assocPanel.setMarkDifferentColor(mark);

			assocPanel.setDataSingleDs(region, null, allPvalues, region.toString());

			grid.addPanel(genePanel, 0, 0);
			grid.addPanel(new SpacerPanel(5, 1), genepanelrows, 0);
			grid.addPanel(assocPanel, genepanelrows + 5, 0);
			grid.addPanel(new SpacerPanel(5, 1), genepanelrows + 5 + assocPanelRows, 0);
			grid.addPanel(annotPanel, genepanelrows + 5 + assocPanelRows + 5, 0);

			grid.draw(options.outfile + region.toString() + ".pdf");

		}

	}

	// TODO: ugh duplicate code...
	public void plotContinuousTrait() throws IOException, DocumentException {

		if (options.geneAnnotationFile == null) {
			System.out.println("Error: provide --gtf");
			System.exit(-1);
		}
		GTFAnnotation geneannotation = new GTFAnnotation(options.geneAnnotationFile);


		// load bed regions to testNormal
		BedFileReader bf = new BedFileReader();
		ArrayList<Feature> regions = bf.readAsList(options.regionFile);


		TextFile tf = new TextFile(options.listOfAnnotations, TextFile.R);

		String[] elems = tf.readLineElems(TextFile.tab);
		ArrayList<String> annotationFiles = new ArrayList<String>();
		ArrayList<String> sampleNames = null;
		while (elems != null) {
			if (elems.length == 1) {
				annotationFiles.add(elems[0]);
			} else {
				if (sampleNames == null) {
					sampleNames = new ArrayList<>();
				}
				sampleNames.add(elems[0]);

				annotationFiles.add(elems[1]);
			}
			elems = tf.readLineElems(TextFile.tab);
		}


		// load posteriors
		BroShifterTask bs = new BroShifterTask(options);


		// load bedgraphs
		BedGraphReader bedGraphReader = new BedGraphReader();
		ArrayList<ArrayList<BedGraphFeature>> bedGraphData = new ArrayList<>();

		for (int i = 0; i < annotationFiles.size(); i++) {
			bedGraphData.add(bedGraphReader.read(annotationFiles.get(i), true, null));
		}

		System.out.println("All annotations loaded. Now iterating " + regions.size() + " loci");
		// iterate regions
		for (int r = 0; r < regions.size(); r++) {

			Feature region = regions.get(r);
			AssociationFile assocFile = new AssociationFile();
			ArrayList<AssociationResult> results = assocFile.read(options.posteriorFile, region);
			ApproximateBayesPosterior ab = new ApproximateBayesPosterior();
			ArrayList<AssociationResult> credibleSet = ab.createCredibleSet(results, options.credibleSetThreshold);
			ArrayList<SNPFeature> snps = new ArrayList<>(results.size());
			ArrayList<SNPFeature> credibleSetSNPFeatures = new ArrayList<>(credibleSet.size());
			HashSet<SNPFeature> credibleSetHash = new HashSet<SNPFeature>();
			for (AssociationResult result : credibleSet) {
				Feature f = result.getSnp();
				SNPFeature f2 = new SNPFeature();
				f2.setChromosome(f.getChromosome());
				f2.setStart(f.getStart());
				f2.setStop(f.getStop());
				f2.setName(f.getName());
				f2.setP(result.getPosterior());
				credibleSetSNPFeatures.add(f2);
				credibleSetHash.add(f2);
			}


			ArrayList<Pair<Integer, Double>> allPvalues = new ArrayList<>();
			boolean[] mark = new boolean[results.size()];
			int ctr = 0;
			for (AssociationResult result : results) {
				Feature f = result.getSnp();

				SNPFeature f2 = new SNPFeature();
				f2.setChromosome(f.getChromosome());
				f2.setStart(f.getStart());
				f2.setStop(f.getStop());
				f2.setName(f.getName());

				f2.setP(result.getPosterior());
				allPvalues.add(new Pair<Integer, Double>(f.getStart(), result.getLog10Pval()));
				snps.add(f2);
				if (credibleSetHash.contains(f2)) {
					mark[ctr] = true;
				}
				ctr++;

			}

			GraphAnnotationPanel graphAnnotationPanel = new GraphAnnotationPanel(annotationFiles.size(), 1);
			ArrayList<ArrayList<BedGraphFeature>> dataForRegion = filterBedGraphData(bedGraphData, region);


			if (sampleNames == null) {
				graphAnnotationPanel.setData(region, dataForRegion, sampleNames.toArray(new String[0]));
			} else {
				graphAnnotationPanel.setData(region, dataForRegion, annotationFiles.toArray(new String[0]));
			}

			// determine size of plotBinaryTrait
			int pixelspertrack = graphAnnotationPanel.getTrackheight() + graphAnnotationPanel.getMarginBetween();
			int genepanelrows = (int) Math.ceil(100 / pixelspertrack);
			int assocPanelRows = (int) Math.ceil(300 / pixelspertrack);
//
			System.out.println(1000 + "x" + pixelspertrack + " grid");
			System.out.println(annotationFiles.size() + genepanelrows + assocPanelRows + " rows ");
//
			Grid grid = new Grid(1000, pixelspertrack, 10 + annotationFiles.size() + genepanelrows + assocPanelRows, 1, 100, 0);


			TreeSet<Gene> genes = geneannotation.getGeneTree();
			Gene geneStart = new Gene("", region.getChromosome(), Strand.POS, region.getStart(), region.getStart());
			Gene geneStop = new Gene("", region.getChromosome(), Strand.POS, region.getStop(), region.getStop());
			SortedSet<Gene> overlappingGenes = genes.subSet(geneStart, true, geneStop, true);

			ArrayList<Gene> overlappingGenesList = new ArrayList<>();
			overlappingGenesList.addAll(overlappingGenes);

			GenePanel genePanel = new GenePanel(genepanelrows, 1);
			genePanel.setData(region, overlappingGenesList);

			AssociationPanel assocPanel = new AssociationPanel(assocPanelRows, 1);
			assocPanel.setMarkDifferentColor(mark);

			assocPanel.setDataSingleDs(region, null, allPvalues, region.toString());

			grid.addPanel(genePanel, 0, 0);
			grid.addPanel(new SpacerPanel(5, 1), genepanelrows, 0);
			grid.addPanel(assocPanel, genepanelrows + 5, 0);
			grid.addPanel(new SpacerPanel(5, 1), genepanelrows + 5 + assocPanelRows, 0);
			grid.addPanel(graphAnnotationPanel, genepanelrows + 5 + assocPanelRows + 5, 0);

			grid.draw(options.outfile + region.toString() + ".pdf");

		}

	}

	private ArrayList<ArrayList<BedGraphFeature>> filterBedGraphData(ArrayList<ArrayList<BedGraphFeature>> bedGraphData, Feature region) {

		ArrayList<ArrayList<BedGraphFeature>> output = new ArrayList<>();
		for (int i = 0; i < bedGraphData.size(); i++) {
			ArrayList<BedGraphFeature> data = bedGraphData.get(i);
			ArrayList<BedGraphFeature> dsOutput = new ArrayList<>();
			for (int q = 0; q < data.size(); q++) {
				BedGraphFeature f = data.get(q);
				if (region.overlaps(f)) {
					dsOutput.add(f);
				}
			}
			output.add(dsOutput);
		}
		return output;
	}

	public void overlapMatrix() throws IOException {
		if (options.geneAnnotationFile == null) {
			System.out.println("Error: provide --gtf");
			System.exit(-1);
		}

		// load bed regions to testNormal
		BedFileReader bf = new BedFileReader();
		ArrayList<Feature> regions = bf.readAsList(options.regionFile);

		Track[] allAnnotations = loadAnnotations(regions);

		// load posteriors
		BroShifterTask bs = new BroShifterTask(options);

		double[][] overlapMatrix = new double[regions.size()][allAnnotations.length];
		int[][] overlapMatrixCount = new int[regions.size()][allAnnotations.length];
		int[] credibleSetSizePerRegion = new int[regions.size()];

		for (int r = 0; r < regions.size(); r++) {

			Feature region = regions.get(r);

			AssociationFile assocFile = new AssociationFile();
			ArrayList<AssociationResult> results = assocFile.read(options.posteriorFile, region);
			ApproximateBayesPosterior ab = new ApproximateBayesPosterior();
			ArrayList<AssociationResult> credibleSet = ab.createCredibleSet(results, options.credibleSetThreshold);
			ArrayList<SNPFeature> snps = new ArrayList<>(results.size());
			ArrayList<SNPFeature> credibleSetSNPFeatures = new ArrayList<>(credibleSet.size());
			HashSet<SNPFeature> credibleSetHash = new HashSet<SNPFeature>();
			for (AssociationResult result : credibleSet) {
				Feature f = result.getSnp();
				SNPFeature f2 = new SNPFeature();
				f2.setChromosome(f.getChromosome());
				f2.setStart(f.getStart());
				f2.setStop(f.getStop());
				f2.setName(f.getName());
				f2.setP(result.getPosterior());
				credibleSetSNPFeatures.add(f2);
				credibleSetHash.add(f2);
			}

			credibleSetSizePerRegion[r] = credibleSet.size();

			System.out.println(snps.size() + " variants loaded");

			for (int i = 0; i < allAnnotations.length; i++) {
				Track annot = allAnnotations[i].getSubset(region.getChromosome(), region.getStart(), region.getStop());
				Pair<Double, Integer> overlap = bs.getOverlap(annot, credibleSetSNPFeatures);
				if (overlap.getRight() > 0) {
					overlapMatrix[r][i] = overlap.getLeft();
					overlapMatrixCount[r][i] = overlap.getRight();
				}
			}
		}

		// remove annotations with zero overlap
		Pair<double[][], ArrayList<String>> pruned = prune(overlapMatrix, allAnnotations, regions);
		overlapMatrix = pruned.getLeft();

		String header = "-";
		TextFile out = new TextFile(options.outfile, TextFile.W);
		for (int r = 0; r < regions.size(); r++) {
			header += "\t" + regions.get(r).toString();
		}
		out.writeln(header);
		for (int i = 0; i < pruned.getRight().size(); i++) {
			String ln = pruned.getRight().get(i);
			for (int r = 0; r < regions.size(); r++) {
				ln += "\t" + overlapMatrix[r][i];
			}
			out.writeln(ln);
		}
		out.close();


	}

	private Pair<double[][], ArrayList<String>> prune(double[][] overlapMatrix, Track[]
			allAnnotations, ArrayList<Feature> regions) {
		// prune the matrix
		boolean[] includeAnnotation = new boolean[allAnnotations.length];
		int nrToInclude = 0;
		ArrayList<String> colNames = new ArrayList<>();
		for (int i = 0; i < allAnnotations.length; i++) {
			double sum = 0;
			for (int r = 0; r < regions.size(); r++) {
				sum += overlapMatrix[r][i];
			}
			if (sum > 0d) {
				includeAnnotation[i] = true;
				colNames.add(allAnnotations[i].getName());
				nrToInclude++;
			}
		}

		double[][] output = new double[overlapMatrix.length][nrToInclude];

		for (int r = 0; r < regions.size(); r++) {
			int ctr = 0;
			for (int i = 0; i < allAnnotations.length; i++) {
				if (includeAnnotation[i]) {
					output[r][ctr] = overlapMatrix[r][i];
					ctr++;
				}
			}
		}
		return new Pair<>(output, colNames);
	}

	private class AssocTripleComparator implements Comparator<Triple<Track, Double, Integer>> {

		@Override
		public int compare(Triple<Track, Double, Integer> o1, Triple<Track, Double, Integer> o2) {
			if (o1.getMiddle() > o2.getMiddle()) {
				return -1;
			} else if (o1.getMiddle() < o2.getMiddle()) {
				return 1;
			} else {
				if (o1.getRight() > o2.getRight()) {
					return -1;
				} else if (o1.getRight() < o2.getRight()) {
					return 1;
				} else {
					return 0;
				}
			}
		}
	}

}
