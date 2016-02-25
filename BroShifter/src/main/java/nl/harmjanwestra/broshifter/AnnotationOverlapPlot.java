package nl.harmjanwestra.broshifter;

import com.itextpdf.text.DocumentException;
import com.sun.deploy.association.Association;
import nl.harmjanwestra.broshifter.CLI.BroShifterOptions;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.*;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.AnnotationPanel;
import nl.harmjanwestra.utilities.graphics.panels.AssociationPanel;
import nl.harmjanwestra.utilities.graphics.panels.GenePanel;
import nl.harmjanwestra.utilities.gtf.GTFAnnotation;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;

import java.io.File;
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
		this.plot();
	}

	public void plot() throws IOException, DocumentException {


		TextFile tf = new TextFile(options.listOfAnnotations, TextFile.R);
		ArrayList<String> annotationFiles = tf.readAsArrayList();
		System.out.println(annotationFiles.size() + " annotation files in: " + options.listOfAnnotations);
		tf.close();

		// load bed regions to test
		BedFileReader bf = new BedFileReader();
		ArrayList<Feature> regions = bf.readAsList(options.regionFile);

		// load annotations
		Track[] allAnnotations = new Track[annotationFiles.size()];
		AnnotationLoader loader = new AnnotationLoader();
		for (int i = 0; i < annotationFiles.size(); i++) {
			allAnnotations[i] = loader.loadAnnotations(annotationFiles.get(i), options.usePeakCenter, options.bpToExtendAnnotation, true, regions);

			// TODO: load in two column file with proper names
			String name = allAnnotations[i].getName();
			File file = new File(name);
			String filename = file.getName();
			allAnnotations[i].setName(filename);
		}

		// load posteriors
		BroShifterTask bs = new BroShifterTask();
		GTFAnnotation geneannotation = new GTFAnnotation(options.getGeneAnnotationFile());

		// iterate regions
		for (int r = 0; r < regions.size(); r++) {

			Feature region = regions.get(r);

			AssociationFile assocFile = new AssociationFile();
			ArrayList<AssociationResult> results = assocFile.read(options.posteriorFile, region);
			ArrayList<SNPFeature> snps = new ArrayList<>(results.size());
			ArrayList<Pair<Integer, Double>> allPvalues = new ArrayList<>();
			for (AssociationResult result : results) {
				Feature f = result.getSnp();
				if (region.overlaps(f)) {
					SNPFeature f2 = new SNPFeature();
					f2.setChromosome(f.getChromosome());
					f2.setStart(f.getStart());
					f2.setStop(f.getStop());
					f2.setName(f.getName());

					f2.setP(result.getPosterior());
					allPvalues.add(new Pair<Integer, Double>(f.getStart(), result.getLog10Pval()));
					snps.add(f2);
				}
			}

			ArrayList<Triple<Track, Double, Integer>> annotationsUnsorted = new ArrayList<>();


			for (int i = 0; i < allAnnotations.length; i++) {
				Track annot = allAnnotations[i].getSubset(region.getChromosome(), region.getStart(), region.getStop());
				Pair<Double, Integer> overlap = bs.getOverlap(annot, snps);
				Triple<Track, Double, Integer> t = new Triple<>(annot, overlap.getLeft(), overlap.getRight());
				annotationsUnsorted.add(t);
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

			// determine size of plot
			AnnotationPanel annotPanel = new AnnotationPanel(annotationFiles.size(), 1);
			annotPanel.setData(region, annotationsSorted);
			annotPanel.setOverlappingFeatures(snps);

			int pixelspertrack = annotPanel.getTrackheight() + annotPanel.getMarginBetween();
			int genepanelrows = (int) Math.ceil(100 / pixelspertrack);
			int assocPanelRows = (int) Math.ceil(300 / pixelspertrack);

			Grid grid = new Grid(1000, pixelspertrack, annotationFiles.size() + genepanelrows + assocPanelRows, 1, 100, 100);


			TreeSet<Gene> genes = geneannotation.getGeneTree();
			Gene geneStart = new Gene("", region.getChromosome(), Strand.POS, region.getStart(), region.getStart());
			Gene geneStop = new Gene("", region.getChromosome(), Strand.POS, region.getStop(), region.getStop());
			SortedSet<Gene> overlappingGenes = genes.subSet(geneStart, true, geneStop, true);

			ArrayList<Gene> overlappingGenesList = new ArrayList<>();
			overlappingGenesList.addAll(overlappingGenes);

			GenePanel genePanel = new GenePanel(genepanelrows, 1);
			genePanel.setData(region, overlappingGenesList);
			grid.addPanel(genePanel);

			AssociationPanel assocPanel = new AssociationPanel(assocPanelRows, 1);


			assocPanel.setDataSingleDs(region, null, allPvalues, region.toString());

			grid.addPanel(assocPanel);


			grid.addPanel(annotPanel);

			grid.draw(options.outfile + region.toString() + ".pdf");

		}

	}

	private class AssocTripleComparator implements Comparator<Triple<Track, Double, Integer>> {

		@Override
		public int compare(Triple<Track, Double, Integer> o1, Triple<Track, Double, Integer> o2) {
			if (o1.getMiddle() == null || Double.isNaN(o1.getMiddle()) || Double.isFinite(o1.getMiddle())) {
				return -1;
			}
			if (o2.getMiddle() == null || Double.isNaN(o2.getMiddle()) || Double.isFinite(o2.getMiddle())) {
				return 1;
			}
			if (o1.getMiddle() > o2.getMiddle()) {
				return 1;
			} else if (o1.getMiddle().equals(o2.getMiddle())) {
				return 0;
			} else {
				return -1;
			}
		}
	}

}
