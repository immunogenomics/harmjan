package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.features.Strand;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.GenePanel;
import nl.harmjanwestra.utilities.gtf.GTFAnnotation;

import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Created by hwestra on 7/16/15.
 */
public class PlotTest {

	public static void main(String[] args) {

		Grid grid = new Grid(250, 250, 2, 2, 50, 50);


//		for (int i = 0; i < 3; i++) {
//			double[][] data = new double[2][10];
//			HistogramPanel.DATASETPLOTTYPE[] types = new HistogramPanel.DATASETPLOTTYPE[2];
//
//			for (int q = 0; q < data.length; q++) {
//				if (q == 1) {
//					types[q] = HistogramPanel.DATASETPLOTTYPE.POLY;
//				} else {
//					types[q] = HistogramPanel.DATASETPLOTTYPE.BAR;
//				}
//			}
//
//			for (int j = 0; j < 10; j++) {
//				for (int q = 0; q < data.length; q++) {
//					data[q][j] = j;
//				}
//			}
//
//			HistogramPanel histogram = new HistogramPanel(1, 1);
//			if (i == 1) {
//				histogram = new HistogramPanel(2, 2);
//			}
//			histogram.setData(data);
//
//			histogram.setMarginX(50);
//			histogram.setMarginY(50);
//
//			if (i == 0) {
//				histogram.setPlotRange(new Range(-10, 0, 10, 20));
//			}
//
//			histogram.setDatasetPlotTypes(types);
//
//			grid.addPanel(histogram);
//
//		}
//
//		try {
//			grid.draw("/Data/tmp/plottest.pdf");
//		} catch (IOException e) {
//			e.printStackTrace();
//		} catch (DocumentException e) {
//			e.printStackTrace();
//		}

		try {

			Feature region = new Feature();
			// Chr1_113931006-114495318
			region.setChromosome(Chromosome.ONE);
			region.setStart(113931006);
			region.setStop(114495318);

			String gtf = "/Data/Annotation/UCSC/genes.gtf";
			GTFAnnotation annot = new GTFAnnotation(gtf);

			TreeSet<Gene> genes = annot.getGeneTree();
			Gene geneStart = new Gene("", region.getChromosome(), Strand.POS, region.getStart(), region.getStart());
			Gene geneStop = new Gene("", region.getChromosome(), Strand.POS, region.getStop(), region.getStop());
			SortedSet<Gene> overlappingGenes = genes.subSet(geneStart, true, geneStop, true);

			grid = new Grid(500, 250, 1, 1, 50, 50);

			GenePanel genePanel = new GenePanel(1, 1);
			genePanel.setGenes(overlappingGenes);
			genePanel.setRegion(region);

			grid.addPanel(genePanel);
			grid.draw("/Data/tmp/plottest-genes.pdf");

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
