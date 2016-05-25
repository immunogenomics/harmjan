package nl.harmjanwestra.gwas;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.gwas.CLI.AssociationPlotterOptions;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.features.Strand;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.AssociationPanel;
import nl.harmjanwestra.utilities.graphics.panels.GenePanel;
import nl.harmjanwestra.utilities.gtf.GTFAnnotation;
import umcg.genetica.containers.Pair;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Created by Harm-Jan on 01/13/16.
 */
public class AssociationPosteriorPlotter {

	AssociationPlotterOptions options;

	public AssociationPosteriorPlotter(AssociationPlotterOptions options) throws IOException, DocumentException {
		this.options = options;
		plotAssociationsAndPosteriors();
	}

	public void plotAssociationsAndPosteriors() throws IOException, DocumentException {

		String associationFiles = options.getAssociationFiles();
		String associationFileNames = options.getAssociationFileNames();
		String annotationfile = options.getAnnotationfile();
		String bedregionfile = options.getBedregionfile();
		String outputprefix = options.getOutputprefix();
		String sequencedRegionsFile = options.getSequencedRegionsFile();
		boolean plotPosteriors = options.isPlotPosterior();

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bedregionfile);
		ArrayList<Feature> sequencedRegionsList = reader.readAsList(sequencedRegionsFile);
		HashSet<Feature> sequencedRegions = new HashSet<Feature>();
		sequencedRegions.addAll(sequencedRegionsList);

		String[] assocNames = associationFileNames.split(",");
		String[] assocFiles = associationFiles.split(",");
		AssociationFile assocFile = new AssociationFile();
		GTFAnnotation annotation = new GTFAnnotation(annotationfile);

		for (Feature region : regions) {
			boolean regionhasvariants = false;

//			Collection<Gene> allgenes = annotation.getGenes();
//			for (Gene g : allgenes) {
//				String geneid = g.getGeneId();
//				if(geneid.equals("PUS10")){
//					System.out.println(g.toString());
//				}
//
//			}

			TreeSet<Gene> genes = annotation.getGeneTree();
//			for (Gene g : genes) {
//				String geneid = g.getGeneId();
//				if(geneid.equals("PUS10")){
//					System.out.println(g.toString());
//				}
//
//			}


			Gene geneStart = new Gene("", region.getChromosome(), Strand.POS, region.getStart(), region.getStart());
			Gene geneStop = new Gene("", region.getChromosome(), Strand.POS, region.getStop(), region.getStop());
			SortedSet<Gene> overlappingGenes = genes.subSet(geneStart, true, geneStop, true);

			for (Gene g : overlappingGenes) {
				System.out.println(g.toString());

			}

			ArrayList<Gene> overlappingGenesList = new ArrayList<>();
			overlappingGenesList.addAll(overlappingGenes);

			Grid grid = new Grid(400, 300, 2, assocFiles.length, 100, 100);
			if (plotPosteriors) {
				grid = new Grid(400, 300, 3, assocFiles.length, 100, 100);
			}

			GenePanel genePanel = new GenePanel(1, 1);
			genePanel.setData(region, overlappingGenesList);
			for (int i = 0; i < assocFiles.length; i++) {
				grid.addPanel(genePanel, 0, i);
			}

			ArrayList<AssociationPanel> allPanels = new ArrayList<>();
			Double maxP = null;
			for (int i = 0; i < assocFiles.length; i++) {

				System.out.println("Reading: " + assocFiles[i]);
				ArrayList<AssociationResult> associations = assocFile.read(assocFiles[i], region);
				HashSet<AssociationResult> credibleSetSet = new HashSet<>();
				boolean[] mark = null;
				if (plotPosteriors) {
					AssociationPanel posteriorPanel = new AssociationPanel(1, 1);
					ArrayList<Pair<Integer, Double>> posteriors = new ArrayList<Pair<Integer, Double>>();
					ApproximateBayesPosterior abp = new ApproximateBayesPosterior();
					ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(associations, 0.95);
					credibleSetSet.addAll(credibleSet);
					mark = new boolean[associations.size()];
					for (int a = 0; a < associations.size(); a++) {
						AssociationResult r = associations.get(a);
						posteriors.add(new Pair<>(r.getSnp().getStart(), r.getPosterior()));
						if (credibleSetSet.contains(r)) {
							mark[a] = true;
						}
					}
					posteriorPanel.setDataSingleDs(region, sequencedRegions, posteriors, "Posteriors");
					posteriorPanel.setMarkDifferentColor(mark);
					posteriorPanel.setMaxPVal(0.99d);
					grid.addPanel(posteriorPanel, 2, i);
				}

				AssociationPanel associationPanel = new AssociationPanel(1, 1);
				ArrayList<Pair<Integer, Double>> pvals = new ArrayList<Pair<Integer, Double>>();

				for (int a = 0; a < associations.size(); a++) {
					AssociationResult r = associations.get(a);

					double p = r.getLog10Pval();
					if (!Double.isNaN(p) && !Double.isInfinite(p)) {
						pvals.add(new Pair<>(r.getSnp().getStart(), r.getLog10Pval()));

						if (maxP == null) {
							maxP = r.getLog10Pval();
						} else if (r.getLog10Pval() > maxP) {
							maxP = r.getLog10Pval();
						}
					} else {
						System.err.println("issue with: " + r.toString());
					}


				}

				if (!pvals.isEmpty()) {
					regionhasvariants = true;
				}
				associationPanel.setDataSingleDs(region, sequencedRegions, pvals, assocNames[i] + " Association P-values");
				associationPanel.setMarkDifferentColor(mark);
				associationPanel.setPlotGWASSignificance(true);
				allPanels.add(associationPanel);
				grid.addPanel(associationPanel, 1, i);
			}

			if (regionhasvariants) {
				for (AssociationPanel p : allPanels) {
					System.out.println("plotting: " + region.toString() + "\tmax pval: " + maxP);
					if (options.getMaxp() != null) {
						p.setMaxPVal(options.getMaxp());
					} else if (maxP != null) {
						p.setMaxPVal(maxP);
					}
				}
				grid.draw(outputprefix + region.toString() + ".pdf");
			} else {
				System.out.println("Region not plotted: " + region.toString() + " since it has no variants.");
			}

		}


	}

}
