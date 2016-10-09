package nl.harmjanwestra.gwas;

import com.itextpdf.text.DocumentException;
import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.gwas.CLI.AssociationPlotterOptions;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Strand;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.AssociationPanel;
import nl.harmjanwestra.utilities.graphics.panels.GenePanel;
import nl.harmjanwestra.utilities.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFTabix;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.*;

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
		String thresholdFile = options.getSignificanceThresholdFile();
		String ldrpefix = options.getLDPrefix();
		String ldlimit = options.getLDLimit();

		HashMap<Feature, Double> regionThresholds = new HashMap<Feature, Double>();
		if (thresholdFile != null) {
			TextFile tf = new TextFile(thresholdFile, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				if (elems.length >= 3) {
					Feature f = Feature.parseFeature(elems[0]);
					double d = Double.parseDouble(elems[2]);
					regionThresholds.put(f, d);
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();

		}

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
			Double threshold = regionThresholds.get(region);
			if (threshold == null) {
				threshold = 5E-8;
			}
//			Collection<Gene> allgenes = annotation.getGenes();
//			for (Gene g : allgenes) {
//				String geneid = g.getGeneId();
//				if(geneid.equals("PUS10")){
//					System.out.println(g.toString());
//				}
//
//			}


//			for (Gene g : genes) {
//				String geneid = g.getGeneId();
//				if(geneid.equals("PUS10")){
//					System.out.println(g.toString());
//				}
//
//			}

			TreeSet<Gene> genes = annotation.getGeneTree();
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
					ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(associations, options.getCredibleSetThreshold());
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
					posteriorPanel.setMarkDifferentShape(mark);
					posteriorPanel.setMaxPVal(0.99d);
					grid.addPanel(posteriorPanel, 2, i);
				}

				AssociationPanel associationPanel = new AssociationPanel(1, 1);
				ArrayList<Pair<Integer, Double>> pvals = new ArrayList<Pair<Integer, Double>>();

				System.out.println(associations.size() + " pvals loaded");


				ArrayList<String> variants = new ArrayList<String>();
				String maxVar = null;

				for (int a = 0; a < associations.size(); a++) {
					AssociationResult r = associations.get(a);
					double p = r.getLog10Pval();
					if (!Double.isNaN(p) && !Double.isInfinite(p)) {
						pvals.add(new Pair<>(r.getSnp().getStart(), r.getLog10Pval()));
						variants.add(r.getSnp().getStart() + "_" + r.getSnp().getName());
						if (maxP == null) {
							maxP = r.getLog10Pval();
							maxVar = r.getSnp().getStart() + "_" + r.getSnp().getName();
						} else if (r.getLog10Pval() > maxP) {
							maxP = r.getLog10Pval();
							maxVar = r.getSnp().getStart() + "_" + r.getSnp().getName();
						}
					} else {
						System.err.println("issue with: " + r.toString());
					}
				}

				double[] ldData = null;
				if (ldrpefix != null) {
					ldData = new double[pvals.size()];
					HashMap<String, Integer> variantsPresentIndex = new HashMap<String, Integer>();
					VCFVariant[] variantArr = new VCFVariant[pvals.size()];
					for (int v = 0; v < variants.size(); v++) {
						variantsPresentIndex.put(variants.get(v), v);
					}

					String tabixfile = ldrpefix.replaceAll("CHR", "" + region.getChromosome().getNumber());
					boolean[] sampleLimit = null;
					if (ldlimit != null) {
						sampleLimit = VCFTabix.getSampleFilter(tabixfile, ldlimit);
					}
					TabixReader.Iterator window = VCFTabix.query(tabixfile, region);

					String next = window.next();
					while (next != null) {
						VCFVariant variant = new VCFVariant(next, VCFVariant.PARSE.HEADER);
						Integer index = variantsPresentIndex.get(variant.getPos() + "_" + variant.getId());
						if (index != null) {
							variantArr[index] = new VCFVariant(next, VCFVariant.PARSE.ALL, sampleLimit);
						}
					}

					Integer maxVarIndex = variantsPresentIndex.get(maxVar);
					ldData[maxVarIndex] = 1d;
					VCFVariant topVariant = variantArr[maxVarIndex];
					DetermineLD ldcalc = new DetermineLD();
					for (int v = 0; v < variantArr.length; v++) {
						if (!maxVarIndex.equals(v)) {
							Pair<Double, Double> ld = ldcalc.getLD(variantArr[v], topVariant);
							ldData[v] = ld.getRight();
						}
					}
					associationPanel.setLDData(ldData);
				}

				if (!pvals.isEmpty()) {
					regionhasvariants = true;
				}

				associationPanel.setDataSingleDs(region, sequencedRegions, pvals, assocNames[i] + " Association P-values");
				associationPanel.setMarkDifferentShape(mark);
				associationPanel.setPlotGWASSignificance(true, threshold);
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
