package nl.harmjanwestra.gwas;

import com.itextpdf.text.DocumentException;
import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.gwas.CLI.AssociationPlotterOptions;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
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


	public static void main(String[] args) {
		String[] arguments = new String[]{
				"--plotposteriors",
				"-a", "/Data/Ref/Annotation/UCSC/genes.gtf",
				"-i",
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/ConditionalOnMeta/META-assoc0.3-COSMO-iter2-merged.txt.gz," +
						"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/ConditionalOnMeta/T1D-assoc0.3-COSMO-iter2-merged.txt.gz," +
						"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/ConditionalOnMeta/RA-assoc0.3-COSMO-iter2-merged.txt.gz",
				"-n", "META,T1D,RA",
				"-o", "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/ConditionalOnMeta/cd28/iter2",
				"-r", "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/ctla4.bed",
				"--ldprefix", "/Data/Ref/beagle_1kg/1kg.phase3.v5a.chrCHR.vcf.gz",
				"--ldlimit", "/Data/Ref/1kg-europeanpopulations.txt.gz",
				"--thresholds", "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/BonferroniThresholds/META.txt," +
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/BonferroniThresholds/T1D.txt," +
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/BonferroniThresholds/RA.txt,"

		};

//		String[] arguments = new String[]{
//				"--plotposteriors",
//				"-a", "/Data/Ref/Annotation/UCSC/genes.gtf",
//				"-i",
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/ConditionalOnMeta/META-assoc0.3-COSMO-merged-posterior.txt.gz," +
//						"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/ConditionalOnMeta/T1D-assoc0.3-COSMO-merged-posterior.txt.gz," +
//						"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/ConditionalOnMeta/RA-assoc0.3-COSMO-merged-posterior.txt.gz",
//				"-n", "META,T1D,RA",
//				"-o", "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/ConditionalOnMeta/cd28/iter0",
//				"-p",
//				"-r", "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/ctla4.bed",
//				"--ldprefix", "/Data/Ref/beagle_1kg/1kg.phase3.v5a.chrCHR.vcf.gz",
//				"--ldlimit", "/Data/Ref/1kg-europeanpopulations.txt.gz",
////				"--thresholds", "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/BonferroniThresholds/RA.txt"
//		};

		//
		AssociationPlotterOptions options = new AssociationPlotterOptions(arguments);
		try {
			AssociationPosteriorPlotter p = new AssociationPosteriorPlotter(options);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
	}

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

		ArrayList<HashMap<Feature, Double>> regionThresholds = new ArrayList<HashMap<Feature, Double>>();
		if (thresholdFile != null) {
			String[] thresholdFiles = thresholdFile.split(",");
			for (int s = 0; s < thresholdFiles.length; s++) {
				System.out.println("Loading significance thresholds: " + thresholdFile);
				TextFile tf = new TextFile(thresholdFiles[s], TextFile.R);
				HashMap<Feature, Double> regionThresholdsHash = new HashMap<Feature, Double>();
				String[] elems = tf.readLineElems(TextFile.tab);
				while (elems != null) {
					if (elems.length >= 3) {
						String[] felems = elems[0].split("_");
						String[] posElems = felems[1].split("-");
						Chromosome chr = Chromosome.parseChr(felems[0]);
						Integer start = Integer.parseInt(posElems[0]);
						Integer stop = Integer.parseInt(posElems[1]);
						Feature f = new Feature(chr, start, stop);
						double d = Double.parseDouble(elems[2]);
						regionThresholdsHash.put(f, d);
					}
					elems = tf.readLineElems(TextFile.tab);
				}
				tf.close();
				regionThresholds.add(regionThresholdsHash);
				System.out.println(regionThresholdsHash.size() + " thresholds loaded.");
			}

		}

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bedregionfile);

		ArrayList<Feature> sequencedRegionsList = null;
		HashSet<Feature> sequencedRegions = null;
		if (sequencedRegionsFile != null) {
			sequencedRegionsList = reader.readAsList(sequencedRegionsFile);
			sequencedRegions = new HashSet<Feature>();
			sequencedRegions.addAll(sequencedRegionsList);
		}


		String[] assocNames = associationFileNames.split(",");
		String[] assocFiles = associationFiles.split(",");
		AssociationFile assocFile = new AssociationFile();
		GTFAnnotation annotation = new GTFAnnotation(annotationfile);


		for (Feature region : regions) {
			boolean regionhasvariants = false;

			TreeSet<Gene> genes = annotation.getGeneTree();
			Gene geneStart = new Gene("", region.getChromosome(), Strand.POS, region.getStart(), region.getStart());
			Gene geneStop = new Gene("", region.getChromosome(), Strand.POS, region.getStop(), region.getStop());
			SortedSet<Gene> overlappingGenes = genes.subSet(geneStart, true, geneStop, true);

			for (Gene g : overlappingGenes) {
				System.out.println(g.toString());
			}

			ArrayList<Gene> overlappingGenesList = new ArrayList<>();
			overlappingGenesList.addAll(overlappingGenes);

			int gridrows = 2;
			if (plotPosteriors) {
				gridrows = 3;
			}
			Grid grid = new Grid(200, 100, gridrows, assocFiles.length, 100, 100);

			GenePanel genePanel = new GenePanel(1, 1);
			genePanel.setData(region, overlappingGenesList);
			for (int i = 0; i < assocFiles.length; i++) {
				grid.addPanel(genePanel, 0, i);
			}

			ArrayList<AssociationPanel> allPanels = new ArrayList<>();
			Double maxP = null;
			for (int i = 0; i < assocFiles.length; i++) {
				Double threshold = null;
				if (regionThresholds != null) {
					threshold = regionThresholds.get(i).get(new Feature(region.getChromosome(), region.getStart(), region.getStop()));
					System.out.println(threshold + " for region: " + region.toString());
				} else {
					threshold = 5E-8;
				}

				Double maxPDs = null;
				String maxVar = null;
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


				for (int a = 0; a < associations.size(); a++) {
					AssociationResult r = associations.get(a);
					double p = r.getLog10Pval();
					if (!Double.isNaN(p) && !Double.isInfinite(p)) {
						pvals.add(new Pair<>(r.getSnp().getStart(), r.getLog10Pval()));
						variants.add("" + r.getSnp().getStart());
						if (maxP == null) {
							maxP = r.getLog10Pval();
						} else if (r.getLog10Pval() > maxP) {
							maxP = r.getLog10Pval();
						}
						if (maxPDs == null) {
							maxPDs = r.getLog10Pval();
							maxVar = "" + r.getSnp().getStart();
						} else if (r.getLog10Pval() > maxPDs) {
							maxPDs = r.getLog10Pval();
							maxVar = "" + r.getSnp().getStart();
						}
					} else {
						System.err.println("issue with: " + r.toString());
					}
				}

				System.out.println(maxVar + " is the max var.");

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
					int found = 0;
					HashSet<String> variantsFound = new HashSet<String>();
					while (next != null) {
						VCFVariant variant = new VCFVariant(next, VCFVariant.PARSE.HEADER);

						Integer index = variantsPresentIndex.get("" + variant.getPos());
						if (index != null) {
							variantArr[index] = new VCFVariant(next, VCFVariant.PARSE.ALL, sampleLimit);
							found++;
							variantsFound.add("" + variant.getPos());
						}
						next = window.next();
					}


					int notfound = 0;
					for (String snp : variantsPresentIndex.keySet()) {
						if (!variantsFound.contains(snp)) {
							System.out.println("Could not find: " + snp);
							notfound++;
						}
					}
					System.out.println(found + " variants found in LD reference");
					System.out.println(notfound + " variants not in LD reference?");
//					System.exit(-1);


					Integer maxVarIndex = variantsPresentIndex.get(maxVar);
					if (maxVarIndex == null) {
						System.out.println(maxVar + " not found for dataset " + assocFiles[i]);
						System.exit(-1);

					}
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
