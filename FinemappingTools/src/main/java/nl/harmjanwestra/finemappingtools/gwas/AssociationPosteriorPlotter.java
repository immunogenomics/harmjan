package nl.harmjanwestra.finemappingtools.gwas;

import com.itextpdf.text.DocumentException;
import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.finemappingtools.gwas.CLI.AssociationPlotterOptions;
import nl.harmjanwestra.utilities.annotation.Annotation;
import nl.harmjanwestra.utilities.annotation.ensembl.EnsemblStructures;
import nl.harmjanwestra.utilities.annotation.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.AssociationPanel;
import nl.harmjanwestra.utilities.graphics.panels.GenePanel;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFTabix;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
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

	public static void main(String[] args) {
		String[] arguments = new String[]{
				"--plotposteriors",
				"-a", "/Data/Ref/Annotation/UCSC/genes.gtf",
				"-i",
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/ConditionalOnMeta/META-assoc0.3-COSMO-iter1-merged.txt.gz," +
						"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/ConditionalOnMeta/T1D-assoc0.3-COSMO-iter1-merged.txt.gz," +
						"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/ConditionalOnMeta/RA-assoc0.3-COSMO-iter1-merged.txt.gz",
				"-n", "META,T1D,RA",
				"-o", "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/testcond",
				"-r", "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/TNFAIP3.bed",
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
//				"-o", "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/tnfaip3",
//				"-p",
//				"-r", "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/TNFAIP3.bed",
//				"--ldprefix", "/Data/Ref/beagle_1kg/1kg.phase3.v5a.chrCHR.vcf.gz",
//				"--ldlimit", "/Data/Ref/1kg-europeanpopulations.txt.gz",
//				"--thresholds", "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/BonferroniThresholds/META.txt," +
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/BonferroniThresholds/T1D.txt," +
//				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/BonferroniThresholds/RA.txt,"
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

	private ArrayList<HashMap<Feature, Double>> loadThresholds(String file, int datasets) throws IOException {
		ArrayList<HashMap<Feature, Double>> regionSignificanceThresholds = new ArrayList<HashMap<Feature, Double>>();
		String[] thresholdFiles = file.split(",");
		if (thresholdFiles.length < datasets) {
			thresholdFiles = new String[datasets];
			for (int i = 0; i < thresholdFiles.length; i++) {
				thresholdFiles[i] = file;
			}
		}
		for (int s = 0; s < thresholdFiles.length; s++) {
			System.out.println("Loading significance thresholds: " + thresholdFiles[s]);
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
					System.out.println(f.toString() + "\t" + d);
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			regionSignificanceThresholds.add(regionThresholdsHash);
			System.out.println(regionThresholdsHash.size() + " thresholds loaded.");

		}
		return regionSignificanceThresholds;
	}

	public void plotAssociationsAndPosteriors() throws IOException, DocumentException {

		String associationFiles = options.getAssociationFiles();
		String associationFileNames = options.getAssociationFileNames();
		String annotationfile = options.getAnnotationfile();
		String bedregionfile = options.getBedregionfile();
		String outputprefix = options.getOutputprefix();
		String sequencedRegionsFile = options.getSequencedRegionsFile();
		boolean plotPosteriors = options.isPlotPosterior();

		String significanceThresholdFile = options.getSignificanceThresholdFile();
		String significanceConditionalThresholdFile = options.getSignificanceConditionalThresholdFile();

		String maxpvalfile = options.getMaxPvalueFile();
		String maxpvalconditonionalfile = options.getMaxConditionalPvalueFile();

		String ldrpefix = options.getLDPrefix();
		String ldlimit = options.getLDLimit();
		double defaultsignificancethreshold = options.getDefaultSignificance();
		double defaultsignificancethresholdconditional = options.getDefaultSignificanceConditional();

		int nriters = options.getNrIters();


		String functionalAnnotation = null;

		String conditionalFiles = options.getConditionalFiles();


		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bedregionfile);
		HashSet<Feature> uniqueRegions = new HashSet<>();
		uniqueRegions.addAll(regions);
		regions = new ArrayList<>();
		regions.addAll(uniqueRegions);

		ArrayList<Feature> sequencedRegionsList = null;
		HashSet<Feature> sequencedRegions = null;
		if (sequencedRegionsFile != null) {
			sequencedRegionsList = reader.readAsList(sequencedRegionsFile);
			sequencedRegions = new HashSet<Feature>();
			sequencedRegions.addAll(sequencedRegionsList);
		}


		String[] assocNames = associationFileNames.split(",");
		String[] assocFiles = associationFiles.split(",");
		String[] assocConditionalFiles = null;
		if (conditionalFiles != null) {
			assocConditionalFiles = conditionalFiles.split(",");
		}

		AssociationFile assocFile = new AssociationFile();
		Annotation annotation = null;
		if (annotationfile.endsWith(".gtf.gz") || annotationfile.endsWith(".gtf")) {
			annotation = new GTFAnnotation(annotationfile);
		} else {
			annotation = new EnsemblStructures(annotationfile);
		}

		ArrayList<HashMap<Feature, Double>> regionSignificanceThresholds = null;
		ArrayList<HashMap<Feature, Double>> regionConditionalSignificanceThresholds = null;
		ArrayList<HashMap<Feature, Double>> regionMaxPThresholds = null;
		ArrayList<HashMap<Feature, Double>> regionMaxConditionalPThresholds = null;

		if (significanceThresholdFile != null && Gpio.exists(significanceThresholdFile)) {
			regionSignificanceThresholds = loadThresholds(significanceThresholdFile, assocNames.length);
		}
		if (significanceConditionalThresholdFile != null && Gpio.exists(significanceConditionalThresholdFile)) {
			regionConditionalSignificanceThresholds = loadThresholds(significanceConditionalThresholdFile, assocNames.length);
		}
		if (maxpvalfile != null && Gpio.exists(maxpvalfile)) {
			System.out.println("Loading max pvals: " + maxpvalfile);
			// make this dependent on iteration?
			regionMaxPThresholds = loadThresholds(maxpvalfile, assocNames.length);
		}

		if (maxpvalconditonionalfile != null && Gpio.exists(maxpvalconditonionalfile)) {
			regionMaxConditionalPThresholds = loadThresholds(maxpvalconditonionalfile, assocNames.length);
		}


		for (Feature region : regions) {
			boolean regionhasvariants = false;

			TreeSet<Gene> genes = annotation.getGeneTree();

			ArrayList<Gene> overlappingGenesList = new ArrayList<>();
			for (Gene g : genes) {
				if (g.overlaps(region)) {
					System.out.println(g.toString());
					overlappingGenesList.add(g);
				}
			}

			System.out.println();

			int gridrows = 2;
			if (plotPosteriors) {
				gridrows = 3;
			}
			if (assocConditionalFiles != null) {
				System.out.println("Adding " + nriters + " extra rows for conditional output..");
				gridrows += (nriters - 1);
			}

			Grid grid = new Grid(200, 100, gridrows, assocFiles.length, 100, 100);

			GenePanel genePanel = new GenePanel(1, 1);
			genePanel.setData(region, overlappingGenesList);
			for (int i = 0; i < assocFiles.length; i++) {
				grid.addPanel(genePanel, 0, i);
			}

			ArrayList<ArrayList<AssociationPanel>> allPanels = new ArrayList<>();

			for (int datasetNr = 0; datasetNr < assocFiles.length; datasetNr++) {
				ArrayList<AssociationPanel> panelsForDs = new ArrayList<>();
				Double threshold = null;


				for (int iter = 0; iter < nriters; iter++) {

					ArrayList<HashMap<Feature, Double>> thresholdHashToUse = regionSignificanceThresholds;
					threshold = defaultsignificancethreshold;
					if (iter != 0) {
						thresholdHashToUse = regionConditionalSignificanceThresholds;
						threshold = defaultsignificancethresholdconditional;
					}
					if (thresholdHashToUse != null) {
						threshold = thresholdHashToUse.get(datasetNr).get(new Feature(region.getChromosome(), region.getStart(), region.getStop()));
						System.out.println(threshold + " for region: " + region.toString());
					}

					String assocFileToRead = assocFiles[datasetNr];
					if (iter > 0) {
						String tmpfilename = assocConditionalFiles[datasetNr];
						tmpfilename = tmpfilename.replaceAll("ITER", "" + iter);
						assocFileToRead = tmpfilename;
					}

					Double maxPDs = null;
					String maxVar = null;
					System.out.println("Reading: " + assocFiles[datasetNr]);
					ArrayList<AssociationResult> associations = new ArrayList<>();
					if (Gpio.exists(assocFileToRead)) {
						associations = assocFile.readRegion(assocFileToRead, region);
					}
					HashSet<AssociationResult> credibleSetSet = new HashSet<>();
					boolean[] mark = null;
					if (plotPosteriors && iter == 0) {
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

						posteriorPanel.setMaxPVal(options.getPosteriorThreshold());
						grid.addPanel(posteriorPanel, 1, datasetNr);
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
					if (!pvals.isEmpty()) {
						regionhasvariants = true;
					}
					double[] ldData = null;
					if (ldrpefix != null && regionhasvariants) {
						ldData = new double[pvals.size()];
						HashMap<String, Integer> variantsPresentIndex = new HashMap<String, Integer>();
						VCFVariant[] variantArr = new VCFVariant[pvals.size()];
						for (int v = 0; v < variants.size(); v++) {
							variantsPresentIndex.put(variants.get(v), v);
						}

						String tabixfile = ldrpefix.replaceAll("CHR", "" + region.getChromosome().getNumber());
						VCFTabix tabix = new VCFTabix(tabixfile);
						boolean[] sampleLimit = null;
						if (ldlimit != null) {
							sampleLimit = tabix.getSampleFilter(ldlimit);
						}
						TabixReader.Iterator window = tabix.query(region);

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
						tabix.close();

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
							System.out.println(maxVar + " not found for dataset " + assocFiles[datasetNr]);
//						System.exit(-1);
						} else {
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
					}

					associationPanel.setDataSingleDs(region, sequencedRegions, pvals, assocNames[datasetNr] + " Association P-values");
					associationPanel.setMarkDifferentShape(mark);
					associationPanel.setPlotGWASSignificance(true, threshold);
					panelsForDs.add(associationPanel);

					if (plotPosteriors) {
						grid.addPanel(associationPanel, 2 + iter, datasetNr);
					} else {
						grid.addPanel(associationPanel, 1 + iter, datasetNr);
					}
				}
				allPanels.add(panelsForDs);
			}

			if (regionhasvariants) {

				// set the max p for each row
				for (int iter = 0; iter < nriters; iter++) {
					double maxP = 0;
					for (int datasetNr = 0; datasetNr < allPanels.size(); datasetNr++) {
						AssociationPanel p = allPanels.get(datasetNr).get(iter);
						ArrayList<ArrayList<Pair<Integer, Double>>> allPvalues = p.getAllPValues();
						for (ArrayList<Pair<Integer, Double>> a : allPvalues) {
							for (Pair<Integer, Double> v : a) {
								if (v.getRight() > maxP) {
									maxP = v.getRight();
								}
							}
						}
					}
					System.out.println("iter:\t" + iter + "\tmaxp: " + maxP);
					for (int datasetNr = 0; datasetNr < allPanels.size(); datasetNr++) {
						AssociationPanel p = allPanels.get(datasetNr).get(0);
						p.setMaxPVal(maxP);
					}
				}

				for (int iter = 0; iter < nriters; iter++) {
					for (int datasetNr = 0; datasetNr < allPanels.size(); datasetNr++) {
						AssociationPanel p = allPanels.get(datasetNr).get(iter);
						ArrayList<HashMap<Feature, Double>> regionmaxthresholdsToUse = regionMaxPThresholds;
						if (iter != 0) {
							System.out.println("Using regional thresholds");
							regionmaxthresholdsToUse = regionMaxConditionalPThresholds;
							System.out.println("" + (regionmaxthresholdsToUse != null));
//							System.exit(-1);
						}

						if (regionmaxthresholdsToUse != null) {
							Double pval = regionmaxthresholdsToUse.get(datasetNr).get(region.newFeatureFromCoordinates());
							if (pval == null) {
								System.out.println("Could not find locus: " + region.toString());
							} else {
								System.out.println("Iter\t" + iter + "\tLocus p: " + region.toString() + "\t" + pval);
							}

							if (pval != null) {
								p.setMaxPVal(-Math.log10(pval));
								p.setRoundUpYAxis(false);
							}
						}

						if (options.getMaxp() != null) {
							p.setMaxPVal(options.getMaxp());
						}
					}
				}
				System.out.println("plotting: " + region.toString() + "\t" + outputprefix + region.toString() + ".pdf");
				grid.draw(outputprefix + region.toString() + ".pdf");
			} else {
				System.out.println("Region not plotted: " + region.toString() + " since it has no variants.");
			}

//			System.exit(-1);
		}


	}

}
