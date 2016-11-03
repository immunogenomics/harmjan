package nl.harmjanwestra.finemapping;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Strand;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.gtf.GTFAnnotation;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Created by hwestra on 8/19/16.
 */
public class Bonferroninator {

	public static void main(String[] args) {

		String assocfile = "";
		String allregionsfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci-woMHC.txt";
		String diseasespecificregionfile = "";
		String out = "";

		// RA
		Bonferroninator b = new Bonferroninator();
		try {

			String[] significantRegionFiles = new String[]{
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/RA-significantloci-75e7.bed",
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/T1D-significantloci-75e7.bed",
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/META-significantloci-75e7.bed",
			};
			String[] assocFiles = new String[]{
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/Old/ConditionalOnMeta/RA-assoc0.3-COSMO-iter0-merged.txt.gz",
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/Old/ConditionalOnMeta/T1D-assoc0.3-COSMO-iter0-merged.txt.gz",
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/Old/ConditionalOnMeta/META-assoc0.3-COSMO-iter0-merged.txt.gz"
			};

			String[] bonferroniOut = new String[]{
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/Old/ConditionalOnMeta/RA.txt",
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/Old/ConditionalOnMeta/T1D.txt",
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/Old/ConditionalOnMeta/META.txt"
			};
			String[] diseaseRegions = new String[]{
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseRALoci.bed",
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DLoci.bed",
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci-woMHC.txt"
			};

			double defaultthreshold = 7.5E-7;
			double significantthreshold = b.determineTopNumberOfVariantsWithinSignificantRegions(significantRegionFiles, assocFiles);

			System.out.println();
			System.out.println("---");
			System.out.println();

			for (int a = 0; a < assocFiles.length; a++) {
				assocfile = assocFiles[a];
				diseasespecificregionfile = diseaseRegions[a];
				// b.determineRegionSignificanceThresholds(assocfile, allregionsfile, diseasespecificregionfile, bonferroniOut[a], defaultthreshold, significantthreshold);
			}

			String[] diseases = new String[]{
					"RA",
					"T1D",
					"META"
			};


			String annot = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/genes.gtf.gz";
			GTFAnnotation annotation = new GTFAnnotation(annot);
			TreeSet<Gene> genes = annotation.getGeneTree();
			int maxiter = 3;
			for (int d = 0; d < diseases.length; d++) {

				HashMap<String, Double> thresholds = new HashMap<String, Double>();

				TextFile tf = new TextFile(bonferroniOut[d], TextFile.R);
				String[] elems = tf.readLineElems(TextFile.tab);
				while (elems != null) {
					thresholds.put(elems[0], Double.parseDouble(elems[2]));
					elems = tf.readLineElems(TextFile.tab);
				}
				tf.close();

				String output = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/Conditional/" + diseases[d] + "-bestAssocPerRegion.txt";
				TextFile outtf = new TextFile(output, TextFile.W);
				outtf.writeln("Iter\tRegion\tGenes\tmaxVariant\tPval\tLog10Pval\tGlobalThreshold\tGlobalSignificant\tNrVariantsInRegion\tLocalThreshold\tLocalSignificant");
				for (int iter = 0; iter < maxiter; iter++) {
					System.out.println(iter + "\t" + d);
					String assoc = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/Conditional/" + diseases[d] + "-assoc0.3-COSMO-gwas-" + iter + "-merged.txt.gz";
					String regions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
					b.determineTopAssocPerRegion(iter, regions, defaultthreshold, assoc, thresholds, genes, outtf);
				}
				outtf.close();

			}


		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void determineTopAssocPerRegion(int a, String regionfile, double defaultthreshold,
											String assoc, HashMap<String, Double> bonferroniThresholds, TreeSet<Gene> geneset, TextFile out) throws IOException {


		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> allRegions = reader.readAsList(regionfile);
		AssociationFile assocF = new AssociationFile();
		for (int r = 0; r < allRegions.size(); r++) {
			Feature region = allRegions.get(r);

			Gene geneStart = new Gene("", region.getChromosome(), Strand.POS, region.getStart(), region.getStart());
			Gene geneStop = new Gene("", region.getChromosome(), Strand.POS, region.getStop(), region.getStop());
			SortedSet<Gene> overlappingGenes = geneset.subSet(geneStart, true, geneStop, true);

			String[] genesArr = new String[overlappingGenes.size()];
			int ctr = 0;
			for (Gene g : overlappingGenes) {
				genesArr[ctr] = g.getName();
				ctr++;
			}

			ArrayList<AssociationResult> results = assocF.read(assoc, region);

			// get top result
			AssociationResult max = null;
			for (AssociationResult result : results) {
				if (max == null) {
					max = result;
				} else if (result.getLog10Pval() > max.getLog10Pval()) {
					max = result;
				}
			}

			// determine if it is significant
			if (max != null) {
				String regionStr = region.toString();

				Double threshold = bonferroniThresholds.get(regionStr);
				if (a == 0) {
					threshold = defaultthreshold;
				}
				boolean significant = false;
				if (max.getPval() < threshold) {
					// write to disk
					significant = true;
				}
				out.writeln(a +
						"\t" + region.toString()
						+ "\t" + Strings.concat(genesArr, Strings.semicolon)
						+ "\t" + max.getSnp().toString()
						+ "\t" + max.getPval()
						+ "\t" + max.getLog10Pval()
						+ "\t" + threshold
						+ "\t" + significant
						+ "\t" + results.size()
						+ "\t" + (0.05 / results.size())
						+ "\t" + (max.getPval() < (0.05 / results.size()))
				);
			}

		}

	}

	public double determineTopNumberOfVariantsWithinSignificantRegions(String[] significantRegionFiles, String[] assocFiles) throws IOException {
		BedFileReader reader = new BedFileReader();
		AssociationFile assocfilereader = new AssociationFile();
		int max = 0;
		for (int q = 0; q < significantRegionFiles.length; q++) {
			ArrayList<Feature> allRegions = reader.readAsList(significantRegionFiles[q]);
			for (Feature region : allRegions) {
				for (int i = 0; i < assocFiles.length; i++) {
					String assocFile = assocFiles[i];
					ArrayList<AssociationResult> results = assocfilereader.read(assocFile, region);
					if (results.size() > max) {
						max = results.size();
					}
				}
			}
		}
		System.out.println(max + " nr of tests within significant regions...");
		return 0.05 / max;

	}


	public void determineRegionSignificanceThresholds(String assocFile, String allRegionsFile, String diseaseSpecificRegionsFile, String out, double origThreshold, double significantThreshold) throws IOException {

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> allRegions = reader.readAsList(allRegionsFile);
		ArrayList<Feature> allICDiseaseRegions = reader.readAsList(diseaseSpecificRegionsFile);
		AssociationFile assocfilereader = new AssociationFile();
//		int nrTotal = 0;

//		for (int i = 0; i < allRegions.size(); i++) {
//			Feature region = allRegions.get(i);
//			ArrayList<AssociationResult> results = assocfilereader.read(assocFile, region);
//			nrTotal += results.size();
//		}


		TextFile outf = new TextFile(out, TextFile.W);

		for (int i = 0; i < allRegions.size(); i++) {
			Feature region = allRegions.get(i);
			boolean significantInIC = isSignificant(region, allICDiseaseRegions);
			ArrayList<AssociationResult> results = assocfilereader.read(assocFile, region);
			double threshold = origThreshold;
			if (significantInIC) {
				threshold = significantThreshold;
				outf.writeln(region.toString() + "\t" + results.size() + "\t" + threshold);
			} else {
				outf.writeln(region.toString() + "\t" + 0 + "\t" + threshold);
			}

		}

		outf.close();

	}

	private boolean isSignificant(Feature region, ArrayList<Feature> allICDiseaseRegions) {
		for (Feature f : allICDiseaseRegions) {
			if (region.overlaps(f)) {
				return true;
			}
		}
		return false;
	}
}
