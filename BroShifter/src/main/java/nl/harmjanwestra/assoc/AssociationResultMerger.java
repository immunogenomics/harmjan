package nl.harmjanwestra.assoc;

import nl.harmjanwestra.assoc.CLI.AssociationResultMergerOptions;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 11/23/15.
 */
public class AssociationResultMerger {


//	public void mergeAssociationResults() throws IOException {
//
//		String sequencedRegionsFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/2015-03-30-allRegionsWith50PercentOfSamplesAbove20x.bed";
//		String assocDir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/";
//		String regionsFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";
//
//		String rsquareDir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/ImpQScores/";
//		String[] refs = new String[]{"1kg", "seq", "1kg-seq-merged", "ImmunoChipStudy"};
//		String gtf = "/Data/Ref/Annotation/UCSC/genes.gtf";
//		String outdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/";
//
//		String condiontalassocDir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/";
//		String conditionaloutdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/Conditional/Plots/";
//
//		String immunoChipT1D = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/T1D-Onengut/hg19_gwas_ic_t1d_onengut_cc_4_18_1.tab";
//		String immunoChipRA = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/RA-Eyre/hg19_gwas_ic_ra_eyre_4_18_0.tab";
//
//		int run = 0;
//
//		// make conditional plots
//		for (int d = 0; d < 2; d++) {
//			String ds = "T1D";
//			String ic = immunoChipT1D;
//			if (d == 1) {
//				ds = "RA";
//				ic = immunoChipRA;
//			}
//		}
//	}

	AssociationResultMergerOptions options;

	public AssociationResultMerger(AssociationResultMergerOptions options) throws IOException {
		this.options = options;

		if (options.isConcat()) {
			concat(options.getInput(), options.getOutputprefix());
		} else {
			BedFileReader reader = new BedFileReader();
			ArrayList<Feature> regions = reader.readAsList(options.getRegions());
			String[] refs = options.getRefStr().split(",");
			String[] input = options.getInput().split(",");
			mergeDatasetForDifferentReferences(options.getOutputprefix(),
					refs,
					input,
					regions,
					options.getBayesthreshold());

		}
	}

	private void concat(String fileStr, String outfile) throws IOException {
		String[] files = fileStr.split(",");
		TextFile out = new TextFile(outfile, TextFile.W);
		boolean headerwritten = false;
		for (String file : files) {
			if (Gpio.exists(file)) {
				TextFile in = new TextFile(file, TextFile.R);
				String ln = in.readLine();
				if (!headerwritten) {
					out.writeln(ln);
				} else {
					ln = in.readLine();
				}
				while (ln != null) {
					out.writeln(ln);
					ln = in.readLine();
				}
				in.close();
			}
		}
		out.close();
	}


	private void mergeDatasetForDifferentReferences(String outprefix,
													String[] refs,
													String[] refFiles,
													ArrayList<Feature> regions,
													double bayesthreshold) throws IOException {
		String outfilename = outprefix + "-MergedTable.txt";
		String credibleSetOut = outprefix + "-CredibleSets-" + bayesthreshold + ".txt";
		String header = "Region\tPosition";
		String credibleSetHeader = "Region";
		for (int ref = 0; ref < refs.length; ref++) {
			String rname = refs[ref];
			header += "\tRsName-" + rname
					+ "\tMAF-" + rname
					+ "\tBeta-" + rname
					+ "\tSE-" + rname
					+ "\tOR-" + rname
					+ "\tPVal-" + rname
					+ "\tBF-" + rname
					+ "\tPosterior-" + rname;
			credibleSetHeader += "\tSize-" + rname
					+ "\tTotalNrVariants-" + rname
					+ "\tVariants-" + rname
					+ "\tPvals-" + rname
					+ "\tORs-" + rname
					+ "\tPosteriors-" + rname;
		}

		TextFile outfile = new TextFile(outfilename, TextFile.W);
		TextFile outfilecs = new TextFile(credibleSetOut, TextFile.W);
		outfile.writeln(header);
		outfilecs.writeln(credibleSetHeader);

		for (int r = 0; r < regions.size(); r++) {
			Feature region = regions.get(r);

			AssociationFile associationFile = new AssociationFile();

			ArrayList<ArrayList<AssociationResult>> allResults = new ArrayList<>();
			ArrayList<ArrayList<AssociationResult>> credibleSetsPerDataset = new ArrayList<>();

			HashSet<Integer> uniquePositions = new HashSet<Integer>();
			ArrayList<Integer> allPositions = new ArrayList<Integer>();

			HashMap<Integer, Integer> posToMaxCount = new HashMap<Integer, Integer>();

			int[] nrAssociationResultsPerDataset = new int[refs.length];
			Chromosome chr = region.getChromosome();
			ApproximateBayesPosterior abp = new ApproximateBayesPosterior();

			if (!chr.equals(Chromosome.X)) {
				for (int refId = 0; refId < refs.length; refId++) {
					// read assoc file
					String ref = refs[refId];


					String regionStr = region.toString();

					String assocFile = refFiles[refId];//indir + "/Conditional/" + ds + "/" + ref + "/" + chr.toString() + "-" + regionStr + "-assoc-" + iter + ".txt";
					ArrayList<AssociationResult> associationResults = null;
					if (assocFile.endsWith(".tab")) {
						associationResults = associationFile.readVariantPValues(assocFile, region);
					} else {
						associationResults = associationFile.loadConditionalAssocData(assocFile, region);
					}

					nrAssociationResultsPerDataset[refId] = associationResults.size();

					ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(associationResults, bayesthreshold);
					credibleSetsPerDataset.add(credibleSet);

					HashMap<Integer, Integer> posToMaxCountDs = new HashMap<Integer, Integer>();

					for (AssociationResult q : associationResults) {
						int pos = q.getSnp().getStart();
						if (!uniquePositions.contains(pos)) {
							uniquePositions.add(pos);
							allPositions.add(pos);
						}
						Integer ct = posToMaxCountDs.get(pos);
						if (ct == null) {
							ct = 1;
						} else {
							ct++;
						}

						posToMaxCountDs.put(pos, ct);
					}
					allResults.add(associationResults);

					for (Integer pos : uniquePositions) {
						Integer ct = posToMaxCountDs.get(pos);
						Integer otherct = posToMaxCount.get(pos);
						if (ct != null) {
							if (otherct != null) {
								if (ct > otherct) {
									posToMaxCount.put(pos, ct);
								}
							} else {
								posToMaxCount.put(pos, ct);
							}
						}
					}
				}

				Collections.sort(allPositions);

				// iterate all positions
				for (int p = 0; p < allPositions.size(); p++) {
					int pos = allPositions.get(p);
					ArrayList<ArrayList<AssociationResult>> resultsForPos = new ArrayList<>();
					int maxNr = 0;

					// get all results for this position
					for (int q = 0; q < allResults.size(); q++) {
						ArrayList<AssociationResult> resultForDsForPos = new ArrayList<>();
						ArrayList<AssociationResult> results = allResults.get(q);

						for (AssociationResult result : results) {
							int pos2 = result.getSnp().getStart();
							if (pos2 == pos) {
								// add result
								resultForDsForPos.add(result);
							}
						}

						resultsForPos.add(resultForDsForPos);
						if (resultForDsForPos.size() > maxNr) {
							maxNr = resultForDsForPos.size();
						}
					}

					// get first result
					for (int q = 0; q < maxNr; q++) {
						String line = region.toString() + "\t" + pos;


						for (int s = 0; s < resultsForPos.size(); s++) {
							ArrayList<AssociationResult> resultForDs = resultsForPos.get(s);
							if (resultForDs.size() > q) {
								AssociationResult result = resultForDs.get(q);

								line += "\t" + result.getSnp().getName()
										+ "\t" + result.getMaf()
										+ "\t" + result.getBeta()
										+ "\t" + result.getSe()
										+ "\t" + result.getOr()
										+ "\t" + result.getPval()
										+ "\t" + result.getBf()
										+ "\t" + result.getPosterior();

							} else {
								// add some nulls
								line += "\t" + null
										+ "\t" + null
										+ "\t" + null
										+ "\t" + null
										+ "\t" + null
										+ "\t" + null
										+ "\t" + null
										+ "\t" + null;

							}
						}


						// write line
						outfile.writeln(line);
					}
				}
			}

			String line = region.toString();
			for (int dataset = 0; dataset < credibleSetsPerDataset.size(); dataset++) {
				ArrayList<AssociationResult> credibleSet = credibleSetsPerDataset.get(dataset);
				double[] csPosteriors = new double[credibleSet.size()];
				double[] csPvals = new double[credibleSet.size()];
				double[] csORs = new double[credibleSet.size()];
				String[] csNames = new String[credibleSet.size()];
				for (int v = 0; v < credibleSet.size(); v++) {
					AssociationResult result = credibleSet.get(v);
					csPosteriors[v] = result.getPosterior();
					csPvals[v] = result.getPval();
					csORs[v] = result.getOr();
					Feature f = result.getSnp();
					csNames[v] = f.getChromosome().toString() + ":" + f.getStart() + "-" + f.getName();
				}

				line += "\t" + credibleSet.size();
				line += "\t" + nrAssociationResultsPerDataset[dataset];
				line += "\t" + Strings.concat(csNames, Strings.semicolon);
				line += "\t" + Strings.concat(csPvals, Strings.semicolon);
				line += "\t" + Strings.concat(csORs, Strings.semicolon);
				line += "\t" + Strings.concat(csPosteriors, Strings.semicolon);
			}
			outfilecs.writeln(line);
		}
		outfile.close();
		outfilecs.close();
	}
}
