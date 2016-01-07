package nl.harmjanwestra.assoc;

import nl.harmjanwestra.assoc.CLI.AssociationResultMergerOptions;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

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

		System.out.println("Concatenating " + files.length + " files into: " + outfile);
		int totalvars = 0;
		for (String file : files) {
			if (Gpio.exists(file)) {
				System.out.println("Concatenating file: " + file);
				AssociationFile assocFile = new AssociationFile();
				ArrayList<AssociationResult> result = assocFile.read(file);
				System.out.println(result.size() + " associations in file..");
				totalvars += result.size();
				if (!headerwritten) {
					out.writeln(assocFile.getHeader());
					System.out.println("Writing header...");
					headerwritten = true;
				}
				for (AssociationResult r : result) {
					out.writeln(r.toString());
				}
			} else {
				System.out.println("Could not find file: " + file);
			}
		}
		out.close();
		System.out.println("Done. " + totalvars + " associations written.");
	}

	private void mergeDatasetDifferentReferences(String outputprefix,
												 String refStr,
												 String refFileStr,
												 String regionfile,
												 double bayesthreshold) throws IOException {

		String[] refs = refStr.split(",");
		String[] refFiles = refFileStr.split(",");
		ArrayList<ArrayList<AssociationResult>> associationResults = new ArrayList<>();
		AssociationFile assocFile = new AssociationFile();
		for (int f = 0; f < refFiles.length; f++) {
			ArrayList<AssociationResult> results = assocFile.read(refFiles[f]);
			System.out.println(results.size() + " associations read from: " + refs[f] + "; " + refFiles[f]);
			associationResults.add(results);
		}

		// get a list of unique variants
		HashMap<Feature, Integer> uniqueVariants = new HashMap<Feature, Integer>();
		ArrayList<Feature> allVariants = new ArrayList<>();
		for (int q = 0; q < associationResults.size(); q++) {
			ArrayList<AssociationResult> list = associationResults.get(q);
			for (AssociationResult r : list) {
				Feature snp = r.getSnp();
				if (!uniqueVariants.containsKey(snp)) {
					uniqueVariants.put(snp, uniqueVariants.size());
					allVariants.add(snp);
				}
			}
			System.out.println(uniqueVariants.size() + " unique variants total after " + refs[q]);
		}


		System.out.println(uniqueVariants.size() + " unique variants total");

		// put them in a matrix for ease of use
		AssociationResult[][] matrix = new AssociationResult[refs.length][uniqueVariants.size()];
		for (int q = 0; q < associationResults.size(); q++) {
			ArrayList<AssociationResult> list = associationResults.get(q);
			for (AssociationResult r : list) {
				Feature snp = r.getSnp();
				Integer index = uniqueVariants.get(snp);
				if (index != null) {
					matrix[q][index] = r;
				}
			}
			System.out.println(uniqueVariants.size() + " unique variants total after " + refs[q]);
		}

		// assign regions to variants
		BedFileReader bf = new BedFileReader();
		ArrayList<Feature> regions = bf.readAsList(regionfile);
		Collections.sort(regions, new FeatureComparator(false));
		Collections.sort(allVariants, new FeatureComparator(false));
		for (Feature f : allVariants) {
			for (Feature f2 : regions) {
				if (f.overlaps(f2)) {
					f.setParent(f2);
					break;
				}
			}
		}


		// headers etc
		String outfilename = outputprefix + "-MergedTable.txt";
		String credibleSetOut = outputprefix + "-CredibleSets-" + bayesthreshold + ".txt";
		String header = "Region\tPosition\tRsName";
		String credibleSetHeader = "Region\tTotalVariantsInRegion";
		for (int ref = 0; ref < refs.length; ref++) {
			String rname = refs[ref];
			header += "\tMAF-" + rname
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
		outfile.writeln(header);

		for (int f = 0; f < allVariants.size(); f++) {
			Feature snp = allVariants.get(f);
			Integer index = uniqueVariants.get(snp);

			StringBuilder ln = new StringBuilder();
			ln.append(snp.getParent().toString());
			ln.append("\t").append(snp.getStart());
			ln.append("\t").append(snp.getName());
			for (int ref = 0; ref < refs.length; ref++) {
				AssociationResult r = matrix[ref][index];
				if (r != null) {


					ln.append("\t").append(r.getMaf());
					ln.append("\t").append(Strings.concat(r.getBeta(), Strings.semicolon));
					ln.append("\t").append(Strings.concat(r.getSe(), Strings.semicolon));
					ln.append("\t").append(Strings.concat(r.getORs(), Strings.semicolon));
					ln.append("\t").append(r.getPval());
					ln.append("\t").append(r.getBf());
					ln.append("\t").append(r.getPosterior());
				} else {
					ln.append("\tnull");
					ln.append("\tnull");
					ln.append("\tnull");
					ln.append("\tnull");
					ln.append("\tnull");
					ln.append("\tnull");
					ln.append("\tnull");
				}
			}
			outfile.writeln(ln.toString());
		}
		outfile.close();


		// determine credible sets per region....
		ApproximateBayesPosterior abp = new ApproximateBayesPosterior();
		TextFile outfilecs = new TextFile(credibleSetOut, TextFile.W);
		outfilecs.writeln(credibleSetHeader);


		for (int r = 0; r < regions.size(); r++) {
			// get variants in region
			Feature region = regions.get(r);
			ArrayList<Integer> variantsInRegion = new ArrayList<>();
			for (Feature f : allVariants) {
				if (f.overlaps(region)) {
					Integer index = uniqueVariants.get(f);
					variantsInRegion.add(index);
				}
			}
			String line = region.toString() + "\t" + variantsInRegion.size();


			ArrayList<ArrayList<AssociationResult>> variantsPerDataset = new ArrayList<>();
			for (int ref = 0; ref < refs.length; ref++) {
				ArrayList<AssociationResult> dsAssociations = new ArrayList<>();
				for (Integer i : variantsInRegion) {
					AssociationResult f = matrix[ref][i];
					if (f != null) {
						dsAssociations.add(f);
					}
				}

				ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(dsAssociations, bayesthreshold);

				double[] csPosteriors = new double[credibleSet.size()];
				double[] csPvals = new double[credibleSet.size()];
				String[] csORs = new String[credibleSet.size()];
				String[] csNames = new String[credibleSet.size()];
				for (int v = 0; v < credibleSet.size(); v++) {
					AssociationResult result = credibleSet.get(v);
					csPosteriors[v] = result.getPosterior();
					csPvals[v] = result.getPval();
					csORs[v] = Strings.concat(result.getORs(), Strings.colon);
					Feature f = result.getSnp();
					csNames[v] = f.getChromosome().toString() + ":" + f.getStart() + "-" + f.getName();
				}

				line += "\t" + credibleSet.size();
				line += "\t" + dsAssociations.size();
				line += "\t" + Strings.concat(csNames, Strings.semicolon);
				line += "\t" + Strings.concat(csPvals, Strings.semicolon);
				line += "\t" + Strings.concat(csORs, Strings.semicolon);
				line += "\t" + Strings.concat(csPosteriors, Strings.semicolon);
			}
			outfilecs.writeln(line);
		}
		outfilecs.close();
	}

}
