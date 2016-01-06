package nl.harmjanwestra.assoc;

import nl.harmjanwestra.assoc.CLI.PosteriorPvalueOptions;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 11/11/15.
 */
public class PosteriorPvalues {
	private final PosteriorPvalueOptions options;

//	public static void main(String[] args) {
//		PosteriorPvalues p = new PosteriorPvalues();
//		try {
////			p.mergeAssociationResults();
////			p.makeplotsOct27();
//
//
//			// calculate posteriors
//			String regionsFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";
//			String t1dDir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/Conditional/T1D/";
//			String[] refs = new String[]{"1kg", "seq", "1kg-seq-merged"};
//			String outfilename = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/Posteriors/T1D/";
//			for (int refId = 0; refId < refs.length; refId++) {
//
//				String out = outfilename + "/" + refs[refId] + "/";
//				Gpio.createDir(out);
//				String assocdir = t1dDir + "/" + refs[refId] + "/";
//				p.determinePosteriors(assocdir, regionsFile, out);
//			}
//
//
//			t1dDir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/Conditional/RA/";
//			outfilename = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/Posteriors/RA/";
//			for (int refId = 0; refId < refs.length; refId++) {
//
//				String out = outfilename + "/" + refs[refId] + "/";
//				Gpio.createDir(out);
//				String assocdir = t1dDir + "/" + refs[refId] + "/";
//				p.determinePosteriors(assocdir, regionsFile, out);
//			}
//
//
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//	}

	public PosteriorPvalues(PosteriorPvalueOptions options) throws IOException {
		this.options = options;
		determinePosteriors(options.getAssocFile(),
				options.getRegionFile(),
				options.getOutputPrefix(),
				options.getBayesThreshold());
	}

	public void determinePosteriors(String assocfile, String regionfile, String output, double bayesthreshold) throws IOException {
		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(regionfile);
		System.out.println(regions.size() + " regions in " + regionfile);

		AssociationFile associationFile = new AssociationFile();
		ApproximateBayesPosterior abp = new ApproximateBayesPosterior();

		ArrayList<AssociationResult> assoc = associationFile.loadConditionalAssocData(assocfile);


		// determine posteriors within defined regions
		TextFile crediblesetout = new TextFile(output + "-credibleSets.txt", TextFile.W);
		crediblesetout.writeln("Region\tN\tNames\tPvals\tORs\tPosteriors");
		for (Feature region : regions) {
			ArrayList<AssociationResult> regionResults = filter(region, assoc);
			abp.calculatePosterior(regionResults);
			ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(regionResults, bayesthreshold);
			crediblesetout.writeln(region.toString() + "\t" + getCredibleSetStr(credibleSet));
		}
		crediblesetout.close();

		// write new output file.

		TextFile outf = new TextFile(output, TextFile.W);

		String header = associationFile.getHeader();
		header += "\tRegion\tBF\tPosterior";
		outf.writeln(header);
		for (int i = 0; i < assoc.size(); i++) {
			AssociationResult r = assoc.get(i);
			String assocStr = r.toString();
			assocStr += "\t" + r.getRegion().toString()
					+ "\t" + r.getBf()
					+ "\t" + r.getPosterior();
			outf.writeln(assocStr);
		}
		outf.close();
	}

	private ArrayList<AssociationResult> filter(Feature region, ArrayList<AssociationResult> assoc) {
		ArrayList<AssociationResult> output = new ArrayList<>();
		for (AssociationResult r : assoc) {
			Feature snp = r.getSnp();
			if (snp.overlaps(region)) {
				r.setRegion(region);
				output.add(r);
			}
		}
		return output;
	}

	public String getCredibleSetStr(ArrayList<AssociationResult> credibleSet) {
		double[] csPosteriors = new double[credibleSet.size()];
		double[] csPvals = new double[credibleSet.size()];
		double[] csORs = new double[credibleSet.size()];
		String[] csNames = new String[credibleSet.size()];
		String line = "";
		for (int v = 0; v < credibleSet.size(); v++) {
			AssociationResult result = credibleSet.get(v);
			csPosteriors[v] = result.getPosterior();
			csPvals[v] = result.getPval();
			csORs[v] = result.getOr();
			Feature f = result.getSnp();
			csNames[v] = f.getChromosome().toString() + ":" + f.getStart() + "-" + f.getName();
		}

		line += "\t" + credibleSet.size();
		line += "\t" + Strings.concat(csNames, Strings.semicolon);
		line += "\t" + Strings.concat(csPvals, Strings.semicolon);
		line += "\t" + Strings.concat(csORs, Strings.semicolon);
		line += "\t" + Strings.concat(csPosteriors, Strings.semicolon);
		return line;
	}
}
