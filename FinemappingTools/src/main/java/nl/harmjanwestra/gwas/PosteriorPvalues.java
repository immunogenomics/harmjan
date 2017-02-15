package nl.harmjanwestra.gwas;

import nl.harmjanwestra.gwas.CLI.PosteriorPvalueOptions;
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

		ArrayList<AssociationResult> assoc = associationFile.read(assocfile);
		System.out.println(assoc.size() + " associations loaded from " + assocfile);


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
		System.out.println("Credible sets are written here: " + output + "-credibleSets.txt");

		// write new output path.

		TextFile outf = new TextFile(output, TextFile.W);
		System.out.println("Writing output here: " + output);

		String header = associationFile.getHeader();
		header += "\tRegion\tBF\tPosterior";
		outf.writeln(header);
		for (int i = 0; i < assoc.size(); i++) {
			AssociationResult r = assoc.get(i);
			String assocStr = r.toString();
			if (r.getRegion() == null) {
				assocStr += "\t" + null;
				assocStr += "\t" + null;
				assocStr += "\t" + null;
			} else {
				assocStr += "\t" + r.getRegion().toString();
				assocStr += "\t" + r.getBf();
				assocStr += "\t" + r.getPosterior();
			}


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
		String[] csORs = new String[credibleSet.size()];
		String[] csNames = new String[credibleSet.size()];
		String line = "";
		for (int v = 0; v < credibleSet.size(); v++) {
			AssociationResult result = credibleSet.get(v);
			csPosteriors[v] = result.getPosterior();
			csPvals[v] = result.getPval();
			csORs[v] = Strings.concat(result.getORs(), Strings.colon);
			Feature f = result.getSnp();
			csNames[v] = f.getChromosome().toString() + ":" + f.getStart() + "-" + f.getName();
		}

		line += credibleSet.size();
		line += "\t" + Strings.concat(csNames, Strings.semicolon);
		line += "\t" + Strings.concat(csPvals, Strings.semicolon);
		line += "\t" + Strings.concat(csORs, Strings.semicolon);
		line += "\t" + Strings.concat(csPosteriors, Strings.semicolon);
		return line;
	}
}
