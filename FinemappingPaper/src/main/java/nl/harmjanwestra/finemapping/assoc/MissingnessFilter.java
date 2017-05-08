package nl.harmjanwestra.finemapping.assoc;

import nl.harmjanwestra.finemapping.annotation.RSIdRewriter;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by Harm-Jan on 05/06/17.
 */
public class MissingnessFilter {

	public static void main(String[] args) {


		MissingnessFilter m = new MissingnessFilter();
		double threshold = 7.5e-7;
		double bayesthreshold = 0.9;
		String tabixprefix = "D:\\Data\\Ref\\beagle1kg\\1kg.phase3.v5a.chrCHR.vcf.gz";
		String samplefilter = null;
		int nrheaderlines = 1;
		int[] columns = new int[]{3};
		int[] rsidcols = new int[]{2};
//		int[] combinedIdCols = new int[]{3};
//		int[] rsidcols = new int[]{2};

		String regionfile = "D:\\Cloud\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2016-09-06-SummaryStats\\NormalHWEP1e4\\AllICLoci-overlappingWithImmunobaseT1DOrRALoci-woMHC.bed";

		BedFileReader reader = new BedFileReader();


		String dir1 = "D:\\ wrongrsids\\";
		String dir2 = "D:\\Cloud\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-03-25-SummaryStats\\normal\\";
		String dir3 = "D:\\Cloud\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-03-25-SummaryStats\\variantstats\\";

		try {
			ArrayList<Feature> regions = reader.readAsList(regionfile);
			String assocfile = dir1 + "T1D-assoc0.3-COSMO-merged-posterior.txt.gz";
			String assocfileout = dir1 + "T1D-assoc0.3-COSMO-merged-missingpfiltered.txt.gz";
			String assocfileoutpost = dir1 + "T1D-assoc0.3-COSMO-merged-posterior-missingpfiltered.txt.gz";
			String assocfileoutpostrsidfix = dir2 + "T1D-assoc0.3-COSMO-merged-posterior.txt.gz";
			String missingfile = dir3 + "T1Dstats.txt";
			String out = dir3 + "T1DstatsComparison.txt";
			m.filter(missingfile, assocfile, assocfileout, threshold);
			m.determinePosteriors(assocfileout, regionfile, assocfileoutpost, bayesthreshold);

			RSIdRewriter rw = new RSIdRewriter();
			rw.run(assocfileoutpost, assocfileoutpostrsidfix, tabixprefix, samplefilter, nrheaderlines, columns, rsidcols, regions);

//			m.compare(missingfile, assocfile, out);

			assocfile = dir1 + "RA-assoc0.3-COSMO-merged-posterior.txt.gz";
			assocfileout = dir1 + "RA-assoc0.3-COSMO-merged-missingpfiltered.txt.gz";
			assocfileoutpost = dir1 + "RA-assoc0.3-COSMO-merged-posterior-missingpfiltered.txt.gz";
			assocfileoutpostrsidfix = dir2 + "RA-assoc0.3-COSMO-merged-posterior.txt.gz";
			missingfile = dir3 + "RAstats.txt";
			out = dir3 + "RAstatsComparison.txt";
			m.filter(missingfile, assocfile, assocfileout, threshold);
			m.determinePosteriors(assocfileout, regionfile, assocfileoutpost, bayesthreshold);
			rw.run(assocfileoutpost, assocfileoutpostrsidfix, tabixprefix, samplefilter, nrheaderlines, columns, rsidcols, regions);
//			m.compare(missingfile, assocfile, out);

			assocfile = dir1 + "META-assoc0.3-COSMO-merged-posterior.txt.gz";
			assocfileout = dir1 + "META-assoc0.3-COSMO-merged-missingpfiltered.txt.gz";
			assocfileoutpost = dir1 + "META-assoc0.3-COSMO-merged-posterior-missingpfiltered.txt.gz";
			assocfileoutpostrsidfix = dir2 + "META-assoc0.3-COSMO-merged-posterior.txt.gz";
			missingfile = dir3 + "METAstats.txt";
			out = dir3 + "METAstatsComparison.txt";
			m.filter(missingfile, assocfile, assocfileout, threshold);
			m.determinePosteriors(assocfileout, regionfile, assocfileoutpost, bayesthreshold);
			rw.run(assocfileoutpost, assocfileoutpostrsidfix, tabixprefix, samplefilter, nrheaderlines, columns, rsidcols, regions);
//			m.compare(missingfile, assocfile, out);

		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	public void filter(String missingnessfile, String assocfile, String out, double threshold) throws IOException {

		AssociationFile assocFile = new AssociationFile();
		ArrayList<AssociationResult> assocresults = assocFile.read(assocfile);
		String header = assocFile.getHeader();
		HashMap<String, Double> missingnessP = readMissingnessFile(missingnessfile);

		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln(header);
		int written = 0;
		int excluded = 0;

		for (int i = 0; i < assocresults.size(); i++) {

			AssociationResult r = assocresults.get(i);
			String id = r.getSnp().getChromosome().getNumber() + "_" + r.getSnp().getStart() + "_" + r.getSnp().getName();
			if (missingnessP.containsKey(id)) {
				double p = missingnessP.get(id);
				if (p > threshold) {
					r.getSnp().setMissingnessP(p);
					outf.writeln(r.toString());
					written++;
				} else {
					excluded++;
				}
			} else {
				excluded++;
			}
		}


		outf.close();
		System.out.println(written + " out of " + assocresults.size() + " written. " + excluded + " excluded.");

	}

	public void compare(String missingnessfile, String assocfile, String out) throws IOException {

		AssociationFile assocFile = new AssociationFile();
		ArrayList<AssociationResult> assocresults = assocFile.read(assocfile);
		String header = assocFile.getHeader();

		// read the missingnessfile
		HashMap<String, Double> missingnessP = readMissingnessFile(missingnessfile);


		// check whether we can map all variants to a missingnessp
		int found = 0;
		int notfound = 0;
		TextFile tf = new TextFile(out, TextFile.W);
		tf.writeln("id\tImpQual\tLog10P\tPosterior\tMissingnessP\tLog10MissingnessP");
		for (AssociationResult r : assocresults) {
//			String id = r.getSnp().toString();
			String id = r.getSnp().getChromosome().getNumber() + "_" + r.getSnp().getStart() + "_" + r.getSnp().getName();
			if (missingnessP.containsKey(id)) {
				found++;
				double p = missingnessP.get(id);
				double logp = -Math.log10(p);
				tf.writeln(id + "\t" + r.getSnp().getImputationQualityScore() + "\t" + r.getLog10Pval() + "\t" + r.getPosterior() + "\t" + p + "\t" + logp);
			} else {
				notfound++;
			}
		}
		tf.close();
		System.out.println("Found: " + found);
		System.out.println("Not found: " + notfound);


	}

	public HashMap<String, Double> readMissingnessFile(String file) throws IOException {
		TextFile tf = new TextFile(file, TextFile.R);
		String header = tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, Double> snpToMissingness = new HashMap<String, Double>();
		int ctr = 0;
		while (elems != null) {
			String snp = elems[3];
			Double p = Double.parseDouble(elems[6]);

			snpToMissingness.put(snp, p);
			elems = tf.readLineElems(TextFile.tab);
			ctr++;
//			if (ctr == 10) {
//				System.exit(-1);
//			} else {
//			System.out.println(snp + "\t" + p + "\t" + elems[6]);
//		}
		}

		tf.close();
		return snpToMissingness;

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
