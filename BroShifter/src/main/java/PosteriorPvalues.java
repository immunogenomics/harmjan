import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 11/11/15.
 */
public class PosteriorPvalues {

	public static void main(String[] args) {
		PosteriorPvalues p = new PosteriorPvalues();
		try {
//			p.mergeAssociationResults();
//			p.makeplotsOct27();


			// calculate posteriors
			String regionsFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";
			String t1dDir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/Conditional/T1D/";
			String[] refs = new String[]{"1kg", "seq", "1kg-seq-merged"};
			String outfilename = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/Posteriors/T1D/";
			for (int refId = 0; refId < refs.length; refId++) {

				String out = outfilename + "/" + refs[refId] + "/";
				Gpio.createDir(out);
				String assocdir = t1dDir + "/" + refs[refId] + "/";
				p.determinePosteriors(assocdir, regionsFile, out);
			}


			t1dDir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/Conditional/RA/";
			outfilename = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/Posteriors/RA/";
			for (int refId = 0; refId < refs.length; refId++) {

				String out = outfilename + "/" + refs[refId] + "/";
				Gpio.createDir(out);
				String assocdir = t1dDir + "/" + refs[refId] + "/";
				p.determinePosteriors(assocdir, regionsFile, out);
			}


		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void determinePosteriors(String assocdir, String regionfile, String output) throws IOException {
		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(regionfile);
		AssociationFile associationFile = new AssociationFile();
		ApproximateBayesPosterior abp = new ApproximateBayesPosterior();
		for (Feature region : regions) {


			String regionStr = region.toString();

			String assocfile = assocdir + region.getChromosome().toString() + "-" + regionStr + "-assoc-" + 0 + ".txt";

			if (Gpio.exists(assocfile)) {
				String regionName = region.getChromosome().toString() + "_" + region.getStart() + "-" + region.getStop();
				TextFile out = new TextFile(output + regionName + ".txt", TextFile.W);
				String header = "Chr\tPos\tId\tCombinedId\tBeta\tSe\tPVal\tPosterior";
				out.writeln(header);

				ArrayList<AssociationResult> assoc = associationFile.loadConditionalAssocData(assocfile, region);
				abp.calculateABF(assoc);

				for (AssociationResult r : assoc) {
//					if (r.getPval() > 0) {
					double abf = r.getAbf();
					if (!Double.isNaN(abf)) {

						Feature snp = r.getSnp();
						String ln = snp.getChromosome().toString()
								+ "\t" + snp.getStart()
								+ "\t" + snp.getName()
								+ "\t" + snp.getChromosome().toString() + ":" + snp.getStart() + "-" + snp.getName()
								+ "\t" + r.getBeta()
								+ "\t" + r.getSe()
								+ "\t" + r.getPval()
								+ "\t" + r.getAbf();
						out.writeln(ln);
					}
//					}
				}

				out.close();
			} else {
				System.out.println("Could not find file: " + assocfile);
			}

		}

	}









}
