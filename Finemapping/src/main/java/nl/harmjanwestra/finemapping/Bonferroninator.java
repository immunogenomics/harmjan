package nl.harmjanwestra.finemapping;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

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

			assocfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/RA-assoc0.3-COSMO-merged.txt.gz";
			diseasespecificregionfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseRALoci.bed";
			out = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/BonferroniThresholds/RA.txt";
			b.run(assocfile, allregionsfile, diseasespecificregionfile, out);

			// T1D
			assocfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/T1D-assoc0.3-COSMO-merged.txt.gz";
			diseasespecificregionfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DLoci.bed";
			out = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/BonferroniThresholds/T1D.txt";
			b.run(assocfile, allregionsfile, diseasespecificregionfile, out);

			// joint
			assocfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/META-assoc0.3-COSMO-merged-posterior.txt.gz";
			diseasespecificregionfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci-woMHC.txt";
			out = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/BonferroniThresholds/META.txt";
			b.run(assocfile, allregionsfile, diseasespecificregionfile, out);

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String assocFile, String allRegionsFile, String diseaseSpecificRegionsFile, String out) throws IOException {

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> allRegions = reader.readAsList(allRegionsFile);
		ArrayList<Feature> allICDiseaseRegions = reader.readAsList(diseaseSpecificRegionsFile);

		int nrTotal = 0;
		AssociationFile assocfilereader = new AssociationFile();
		for (int i = 0; i < allRegions.size(); i++) {
			Feature region = allRegions.get(i);
			ArrayList<AssociationResult> results = assocfilereader.read(assocFile, region);
			nrTotal += results.size();
		}


		TextFile outf = new TextFile(out, TextFile.W);

		for (int i = 0; i < allRegions.size(); i++) {
			Feature region = allRegions.get(i);
			boolean significantInIC = isSignificant(region, allICDiseaseRegions);
			ArrayList<AssociationResult> results = assocfilereader.read(assocFile, region);
			double threshold = 0.05 / nrTotal;
			if (significantInIC) {
				threshold = 0.05 / results.size();
				outf.writeln(region.toString() + "\t" + results.size() + "\t" + threshold);
			} else {
				outf.writeln(region.toString() + "\t" + nrTotal + "\t" + threshold);
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
