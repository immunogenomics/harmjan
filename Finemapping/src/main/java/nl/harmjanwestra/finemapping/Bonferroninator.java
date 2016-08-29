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

			String[] significantRegionFiles = new String[]{
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/RA-assoc0.3-COSMO-significantregions-75e7.bed",
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/T1D-assoc0.3-COSMO-significantregions-75e7.bed",
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/META-assoc0.3-COSMO-significantregions-75e7.bed",
			};
			String[] assocFiles = new String[]{
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/RA-assoc0.3-COSMO-merged.txt.gz",
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/T1D-assoc0.3-COSMO-merged.txt.gz",
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/META-assoc0.3-COSMO-merged-posterior.txt.gz"
			};

			double defaultthreshold = 7.5E-7;
			double significantthreshold = b.determineTopNumberOfVariantsWithinSignificantRegions(significantRegionFiles, assocFiles);

			System.out.println();
			System.out.println("---");
			System.out.println();

			assocfile = assocFiles[0];
			diseasespecificregionfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseRALoci.bed";
			out = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/BonferroniThresholds/RA.txt";
			b.run(assocfile, allregionsfile, diseasespecificregionfile, out, defaultthreshold, significantthreshold);

			// T1D
			assocfile = assocFiles[1];
			diseasespecificregionfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DLoci.bed";
			out = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/BonferroniThresholds/T1D.txt";
			b.run(assocfile, allregionsfile, diseasespecificregionfile, out, defaultthreshold, significantthreshold);

			// joint
			assocfile = assocFiles[2];
			diseasespecificregionfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci-woMHC.txt";
			out = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/BonferroniThresholds/META.txt";
			b.run(assocfile, allregionsfile, diseasespecificregionfile, out, defaultthreshold, significantthreshold);


		} catch (IOException e) {
			e.printStackTrace();
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


	public void run(String assocFile, String allRegionsFile, String diseaseSpecificRegionsFile, String out, double origThreshold, double significantThreshold) throws IOException {

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
