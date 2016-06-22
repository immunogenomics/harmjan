package nl.harmjanwestra.finemapping;

import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by hwestra on 6/21/16.
 */
public class LocusMatcher {


	public static void main(String[] args) {

		String icregionfile = "";
		String imbasefile = "";
		String outbed = "";


		LocusMatcher m = new LocusMatcher();
		try {

			icregionfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci.bed";
			imbasefile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ImmunoBase/Hs_GRCh37-RA-assoc_tableBED";
			outbed = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseRALoci.bed";
			m.filter(icregionfile, imbasefile, outbed, false);

			imbasefile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseRALoci.bed";
			icregionfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-19-SummaryStats/RA-assoc0.3-COSMO-significantregions-1e5.bed";
			outbed = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseRALoci-overlappingWithSignificantFinemapLoci.bed";
			m.filter(icregionfile, imbasefile, outbed, true);

			icregionfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci.bed";
			imbasefile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ImmunoBase/Hs_GRCh37-T1D-assoc_tableBED";
			outbed = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DLoci.bed";
			m.filter(icregionfile, imbasefile, outbed, false);

			imbasefile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DLoci.bed";
			icregionfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-19-SummaryStats/T1D-assoc0.3-COSMO-significantregions-1e5.bed";
			outbed = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DLoci-overlappingWithSignificantFinemapLoci.bed";
			m.filter(icregionfile, imbasefile, outbed, true);

			//
			String bed1 = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseRALoci.bed";
			String bed2 = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DLoci.bed";
			String outb = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
			m.merge(bed1, bed2, outb);

		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void filter(String icregionfile, String imbasefile, String outbed, boolean printnonoverlap) throws IOException {

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regionsImmunobase = reader.readAsList(imbasefile);
		System.out.println(regionsImmunobase.size() + " regions in file: " + imbasefile);
		ArrayList<Feature> regionsIC = reader.readAsList(icregionfile);
		System.out.println(regionsIC.size() + " regions in file: " + icregionfile);

		TextFile out = new TextFile(outbed, TextFile.W);
		int nroverlap = 0;
		for (int i = 0; i < regionsIC.size(); i++) {
			Feature f = regionsIC.get(i);
			boolean overlap = false;
			for (Feature f2 : regionsImmunobase) {
				if (f.overlaps(f2)) {
					overlap = true;
				}
			}
			if (overlap) {
				out.writeln(f.toBedString());
				nroverlap++;
			} else if (printnonoverlap) {
				System.out.println("non overlapping: " + f.toBedString());
			}
		}

		System.out.println(nroverlap + " overlap");
		out.close();

	}


	public void merge(String bed1, String bed2, String bedout) throws IOException {

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> r1 = reader.readAsList(bed1);
		ArrayList<Feature> r2 = reader.readAsList(bed2);

		HashSet<Feature> f = new HashSet<Feature>();
		f.addAll(r1);
		f.addAll(r2);

		TextFile out = new TextFile(bedout, TextFile.W);
		for (Feature a : f) {
			out.writeln(a.toBedString());
		}
		out.close();

	}

}
