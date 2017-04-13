package nl.harmjanwestra.finemapping.genotypesandimputation;

import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

/**
 * Created by hwestra on 6/21/16.
 */
public class SupTable2LocusMatcher {


	public static void main(String[] args) {

		String icregionfile = "";
		String imbasefile = "";
		String outbed = "";


		SupTable2LocusMatcher m = new SupTable2LocusMatcher();
		try {

			icregionfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci.bed";
//			imbasefile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ImmunoBase/Hs_GRCh37-RA-assoc_tableBED";
//			outbed = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseRALoci.bed";
//			m.filter(icregionfile, imbasefile, outbed, false);
//
//			imbasefile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseRALoci.bed";
//			icregionfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-19-SummaryStats/RA-assoc0.3-COSMO-significantregions-1e5.bed";
//			outbed = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseRALoci-overlappingWithSignificantFinemapLoci.bed";
//			m.filter(icregionfile, imbasefile, outbed, true);
//
//			icregionfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci.bed";
//			imbasefile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ImmunoBase/Hs_GRCh37-T1D-assoc_tableBED";
//			outbed = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DLoci.bed";
//			m.filter(icregionfile, imbasefile, outbed, false);
//
//			imbasefile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DLoci.bed";
//			icregionfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-19-SummaryStats/T1D-assoc0.3-COSMO-significantregions-1e5.bed";
//			outbed = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DLoci-overlappingWithSignificantFinemapLoci.bed";
//			m.filter(icregionfile, imbasefile, outbed, true);
//
//			//
//			String bed1 = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseRALoci.bed";
//			String bed2 = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DLoci.bed";
//			String outb = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
//			m.merge(bed1, bed2, outb);

			String out = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.txt";
			m.combine(icregionfile,
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ImmunoBase/Hs_GRCh37-RA-assoc_tableBED",
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ImmunoBase/Hs_GRCh37-T1D-assoc_tableBED",
					out);


		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private void combine(String icregionfile, String bed1, String bed2, String out) throws IOException {

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> icregions = reader.readAsList(icregionfile);
		ArrayList<Feature> bed1regions = reader.readAsList(bed1);
		Collections.sort(bed1regions, new FeatureComparator());
		ArrayList<Feature> bed2regions = reader.readAsList(bed2);
		Collections.sort(bed2regions, new FeatureComparator());

		TextFile outf = new TextFile(out, TextFile.W);

		outf.writeln("ImmunoChip region\t\t\tRA associated region\t\t\tT1D associated region");
		outf.writeln("Chr\tStart\tEnd\tChr\tStart\tEnd\tChr\tStart\tEnd");

		for (int f = 0; f < icregions.size(); f++) {

			Feature region0 = icregions.get(f);
			ArrayList<Feature> overlap1 = new ArrayList<>();
			ArrayList<Feature> overlap2 = new ArrayList<>();
			for (int i = 0; i < bed1regions.size(); i++) {
				Feature region1 = bed1regions.get(i);
				if (region1.overlaps(region0)) {
					overlap1.add(region1);
				}
			}
			for (int i = 0; i < bed2regions.size(); i++) {
				Feature region1 = bed2regions.get(i);
				if (region1.overlaps(region0)) {
					overlap2.add(region1);
				}
			}

			if (overlap1.isEmpty() && overlap2.isEmpty()) {
			} else {
				String ln = region0.getChromosome().getName() + "\t" + region0.getStart() + "\t" + region0.getStop();

				int max = 0;
				if (overlap1.size() > max) {
					max = overlap1.size();
				}
				if (overlap2.size() > max) {
					max = overlap2.size();
				}

				String[] overlap1str = new String[max];
				String[] overlap2str = new String[max];


				if (!overlap1.isEmpty()) {

					for (int i = 0; i < overlap1.size(); i++) {
						Feature region1 = overlap1.get(i);
						overlap1str[i] = region1.getChromosome().getName() + "\t" + region1.getStart() + "\t" + region1.getStop();
					}
				}

				if (!overlap2.isEmpty()) {
					for (int i = 0; i < overlap2.size(); i++) {
						Feature region1 = overlap2.get(i);
						overlap2str[i] = region1.getChromosome().getName() + "\t" + region1.getStart() + "\t" + region1.getStop();
					}
				}

				for (int i = 0; i < max; i++) {
					String str1 = overlap1str[i];
					String str2 = overlap2str[i];
					if (str1 == null) {
						str1 = "\t\t";
					}
					if (str2 == null) {
						str2 = "\t\t";
					}

					if (i == 0) {
						outf.writeln(ln + "\t" + str1 + "\t" + str2);
					} else {
						outf.writeln("\t\t\t" + str1 + "\t" + str2);
					}
				}


			}
		}
		outf.close();
	}

	public void filter(String icregionfile, String imbasefile, String outbed, boolean printnonoverlap) throws IOException {

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regionsImmunobase = reader.readAsList(imbasefile);
		System.out.println(regionsImmunobase.size() + " regions in path: " + imbasefile);
		ArrayList<Feature> regionsIC = reader.readAsList(icregionfile);
		System.out.println(regionsIC.size() + " regions in path: " + icregionfile);

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
