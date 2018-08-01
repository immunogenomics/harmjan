package nl.harmjanwestra.finemapping.assoc;

import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by Harm-Jan on 04/21/16.
 */
public class DetermineBestAssocPerLocus {

	public static void main(String[] args) {
		DetermineBestAssocPerLocus d = new DetermineBestAssocPerLocus();
		try {
//			String t1d = "D:\\upload\\ic\\ImmunoChip\\hg19_gwas_ic_t1d_onengut_cc_4_19_1.tab";
//			String ra = "D:\\upload\\ic\\ImmunoChip\\hg19_gwas_ic_ra_eyre_4_19_1.tab";
//			String bedfile = "D:\\tmp\\2016-04-21\\AllICLoci.bed";
//			String outt1d = "D:\\tmp\\2016-04-21\\topt1deffects.txt";
//			String outra = "D:\\tmp\\2016-04-21\\topraeffects.txt";

			String t1d = "/Data/ImmunoChip/ImmunoBase/ImmunoChip/hg19_gwas_ic_t1d_onengut_cc_4_19_1.tab";
			String ra = "/Data/ImmunoChip/ImmunoBase/ImmunoChip/hg19_gwas_ic_ra_eyre_4_19_1.tab";
			String bedfile = "/Data/tmp/2016-04-21/AllICLoci.bed";
			String outt1d = "/Data/tmp/2016-04-21/onengut-topt1deffects.txt";
			String outra = "/Data/tmp/2016-04-21/eyre-topraeffects.txt";

			d.getSignificantLociForICStudies(bedfile, t1d, outt1d, true);

			d.getSignificantLociForICStudies(bedfile, ra, outra, true);

			d.merge("/Data/tmp/2016-04-21/RA-assoc0.3-eur-merged-topvariantperlocus.bed-topAssociations.txt",
					outra,
					"/Data/tmp/2016-04-21/RA-assoc0.3-eur-merged-topvariantperlocus.bed-topAssociations-snppairs.txt");

			d.merge("/Data/tmp/2016-04-21/RA-assoc0.3-merged-topvariantperlocus.bed-topAssociations.txt",
					outra,
					"/Data/tmp/2016-04-21/RA-assoc0.3-merged-topvariantperlocus.bed-topAssociations-snppairs.txt");

			d.merge("/Data/tmp/2016-04-21/T1D-assoc0.3-eur-merged-topvariantperlocus.bed-topAssociations.txt",
					outt1d,
					"/Data/tmp/2016-04-21/T1D-assoc0.3-eur-merged-topvariantperlocus.bed-topAssociations-snppairs.txt");

			d.merge("/Data/tmp/2016-04-21/T1D-assoc0.3-merged-topvariantperlocus.bed-topAssociations.txt",
					outt1d,
					"/Data/tmp/2016-04-21/T1D-assoc0.3-merged-topvariantperlocus.bed-topAssociations-snppairs.txt");


//			d.merge("D:\\tmp\\2016-04-21\\T1D-assoc0.3-merged-significantloci-3.23E-7.bed-topAssociations.txt",
//					outt1d,
//					"D:\\tmp\\2016-04-21\\T1D-assoc0.3-merged-significantloci-3.23E-7.bed-topAssociations-snppairs.txt");
//
//			d.merge("D:\\tmp\\2016-04-21\\T1D-assoc0.3-eur-merged-significantloci-3.23E-7.bed-topAssociations.txt",
//					outt1d,
//					"D:\\tmp\\2016-04-21\\T1D-assoc0.3-eur-merged-significantloci-3.23E-7.bed-topAssociations-snppairs.txt");
//
//			d.merge("D:\\tmp\\2016-04-21\\RA-assoc0.3-merged-significantloci-3.23E-7.bed-topAssociations.txt",
//					outra,
//					"D:\\tmp\\2016-04-21\\RA-assoc0.3-merged-significantloci-3.23E-7.bed-topAssociations-snppairs.txt");
//
//			d.merge("D:\\tmp\\2016-04-21\\RA-assoc0.3-eur-merged-significantloci-3.23E-7.bed-topAssociations.txt",
//					outra,
//					"D:\\tmp\\2016-04-21\\RA-assoc0.3-eur-merged-significantloci-3.23E-7.bed-topAssociations-snppairs.txt");

//			String file1 = "/Data/tmp/2016-04-21/RA-assoc0.3-merged-significantloci-3.23E-7.bed-topAssociations.txt";
//			String file2 = "/Data/tmp/2016-04-21/T1D-assoc0.3-merged-significantloci-3.23E-7.bed-topAssociations.txt";
//			String out = "/Data/tmp/2016-04-21/overlappingregions-cosmo.txt";
//			d.overlapRegions(file1, file2, out);

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void getSignificantLociForICStudies(String bedfile, String tabfile, String out, boolean onlysignificantloci) throws IOException {

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bedfile);

		TextFile tf = new TextFile(tabfile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<Feature, Pair<Feature, Double>> bestAssocPerRegion = new HashMap<>();
		while (elems != null) {
			// rs61733845	1	1118275	0.3552	1.0406	\N	\N	G>A
			String snp = elems[0];
			Chromosome chr = Chromosome.parseChr(elems[1]);
			Integer pos = Integer.parseInt(elems[2]);
			Double pval = -Math.log10(Double.parseDouble(elems[3]));
			Feature f = new Feature();
			f.setChromosome(chr);
			f.setStart(pos);
			f.setStop(pos);
			f.setName(snp);

			boolean overlap = false;
			Feature regionOfOverlap = null;
			for (Feature region : regions) {
				if (region.overlaps(f)) {
					overlap = true;
					regionOfOverlap = region;
					break;
				}
			}

			if (overlap) {
				Pair<Feature, Double> d = bestAssocPerRegion.get(regionOfOverlap);
				if (d == null) {
					Pair<Feature, Double> p = new Pair<Feature, Double>(f, pval);
					bestAssocPerRegion.put(regionOfOverlap, p);
				} else {
					double bestpval = d.getRight();

					if (pval > bestpval) {
						Pair<Feature, Double> p = new Pair<Feature, Double>(f, pval);
						bestAssocPerRegion.put(regionOfOverlap, p);
					}
				}
			}

			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();

		TextFile outf = new TextFile(out, TextFile.W);

		double t = -Math.log10(3.23E-7);
		for (Feature f : regions) {
			Pair<Feature, Double> d = bestAssocPerRegion.get(f);
			if (d != null) {
				Feature snp = d.getLeft();

				if (!onlysignificantloci || (d.getRight() > t)) {
					outf.writeln(f.toString() + "\t" + snp.getName() + "\t" + snp.getChromosome().toString() + "\t" + snp.getStart() + "\t" + d.getRight() + "\t" + (d.getRight() > t));
				}
			}
		}
		outf.close();
	}

	public void merge(String assocfile, String icfile, String out) throws IOException {

		TextFile in1 = new TextFile(assocfile, TextFile.R);
		String[] elems = in1.readLineElems(TextFile.tab);
		HashMap<String, String> regionToSNP = new HashMap<>();
		HashMap<String, String> regionToP = new HashMap<>();
		while (elems != null) {

			String region = elems[0];
			String snp = elems[1] + "_" + elems[2] + "_" + elems[3];
			regionToSNP.put(region, snp);
			regionToP.put(region, elems[elems.length - 1]);

			elems = in1.readLineElems(TextFile.tab);

		}
		in1.close();

		TextFile outf = new TextFile(out, TextFile.W);

		TextFile in2 = new TextFile(icfile, TextFile.R);
		elems = in2.readLineElems(TextFile.tab);
		while (elems != null) {
			String region = elems[0];
			String snp = elems[2] + "_" + elems[3] + "_" + elems[1];
			if (regionToSNP.containsKey(region)) {
				outf.writeln(regionToSNP.get(region) + "\t" + snp + "\t" + regionToP.get(region) + "\t" + elems[4]);
			} else {
				System.out.println(region.toString() + " not found?");
			}
			elems = in2.readLineElems(TextFile.tab);

		}
		in1.close();
		outf.close();


	}

	public void overlapRegions(String file1, String file2, String out) throws IOException {
		TextFile in1 = new TextFile(file1, TextFile.R);
		String[] elems = in1.readLineElems(TextFile.tab);
		HashMap<String, String> regionToSNP = new HashMap<>();
		HashMap<String, String> regionToP = new HashMap<>();
		while (elems != null) {

			String region = elems[0];
			String snp = elems[1] + "_" + elems[2] + "_" + elems[3];
			regionToSNP.put(region, snp);
			regionToP.put(region, elems[elems.length - 1]);

			elems = in1.readLineElems(TextFile.tab);

		}
		in1.close();


		TextFile outf = new TextFile(out, TextFile.W);
		TextFile in2 = new TextFile(file2, TextFile.R);
		elems = in2.readLineElems(TextFile.tab);

		while (elems != null) {

			String region = elems[0];
			String snp = elems[1] + "_" + elems[2] + "_" + elems[3];

			if (regionToSNP.containsKey(region)) {
				outf.writeln(regionToSNP.get(region) + "\t" + snp);
			}

			elems = in2.readLineElems(TextFile.tab);

		}
		in1.close();
		outf.close();


	}

}
