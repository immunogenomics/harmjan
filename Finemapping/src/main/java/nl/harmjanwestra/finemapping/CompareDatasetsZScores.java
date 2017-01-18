package nl.harmjanwestra.finemapping;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by hwestra on 12/15/16.
 */
public class CompareDatasetsZScores {

	public static void main(String[] args) {


		String bedregions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/ZScoreComparisons/regionsToCompare.bed";
		bedregions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		String outloc = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/ZScoreComparisons/";
		String[] assocfiles = new String[]{
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/RA-assoc0.3-COSMO-merged-posterior.txt.gz",
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
		};

		String[] assocfilenames = new String[]{
				"RA",
				"T1D",
		};

		String regionToGenesFile = "";

		CompareDatasetsZScores z = new CompareDatasetsZScores();
		try {
			z.run(assocfiles, assocfilenames, bedregions, regionToGenesFile, outloc);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}


	public void run(String[] associationFiles, String[] assocfilenames, String bedregions, String regionToGeneFile, String outloc) throws IOException {


		HashMap<String, String> regionToGene = new HashMap<String, String>();

		TextFile tf = new TextFile(regionToGeneFile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String region = elems[0];
			String genes = elems[1];
			regionToGene.put(region, genes);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bedregions);

		AssociationFile f = new AssociationFile();
		ArrayList<ArrayList<AssociationResult>> results = new ArrayList<>();
		for (int r = 0; r < associationFiles.length; r++) {
			results.add(f.read(associationFiles[r]));
		}

		System.out.println(regions.size() + " regionss");
		System.out.println(results.size() + " datasets");


		TextFile outall = new TextFile(outloc + "allcorrel.txt", TextFile.W);
		outall.writeln("Ds1\tDs2\tregion\tn\tcorrel\tspearman");
		TextFile out = new TextFile(outloc + "allAsssocMeta.txt", TextFile.W);
		out.writeln("Region\tgenes\tpos\trsid1\tZ1\tBeta1\tse1\tn1\tor1\tAlleles1\tminorAllele1"
				+ "\trsid2\tZ2\tBeta2\tse2\tn2\tor2\tAlleles2\tminorAllele2\tMetaZ\tMetaP\tMetaZAbs\tMetaPAbs");
		for (Feature region : regions) {
			System.out.println(region.toString());
			for (int a = 0; a < results.size(); a++) {
				ArrayList<AssociationResult> ares = results.get(a);
				ares = filter(ares, region);
				System.out.println(ares.size() + " in A dataset for region");
				HashMap<String, AssociationResult> map1 = hash(ares);
				for (int b = a + 1; b < results.size(); b++) {

					// find shared assocs
					System.out.println(a + " vs " + b);
					ArrayList<AssociationResult> bres = results.get(b);
					bres = filter(bres, region);
					System.out.println(bres.size() + " in B dataset for region");
					ArrayList<AssociationResult> shared = new ArrayList<>();
					HashMap<String, Integer> sharedIndex = new HashMap<>();

					for (AssociationResult r : bres) {
						String str = r.getSnp().toString();
						if (map1.containsKey(str)) {
							if (!sharedIndex.containsKey(str)) {
								sharedIndex.put(str, shared.size());
								shared.add(r);
							}
						}
					}
					AssociationResult[] resultA = new AssociationResult[shared.size()];

					double[] z1 = new double[shared.size()];

					for (int q = 0; q < ares.size(); q++) {
						Integer id = sharedIndex.get(ares.get(q).getSnp().toString());
						if (id != null) {
							z1[id] = ares.get(q).getZ();
							resultA[id] = ares.get(q);
						}

					}
					double[] z2 = new double[shared.size()];

					for (int q = 0; q < shared.size(); q++) {
						z2[q] = shared.get(q).getZ();
					}

//					String fname = outloc + assocfilenames[a] + "-" + assocfilenames[b] + "-" + region.toString() + ".txt";

					for (int i = 0; i < z1.length; i++) {
						String allelesA = Strings.concat(resultA[i].getSnp().getAlleles(), Strings.comma);
						String minorAlleleA = resultA[i].getSnp().getMinorAllele();
						double betaA = resultA[i].getBeta()[0];
						double orA = resultA[i].getORs()[0];
						String rsA = resultA[i].getSnp().getName();

						double seA = resultA[i].getSe()[0];
						int nA = resultA[i].getN();

						String allelesB = Strings.concat(shared.get(i).getSnp().getAlleles(), Strings.comma);
						String minorAlleleB = shared.get(i).getSnp().getMinorAllele();
						double betaB = shared.get(i).getBeta()[0];
						double orB = shared.get(i).getORs()[0];
						String rsB = shared.get(i).getSnp().getName();

						double seB = shared.get(i).getSe()[0];
						int nB = shared.get(i).getN();

						double metaZ = ZScores.getWeightedZ(new double[]{z1[i], z2[i]}, new int[]{nA, nB});
						double metaP = ZScores.zToP(metaZ);

						double metaZAbs = ZScores.getWeightedZ(new double[]{Math.abs(Math.abs(z1[i])), Math.abs(z2[i])}, new int[]{nA, nB});
						double metaPAbs = ZScores.zToP(metaZAbs);

						out.writeln(region.toString()
								+ "\t" + regionToGene.get(region.toString())
								+ "\t" + shared.get(i).getSnp().getStart()
								+ "\t" + rsA
								+ "\t" + z1[i]
								+ "\t" + betaA
								+ "\t" + seA
								+ "\t" + nA
								+ "\t" + orA
								+ "\t" + allelesA
								+ "\t" + minorAlleleA
								+ "\t" + rsB
								+ "\t" + z2[i]
								+ "\t" + betaB
								+ "\t" + seB
								+ "\t" + nB
								+ "\t" + orB
								+ "\t" + allelesB
								+ "\t" + minorAlleleB
								+ "\t" + metaZ
								+ "\t" + metaP
								+ "\t" + metaZAbs
								+ "\t" + metaPAbs);
					}

					double corr = Correlation.correlate(z1, z2);
					SpearmansCorrelation c = new SpearmansCorrelation();
					double spear = c.correlation(z1, z2);
					System.out.println(assocfilenames[a] + "-" + assocfilenames[b] + "\t" + region.toString() + "\t" + corr + "\t" + spear);

					outall.writeln(assocfilenames[a] + "-" + assocfilenames[b] + "\t" + region.toString() + "\t" + z1.length + "\t" + corr + "\t" + spear);


				}
			}
		}
		outall.close();
		out.close();


	}

	private ArrayList<AssociationResult> filter(ArrayList<AssociationResult> res, Feature region) {
		ArrayList<AssociationResult> output = new ArrayList<>();
		for (AssociationResult r : res) {
			if (r.getSnp().overlaps(region)) {
				output.add(r);
			}
		}
		return output;
	}

	private HashMap<String, AssociationResult> hash(ArrayList<AssociationResult> associationResults) {
		HashMap<String, AssociationResult> out = new HashMap<>();
		for (AssociationResult r : associationResults) {
			out.put(r.getSnp().toString(), r);
		}
		return out;
	}
}
