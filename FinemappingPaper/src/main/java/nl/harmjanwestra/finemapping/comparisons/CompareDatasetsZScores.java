package nl.harmjanwestra.finemapping.comparisons;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.util.Primitives;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.Correlation;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.ZScores;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by hwestra on 12/15/16.
 */
public class CompareDatasetsZScores {
	
	public static void main(String[] args) {


//		String bedregions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/ZScoreComparisons/regionsToCompare.bed";
//		bedregions = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
//		String outloc = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/wosharedcontrols/";
//		String[] assocfiles = new String[]{
//				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/wosharedcontrols/RA-assoc0.3-COSMO-merged-posterior.txt.gz",
//				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/wosharedcontrols/T1D-assoc0.3-COSMO-merged-posterior.txt.gz"
//		};
//
//		String[] assocfilenames = new String[]{
//				"RA",
//				"T1D",
//		};
//
//		String regionToGenesFile = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/AllLoci-GenesPerLocus.txt";
		
		
		String disk = "d:";
		String bedregions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/ZScoreComparisons/regionsToCompare.bed";
		bedregions = disk + "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		String outloc = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\woSharedSamples\\";
		String[] assocfiles = new String[]{
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\woSharedSamples\\RA-assoc0.3-COSMO-merged.txt.gz",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\woSharedSamples\\T1D-assoc0.3-COSMO-merged.txt.gz"
		};
//
		String[] assocfilenames = new String[]{
				"RA",
				"T1D",
		};
		
		String regionToGenesFile = disk + "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/AllLoci-GenesPerLocus.txt";
		
		bedregions = "D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\META-assoc0.3-COSMO-merged-significantregions-75e7.bed";
		assocfiles = new String[]{
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\META-assoc0.3-COSMO-merged-posterior.txt.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\multinomial\\METAPCA-assoc0.3-COSMO-merged.txt.gz"
		};
		assocfilenames = new String[]{
				"META-covars",
				"META-PCA",
		};
		
		outloc = disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\multinomial\\pcanonpcacomp\\pcanonpcacomp";
		
		CompareDatasetsZScores z = new CompareDatasetsZScores();
		
		bedregions = disk + "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		assocfiles = new String[]{
				"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged.txt.gz",
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\hrcimputedassoc\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
		};
		outloc = disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\hrcimputedassoc\\merge\\RAZscoreComp";
		try {
			z.run(assocfiles, assocfilenames, bedregions, regionToGenesFile, outloc);
		} catch (IOException e) {
			e.printStackTrace();
		}
		assocfiles = new String[]{
				"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged.txt.gz",
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\hrcimputedassoc\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz"
		};
		outloc = disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\hrcimputedassoc\\merge\\T1DZscoreComp";
		
		
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
			if (elems.length == 2) {
				String region = elems[0];
				String genes = elems[1];
				regionToGene.put(region, genes);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println(regionToGene.size() + " region to genes");
		
		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bedregions);
		
		AssociationFile f = new AssociationFile();
		ArrayList<ArrayList<AssociationResult>> resultsPerDataset = new ArrayList<>();
		for (int r = 0; r < associationFiles.length; r++) {
			resultsPerDataset.add(f.read(associationFiles[r]));
		}
		
		System.out.println(regions.size() + " regionss");
		System.out.println(resultsPerDataset.size() + " datasets");
		
		
		TextFile outall = new TextFile(outloc + "allcorrel.txt", TextFile.W);
		outall.writeln("Ds1-Ds2\tregion\tgenes\tn\tcorrel\tspearman\tLowestPInDs1\tLowestPInDs2\tSignificantInBoth");
		TextFile out = new TextFile(outloc + "allAsssocMeta.txt", TextFile.W);
		out.writeln("Region" +
				"\tgenes" +
				"\tpos" +
				"\trsid1" +
				
				"\tBeta1" +
				"\tse1" +
				"\tZ1" +
				"\tp1" +
				"\tn1" +
				"\tor1" +
				"\tAlleles1" +
				"\tminorAllele1" +
				"\trsid2" +
				"\tBeta2" +
				"\tse2" +
				"\tZ2" +
				"\tp2" +
				"\tn2" +
				"\tor2" +
				"\tAlleles2" +
				"\tminorAllele2" +
				"\tMetaZ" +
				"\tMetaP" +
				"\tMetaZAbs" +
				"\tMetaPAbs");
		
		double pvalthreshold = 7.5E-7;
		for (Feature region : regions) {
			
			
			boolean significantInBoth = false;
			double lowestPInA = 1;
			double lowestPInB = 1;
			
			System.out.println(region.toString());
			for (int a = 0; a < resultsPerDataset.size(); a++) {
				ArrayList<AssociationResult> ares = resultsPerDataset.get(a);
				ares = filter(ares, region);
				for (AssociationResult r : ares) {
					if (r.getPval() < lowestPInA) {
						lowestPInA = r.getPval();
					}
				}
				System.out.println(ares.size() + " in A dataset for region");
				HashMap<String, AssociationResult> map1 = hash(ares);
				for (int b = a + 1; b < resultsPerDataset.size(); b++) {
					
					// find shared assocs
					System.out.println(a + " vs " + b);
					ArrayList<AssociationResult> bres = resultsPerDataset.get(b);
					bres = filter(bres, region);
					for (AssociationResult r : bres) {
						if (r.getPval() < lowestPInB) {
							lowestPInB = r.getPval();
						}
					}
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
					
					if (!shared.isEmpty()) {
						
						Pair<double[], double[]> zs = removeNaN(z1, z2);
						z1 = zs.getLeft();
						z2 = zs.getRight();
						
						for (int i = 0; i < z1.length; i++) {
							String allelesA = Strings.concat(resultA[i].getSnp().getAlleles(), Strings.comma);
							String minorAlleleA = resultA[i].getSnp().getMinorAllele();
							double betaA = resultA[i].getBeta()[0][0];
							double orA = resultA[i].getORs()[0][0];
							String rsA = resultA[i].getSnp().getName();
							double seA = resultA[i].getSe()[0][0];
							int nA = resultA[i].getN();

//							if (betaA < 0 && z1[i] > 0 || z1[i] < 0 && betaA > 0) {
//								double z = resultA[i].getZ();
//								System.out.println(betaA + "\t" + seA + "\t" + z + "\t" + z1[i]);
//								System.exit(-1);
//							}
							
							
							String allelesB = Strings.concat(shared.get(i).getSnp().getAlleles(), Strings.comma);
							String minorAlleleB = shared.get(i).getSnp().getMinorAllele();
							double betaB = shared.get(i).getBeta()[0][0];
							double orB = shared.get(i).getORs()[0][0];
							String rsB = shared.get(i).getSnp().getName();
							
							double seB = shared.get(i).getSe()[0][0];
							int nB = shared.get(i).getN();
							
							double metaZ = ZScores.getWeightedZ(new double[]{z1[i], z2[i]}, new int[]{nA, nB});
							double metaP = ZScores.zToP(metaZ);
							
							double metaZAbs = ZScores.getWeightedZ(new double[]{Math.abs(Math.abs(z1[i])), Math.abs(z2[i])}, new int[]{nA, nB});
							double metaPAbs = ZScores.zToP(metaZAbs);
							
							
							out.writeln(region.toString()
									+ "\t" + regionToGene.get(region.toString())
									+ "\t" + shared.get(i).getSnp().getStart()
									+ "\t" + rsA
									+ "\t" + betaA
									+ "\t" + seA
									+ "\t" + z1[i]
									+ "\t" + resultA[i].getPval()
									+ "\t" + nA
									+ "\t" + orA
									+ "\t" + allelesA
									+ "\t" + minorAlleleA
									+ "\t" + rsB
									+ "\t" + betaB
									+ "\t" + seB
									+ "\t" + z2[i]
									+ "\t" + shared.get(i).getPval()
									+ "\t" + nB
									+ "\t" + orB
									+ "\t" + allelesB
									+ "\t" + minorAlleleB
									+ "\t" + metaZ
									+ "\t" + metaP
									+ "\t" + metaZAbs
									+ "\t" + metaPAbs
							);
						}
						
						
						double corr = Correlation.correlate(z1, z2);
						SpearmansCorrelation c = new SpearmansCorrelation();
						double spear = c.correlation(z1, z2);
						System.out.println(assocfilenames[a] + "-" + assocfilenames[b] + "\t" + region.toString() + "\t" + corr + "\t" + spear + "\t" + significantInBoth);
						
						outall.writeln(assocfilenames[a] + "-" + assocfilenames[b]
								+ "\t" + region.toString()
								+ "\t" + regionToGene.get(region.toString())
								+ "\t" + z1.length
								+ "\t" + corr
								+ "\t" + spear
								+ "\t" + lowestPInA
								+ "\t" + lowestPInB
								+ "\t" + significantInBoth);
						
						
					}
					
					
				}
			}
		}
		outall.close();
		out.close();
		
		
	}
	
	private Pair<double[], double[]> removeNaN(double[] z1, double[] z2) {
		ArrayList<Double> d1 = new ArrayList<>();
		ArrayList<Double> d2 = new ArrayList<>();
		
		for (int i = 0; i < z1.length; i++) {
			if (!Double.isNaN(z1[i]) && !Double.isNaN(z2[i])) {
				d1.add(z1[i]);
				d2.add(z2[i]);
			}
		}
		double[] o1 = Primitives.toPrimitiveArr(d1.toArray(new Double[0]));
		double[] o2 = Primitives.toPrimitiveArr(d2.toArray(new Double[0]));
		return new Pair<double[], double[]>(o1, o2);
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
