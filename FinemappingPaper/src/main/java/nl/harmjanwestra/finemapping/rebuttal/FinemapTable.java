package nl.harmjanwestra.finemapping.rebuttal;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class FinemapTable {
	
	public void filter(String bedin, String bedout, String assoc) throws IOException {
		AssociationFile f = new AssociationFile();
		ArrayList<AssociationResult> assocs = f.read(assoc);
		BedFileReader r = new BedFileReader();
		
		ArrayList<Feature> regions = r.readAsList(bedin);
		TextFile tf = new TextFile(bedout, TextFile.W);
		int nrsig = 0;
		for (Feature region : regions) {
			ArrayList<AssociationResult> ass = getRegionAssocs(assocs, region);
			boolean sig = false;
			for (AssociationResult a : ass) {
				if (a.getPval() < 7.5E-7) {
					sig = true;
				}
			}
			if (sig) {
				nrsig++;
				System.out.println(region.toString());
				tf.writeln(region.toBedString());
			}
		}
		tf.close();
		System.out.println(nrsig);
	}
	
	public static void main(String[] args) {
		
		String allregions = "c:/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		String out = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\META-assoc0.3-COSMO-merged-significantregions-75e7.bed";
		String assoc = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\META-assoc0.3-COSMO-merged-posterior.txt.gz";

//		FinemapTable p = new FinemapTable();
//		try {
//			p.filter(allregions, out, assoc);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

//		System.exit(-1);
//		String[] regions = new String[]{
//				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-significantregions-75e7.bed",
//				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-significantregions-75e7.bed",
//				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\META-assoc0.3-COSMO-merged-significantregions-75e7.bed",
//		};
		String[] regions = new String[]{
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\finemap\\raloci.txt",
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\finemap\\t1dloci.txt",
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\finemap\\metaloci.txt",
		};
//		String allregionsfile = "";
		String[] conditionalFile = new String[]{
				"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\finemap\\2017-11-14-conditionalVariantsRA.txt",
				"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\finemap\\2017-11-14-conditionalVariantsT1D.txt",
				"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\finemap\\2017-11-14-conditionalVariants.txt"
		};
		String[] datasets = new String[]{
				"RA",
				"T1D",
				"META"
		};
		String[] assocfiles = new String[]{
				"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
				"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
				"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\META-assoc0.3-COSMO-merged-posterior.txt.gz"
		};
		String finemapdir = "d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\finemap\\";
		int maxk = 5;
		String output = "d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\finemap\\output\\";
		
		FinemapTable t = new FinemapTable();
		try {
			t.run(regions, conditionalFile, assocfiles, datasets, finemapdir, maxk, output);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String[] regions, String[] conditionalFile, String[] assocfiles, String[] datasets, String finemapdir, int maxk, String output) throws IOException {
		
		
		ArrayList<ArrayList<AssociationResult>> assocs = new ArrayList<>();
		for (int d = 0; d < assocfiles.length; d++) {
			AssociationFile f = new AssociationFile();
			assocs.add(f.read(assocfiles[d]));
		}

//		HashSet<Feature> uniqueRegions = new HashSet<Feature>();
		ArrayList<ArrayList<Feature>> featuresPerDataset = new ArrayList<>();
		BedFileReader reader = new BedFileReader();

//		ArrayList<Feature> allregions = reader.readAsList(allregionsfile);

//		for (int d = 0; d < regions.length; d++) {
//			ArrayList<Feature> diseaseregions = reader.readAsList(regions[d]);
//			ArrayList<Feature> tmp = new ArrayList<>();
//			for (Feature f : diseaseregions) {
//				if (f.getChromosome().equals(Chromosome.TWO)) {
//					tmp.add(f);
//				}
//			}
//			featuresPerDataset.add(tmp);
////			uniqueRegions.addAll(diseaseregions);
//		}

//		ArrayList<Feature> uniqueRegionsArr = new ArrayList<>();
//		uniqueRegionsArr.addAll(uniqueRegions);
//		Collections.sort(allregions, new FeatureComparator(true));
//		Collections.sort(uniqueRegionsArr, new FeatureComparator(true));
		
		// loop over regions
		ApproximateBayesPosterior abf = new ApproximateBayesPosterior();
		for (int d = 0; d < datasets.length; d++) {
			
			TextFile out = new TextFile(output + datasets[d] + ".txt", TextFile.W);
			System.out.println("Output will be here: " + output + datasets[d] + ".txt");
			ArrayList<Feature> diseaseregions = reader.readAsList(regions[d]);

//			ArrayList<Feature> tmp = new ArrayList<>();
//			for (Feature f : diseaseregions) {
//				if (f.getChromosome().equals(Chromosome.TWO)) {
//					tmp.add(f);
//				}
//			}
//			diseaseregions = tmp;
			
			String header = "Locus\tTopAssoc\tTopAssocPval\tTopAssocPosterior\tSizeOfCredibleSet\tNConditionalVars\tConditionalVars";
			for (int k = 0; k < maxk; k++) {
				int printk = k + 1;
				header += "\t" + printk + "p(k)\t"
						+ printk + "max(pk)\t"
						+ printk + "rank(ConditionalVariants)\t"
						+ printk + "pfinemap(ConditionalVariants)\t" +
						+printk + "topConfigFINEMAP\t" +
						+printk + "p(topconfigFINEMAP)";
			}
			out.writeln(header);
			
			// loop over diseases
			for (int r = 0; r < diseaseregions.size(); r++) {
				
				// print conditional combo, if any
				Feature region = diseaseregions.get(r);
				
				// loop over k (columns)
				// META-Chr1_2406887-2785671.config
				String outln = region.toString();
				
				ApproximateBayesPosterior abp = new ApproximateBayesPosterior();
				ArrayList<AssociationResult> regionAssocs = getRegionAssocs(assocs.get(d), region);
				abp.calculatePosterior(regionAssocs);
				if (regionAssocs == null || regionAssocs.isEmpty()) {
					System.out.println("No assocs for dataset: " + datasets[d] + " in region " + region.toString());
					System.out.println(regionAssocs.size());
					System.exit(-1);
				}
				AssociationResult topassoc = getTopAssoc(regionAssocs, region, true);
				ArrayList<AssociationResult> cs = abf.createCredibleSet(regionAssocs, 0.95);
				
				if (cs == null) {
					System.out.println("cs == null for " + region);
					System.exit(-1);
				}
				
				outln += "\t" +
						topassoc.getSnp().toString() + "\t" +
						topassoc.getPval() + "\t" +
						topassoc.getPosterior() + "\t" +
						cs.size();
				
				HashMap<String, ArrayList<String>> conditionalVariants = loadConditional(conditionalFile[d]);
				ArrayList<String> conditionalVars = conditionalVariants.get(region.toString());
				
				if (conditionalVars == null) {
					// get the variant from the assoc file, I guess. leave empty for now
					outln += "\t0\t-";
				} else {
					outln += "\t" + conditionalVars.size() + "\t" + Strings.concat(conditionalVars, Strings.semicolon);
				}
				if (region.getChromosome().equals(Chromosome.TWO)) {
					if (conditionalVars == null) {
						System.out.println("No conditional vars!");
					} else {
						System.out.println("bliep! " + region + "\t" + conditionalVars.size());
					}
					System.out.println(outln);
				} else if (region.getChromosome().equals(Chromosome.THREE)) {
//					System.exit(-1);
				}
				
				for (int k = 0; k < maxk; k++) {
					String fileStr = finemapdir + (k + 1) + "/" + datasets[d] + "-" + region.toString();
					System.out.println(fileStr);
					// -- get p(k) from log
					FinemapLog log = parseLog(fileStr + ".log.gz");
					ArrayList<FinemapConfig> cfg = parseConfig(fileStr + ".config.gz");
					
					double[] pks = new double[k + 2];
					
					int maxpk = 0;
					double maxp = 0;
					for (int q = 0; q < log.pks.size(); q++) {
						pks[q] = log.pks.get(q).getRight();
						if (pks[q] > maxp) {
							maxpk = q;
							maxp = pks[q];
						}
					}
					outln += "\t" + Strings.concat(pks, Strings.semicolon) + "\t" + maxpk;
					
					// find conditional variants, if any
					if (conditionalVars != null) {
						// check configs for conditional variants
						FinemapConfig select = null;
						for (FinemapConfig c : cfg) {
							if (c.variants.size() == conditionalVars.size()) {
								int nrequal = 0;
								for (String v : c.variants) {
									for (String v2 : conditionalVars) {
										if (v2.equals(v)) {
											nrequal++;
										}
									}
								}
								if (nrequal == conditionalVars.size()) {
									// found it
									select = c;
								}
							}
						}
						if (select != null) {
							outln += "\t" + select.rank + "\t" + select.configprob;
						} else {
							outln += "\t-\t-";
						}
					} else {
						outln += "\t-\t-";
					}
					
					outln += "\t" + Strings.concat(cfg.get(0).variants, Strings.semicolon) + "\t" + cfg.get(0).configprob;
					
				}
				out.writeln(outln);
				
			}
			out.close();
		}
		
	}
	
	public ArrayList<AssociationResult> getRegionAssocs(ArrayList<AssociationResult> associationResults, Feature region) {
		ArrayList<AssociationResult> out = new ArrayList<>();
		for (AssociationResult r : associationResults) {
			if (r.getSnp().overlaps(region)) {
				out.add(r);
			}
		}
		return out;
	}
	
	private AssociationResult getTopAssoc(ArrayList<AssociationResult> associationResults, Feature region, boolean matchonposterior) {
		double max = 0;
		AssociationResult maxsnp = null;
		for (AssociationResult r : associationResults) {
			if (r.getSnp().overlaps(region)) {
				if (matchonposterior) {
					if (r.getPosterior() > max) {
						max = r.getPosterior();
						maxsnp = r;
					}
				} else {
					if (r.getLog10Pval() > max) {
						max = r.getLog10Pval();
						maxsnp = r;
					}
				}
				
			}
		}
		return maxsnp;
	}
	
	private HashMap<String, ArrayList<String>> loadConditional(String conditionalFile) throws IOException {
		
		HashMap<String, ArrayList<String>> output = new HashMap<>();
		TextFile tf = new TextFile(conditionalFile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length > 2) {
				Feature region = Feature.parseFeature(elems[0]);
				String variant = elems[2];
				ArrayList<String> vars = output.get(region.toString());
				if (vars == null) {
					vars = new ArrayList<>();
				}
				vars.add(variant);
				output.put(region.toString(), vars);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		return output;
		
	}
	
	private class FinemapLog {
		ArrayList<Pair<Integer, Double>> pks;
		public double psinglevar;
	}
	
	private class FinemapConfig {
		ArrayList<String> variants;
		int rank;
		double configprob;
		double configprobelog10bf;
	}
	
	
	private ArrayList<FinemapConfig> parseConfig(String s) throws IOException {
		System.out.println("Parsing: " + s);
		ArrayList<FinemapConfig> cfg = new ArrayList<>();
		
		TextFile tf = new TextFile(s, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {
			// rank config config_prob config_log10bf
			Integer rank = Integer.parseInt(elems[0]);
			String config = elems[1];
			String[] configElems = config.split(",");
			ArrayList<String> variants = new ArrayList<>();
			for (String e : configElems) {
				String[] meh = e.split("_");
				String variant = meh[0] + "_" + meh[1] + "_" + meh[2];
				variants.add(variant);
			}
			Double p0 = Double.parseDouble(elems[2]);
			Double p1 = Double.parseDouble(elems[3]);
			
			FinemapConfig cf = new FinemapConfig();
			cf.rank = rank;
			cf.configprob = p0;
			cf.configprobelog10bf = p1;
			cf.variants = variants;
			
			
			cfg.add(cf);
			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();
		return cfg;
		
	}
	
	private FinemapLog parseLog(String s) throws IOException {
		System.out.println("Parsing: " + s);
		TextFile tf = new TextFile(s, TextFile.R);
		String ln = tf.readLine();
		boolean probsfound = false;
		ArrayList<Pair<Integer, Double>> pk = new ArrayList<>();
		FinemapLog out = new FinemapLog();
		while (ln != null) {
			if (ln.startsWith("- log10(Evidence of at least 1 causal SNP)")) {
				
				while (ln.contains(" ")) {
					ln = ln.replaceAll(" ", "");
				}
				String[] pelems = ln.split(":");
				double psinglevar = Double.parseDouble(pelems[1]);
				out.psinglevar = psinglevar;
			} else if (ln.startsWith("- Writing output")) {
				// stop parsing
				probsfound = false;
			} else if (ln.startsWith("- Post-Pr(# of causal SNPs is k)")) {
				// parse next line ;)
				probsfound = true;
			} else if (probsfound) {
				
				// strip whitespace
				while (ln.contains(" ")) {
					ln = ln.replaceAll(" ", "");
				}
				String[] kelems = ln.split("->");
				int k = Integer.parseInt(kelems[0]);
				double p = Double.parseDouble(kelems[1]);
				Pair<Integer, Double> q = new Pair<>(k, p);
				pk.add(q);
			}
			ln = tf.readLine();
		}
		tf.close();
		
		out.pks = pk;
		return out;
	}
}
