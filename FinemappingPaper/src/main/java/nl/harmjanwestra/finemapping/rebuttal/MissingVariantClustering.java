package nl.harmjanwestra.finemapping.rebuttal;


import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.Descriptives;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.TTest;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.WilcoxonMannWhitney;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import nl.harmjanwestra.utilities.legacy.genetica.util.Primitives;
import nl.harmjanwestra.utilities.vcf.VCFVariant;


import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class MissingVariantClustering {
	
	public static void main(String[] args) {
		MissingVariantClustering c = new MissingVariantClustering();
		
		
		String[] imputedVCFs = new String[]{
				"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\INFO\\RA-COSMO.vcf.gz",
				"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\INFO\\T1D-COSMO.vcf.gz"
		};
		
		
		String disk = "d:";
		
		imputedVCFs = new String[]{
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\INFO\\RA-COSMO.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\INFO\\RA-COSMO-EAGLE-PBWT.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\INFO\\RA-EUR.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\INFO\\RA-HRC-EAGLE.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\INFO\\T1D-COSMO.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\INFO\\T1D-COSMO-EAGLE-PBWT.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\INFO\\T1D-EUR.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\INFO\\T1D-HRC-EAGLE.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\INFO\\T1D-HRC-SHAPEIT.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\INFO\\T1D-HRC-Michigan.vcf.gz",
		};
		
		String[] diseaseassoc = new String[]{
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
		};
		
		String[] diseaseout = new String[]{
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2-VariantClustering\\ra-COSMO.txt",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2-VariantClustering\\ra-COSMO-EAGLE-PBWT.txt",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2-VariantClustering\\ra-EUR.txt",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2-VariantClustering\\ra-HRC.txt",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2-VariantClustering\\t1d-COSMO.txt",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2-VariantClustering\\t1d-COSMO-EAGLE-PBWT.txt",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2-VariantClustering\\t1d-EUR.txt",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2-VariantClustering\\t1d-HRC-EAGLE.txt",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2-VariantClustering\\t1d-HRC-SHAPEIT.txt",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2-VariantClustering\\t1d-HRC-Michigan.txt",
		};
		
		String stat1kgfile = disk + "\\Data\\Ref\\1kg-maf\\stats.full.eur.txt.gz";
//		String stat1kgfile = "C:\\Data\\Ref\\1kg-maf\\stats.full.eur.txt.gz";
		String regions = disk + "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		
		
		int ldthreshold = 8;
		double impqualthreshold = 0.3;
		double mafthresholdref = 0.01;
		double mafthresholdds = 0.01;
		int nriter = 100;
		boolean samplegenomewide = false;
		boolean considerImputedButNotTestedAsMissing = false;
		
		try {
//			c.countindels("C:\\Data\\tmp\\outputsamtools.vcf.gz", regions);
//			System.exit(-1);
			c.determineIfMissingVariantsCluster(stat1kgfile,
					imputedVCFs,
					diseaseassoc,
					diseaseout,
					ldthreshold,
					mafthresholdref,
					mafthresholdds,
					impqualthreshold,
					nriter,
					regions,
					samplegenomewide,
					considerImputedButNotTestedAsMissing);
			System.exit(0);
//
			
			imputedVCFs = new String[]{
					disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\RA-cosmo-HC.txt",
					disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\RA-PBWT-HC.txt",
					disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\RA-eur-HC.txt",
					disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\RA-HRC-HC.txt",
					disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-cosmo-HC.txt",
					disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-PBWT-HC.txt",
					disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-eur-HC.txt",
					disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-HRC-HC.txt",
					disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-Michigan-HRC-HC.txt",
					disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\T1D-SHAPEIT-HRC-HC.txt",
			};
			
			String[] diseasenames = new String[]{
					"RA-COSMO",
					"RA-EAGLE-PBWT",
					"RA-EUR",
					"RA-HRC",
					"T1D-COSMO",
					"T1D-EAGLE-PBWT",
					"T1D-EUR",
					"T1D-HRC",
					"T1D-HRC-SHAPEIT",
					"T1D-HRC-Michigan",
			};
			
			String[] seqpanelvcf = new String[]{
					disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\panels\\unifiedgenotyper-maf0005-cr0950-rd10-gq30.vcf.gz-samplenamefix-mixupfix-nonmatchingremoved.vcf.gz",
					disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\panels\\hapcaller-maf0005-cr0950-rd10-gq30.vcf.gz-samplenamefix-mixupfix-nonmatchingremoved.vcf.gz"
			};
			String[] seqpanelnames = new String[]{
					"UnifiedGenotyper",
					"HaplotypeCaller"
			};
			String out = disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2-VariantMissing\\";
			
			mafthresholdds = 0;
			impqualthreshold = 0;
			boolean useAllelesForComp = false;
			boolean useIdForComp = false;
			
			c.determineMissingVariantTypes(regions,
					imputedVCFs,
					diseaseassoc,
					diseasenames,
					considerImputedButNotTestedAsMissing,
					useIdForComp,
					useAllelesForComp,
					mafthresholdref,
					mafthresholdds,
					impqualthreshold,
					seqpanelvcf,
					seqpanelnames,
					out);
//
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	private void countindels(String s, String regionsFile) throws IOException {
		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(regionsFile);
		
		TextFile tf = new TextFile(s, TextFile.R);
		String ln = tf.readLine();
		int nrindel = 0;
		int nrlns = 0;
		int nrinregion = 0;
		while (ln != null) {
			
			if (!ln.startsWith("#")) {
				VCFVariant v = new VCFVariant(ln, VCFVariant.PARSE.ALL);
				if (v.asFeature().overlaps(regions)) {
					if (v.isIndel()) {
						nrindel++;
					}
					nrinregion++;
				}
				
			}
			nrlns++;
			if (nrlns % 1000 == 0) {
				System.out.println(nrlns + "\t" + nrindel + "\t" + nrinregion);
			}
			ln = tf.readLine();
		}
		tf.close();
		
		System.out.println(nrindel);
	}
	
	public void determineImputationOutputInfoScores(String[] refpanels,
													String[] ds,
													String[] dsinputfromIC,
													String[] dsnames,
													double mafthresholdref,
													double mafthresholdds,
													double infothreshold,
													String outputfileloc
	) {
		
		// gtcaller ds referencevariants referencevariantsnotonIC variantsinds info>0.8 maf>1% info>0.8+maf>1%
		for (int d = 0; d < ds.length; d++) {
			// get a list of variants present in this dataset
			
			
			// get a list of variants for the input
			
			for (int r = 0; r < refpanels.length; r++) {
				// get a list of variants present in this reference panel
				
				// determine how many variants total
				// determine how many variants maf>threshold
				// determine how many variants info>threshold
				// determine how many variants maf>threshold+info>threshold
				// repeat for variants not on IC
				
				
			}
		}
		
	}
	
	
	protected int idcol = 0;
	protected int minorallele1col = 1;
	protected int aleleles1col = 2;
	protected int maf1col = 3;
	protected int cr1col = 4;
	protected int minorallele2col = 5;
	protected int aleleles2col = 6;
	protected int maf2col = 7;
	protected int cr2col = 8;
	protected int dfcol = 9;
	protected int samplecol = 10;
	protected int rcol = 11;
	protected int rsqlcol = 12;
	protected int betacol = 13;
	protected int secol = 14;
	protected int impqual1 = 15;
	protected int impqual2 = 16;
	
	public void determineMissingVariantTypes(String regionsFile,
											 String[] imputedVCFFiles,
											 String[] diseaseAssocFiles,
											 String[] diseasenames,
											 boolean considerImputedButNotTestedAsMissing,
											 boolean includeIdForComparison,
											 boolean includeAllelesForComparison,
											 double mafthresholdref,
											 double mafthresholdds,
											 double impqualthreshold,
											 String[] referencepanels,
											 String[] refpanelnames,
											 String outfilename
	) throws IOException {
		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(regionsFile);
		
		// determine which variants are included
		ArrayList<ArrayList<KgVariant>> imputedVariants = new ArrayList<>();
		for (int d = 0; d < imputedVCFFiles.length; d++) {
			ArrayList<KgVariant> vars = new ArrayList<>();
			TextFile tf = new TextFile(imputedVCFFiles[d], TextFile.R);
			System.out.println("parsing: " + imputedVCFFiles[d]);
			String ln = tf.readLine();
			
			while (ln != null) {
				if (!ln.startsWith("#")) {
					VCFVariant v = new VCFVariant(ln, VCFVariant.PARSE.HEADER);
					
					String af = v.getInfo().get("AF");
					if (af != null) {
						String[] split = af.split(",");
						if (split.length == 1) {
							Double maf = Double.parseDouble(af);
							if (maf > 0.5) {
								maf = 1 - maf;
							}
							double impqual = v.getImputationQualityScore();
							if (maf > mafthresholdds && impqual > impqualthreshold && v.asFeature().overlaps(regions)) {
								KgVariant vs = new KgVariant();
								vs.f = v.asSNPFeature();
								vs.f.useNameForComparison(false);
								vs.maf = maf;
								String hweps = v.getInfo().get("HWEP");
								if (hweps != null) {
									vs.hwep = Double.parseDouble(hweps);
								}
								if (vs.f.isMultiAllelic()) {
									vs.f.useAllelesForComparison(false);
								}
								
								if (maf > mafthresholdds) {
									if (!includeAllelesForComparison) {
										vs.f.useAllelesForComparison(false);
									}
									if (!includeIdForComparison) {
										vs.f.setName(null);
									}
									vars.add(vs);
								}
								
								
							}
						} else {
							// System.out.println("skipping multi-allelic variant " + v.getId() + " \t " + Strings.concat(v.getAlleles(), Strings.comma) + "\t" + af);
						}
						
					}
				}
				ln = tf.readLine();
			}
			tf.close();
			System.out.println(vars.size() + " variants included after imputation for disease " + d);
			imputedVariants.add(vars);
		}
		
		if (considerImputedButNotTestedAsMissing) {
			ArrayList<ArrayList<KgVariant>> tmpimputedVariants = new ArrayList<>();
			for (int d = 0; d < diseaseAssocFiles.length; d++) {
				AssociationFile f = new AssociationFile();
				ArrayList<KgVariant> variants = imputedVariants.get(d);
				ArrayList<AssociationResult> results = f.read(diseaseAssocFiles[d]);
				HashSet<String> assocvars = new HashSet<>();
				for (AssociationResult r : results) {
					String var = r.getSnp().getChromosome().toString() + "_" + r.getSnp().getStart();
					assocvars.add(var);
				}
				
				ArrayList<KgVariant> vars = new ArrayList<>();
				for (KgVariant v : variants) {
					String var = v.f.getChromosome().toString() + "_" + v.f.getStart();
					v.f.useNameForComparison(false);
					if (assocvars.contains(var)) {
						vars.add(v);
					}
				}
				System.out.println(variants.size() + " variants before comparing to assoc file; " + vars.size() + " after assoc filter for disease " + d);
				tmpimputedVariants.add(vars);
			}
			imputedVariants = tmpimputedVariants;
		}
		
		// loop through the sequenced variants
		// determine all possible variants in the listed regions
		ArrayList<ArrayList<KgVariant>> allRefVariants = new ArrayList<>();
		HashSet<KgVariant> allRefVariantsSet = new HashSet<>();
		for (int d = 0; d < referencepanels.length; d++) {
			String reffile = referencepanels[d];
			
			TextFile tf = new TextFile(reffile, TextFile.R);
			String ln = tf.readLine();
			
			ArrayList<KgVariant> refVariants = new ArrayList<>();
			while (ln != null) {
				if (!ln.startsWith("#")) {
					VCFVariant v = new VCFVariant(ln, VCFVariant.PARSE.ALL);
					if (v.asFeature().getChromosome().isAutosome() && v.getMAF() > mafthresholdref && v.asFeature().overlaps(regions)) {
						KgVariant var = new KgVariant();
						var.f = v.asSNPFeature();
						
						var.maf = v.getMAF();
						var.hwep = v.getHwep();
						
						var.f.useNameForComparison(false);
						if (var.f.isMultiAllelic()) {
							var.f.useAllelesForComparison(false);
						}
						if (!includeAllelesForComparison) {
							var.f.useAllelesForComparison(false);
						}
						if (!includeIdForComparison) {
							var.f.setName(null);
						}
						refVariants.add(var);
						allRefVariantsSet.add(var);
					}
				}
				ln = tf.readLine();
			}
			tf.close();
			allRefVariants.add(refVariants);
			System.out.println(refVariants.size() + " variants in " + refpanelnames[d]);
		}
		
		// count number of variants of each type per reference panel
		int[][] nrMultiAllelic = new int[referencepanels.length + 1][2];
		int[][] nrIndels = new int[referencepanels.length + 1][2];
		int[][] nrSNPs = new int[referencepanels.length + 1][2];
		
		for (int r1 = 0; r1 < referencepanels.length; r1++) {
			ArrayList<KgVariant> list1 = allRefVariants.get(r1);
			HashSet<KgVariant> presentInOtherRefs = new HashSet<>();
			for (int r2 = 0; r2 < referencepanels.length; r2++) {
				if (r2 != r1) {
					presentInOtherRefs.addAll(allRefVariants.get(r2));
				}
			}
			
			for (KgVariant v : list1) {
				if (v.f.isMultiAllelic()) {
					nrMultiAllelic[r1][0]++;
					if (!presentInOtherRefs.contains(v)) {
						nrMultiAllelic[r1][1]++;
						
					}
				} else if (v.f.isIndel()) {
					nrIndels[r1][0]++;
					if (!presentInOtherRefs.contains(v)) {
						nrIndels[r1][1]++;
					}
				} else {
					nrSNPs[r1][0]++;
					if (!presentInOtherRefs.contains(v)) {
						nrSNPs[r1][1]++;
					}
				}
			}
		}
		
		for (KgVariant v : allRefVariantsSet) {
			if (v.f.isMultiAllelic()) {
				nrMultiAllelic[referencepanels.length][0]++;
			} else if (v.f.isIndel()) {
				nrIndels[referencepanels.length][0]++;
			} else {
				nrSNPs[referencepanels.length][0]++;
			}
		}
		
		
		// compare reference panels
		int[][] refComparisonsSNPs = new int[referencepanels.length + 1][referencepanels.length + 1];
		int[][] refComparisonsIndels = new int[referencepanels.length + 1][referencepanels.length + 1];
		int[][] refComparisonsMultiAllelic = new int[referencepanels.length + 1][referencepanels.length + 1];
		for (int r1 = 0; r1 < referencepanels.length; r1++) {
			HashSet<KgVariant> set1 = new HashSet<>();
			HashSet<KgVariant> set1Indels = new HashSet<>();
			HashSet<KgVariant> set1MultiAllelic = new HashSet<>();
			
			ArrayList<KgVariant> list1 = allRefVariants.get(r1);
			for (KgVariant v1 : list1) {
				if (v1.f.isMultiAllelic()) {
					set1MultiAllelic.add(v1);
				} else if (v1.f.isIndel()) {
					set1Indels.add(v1);
				} else {
					set1.add(v1);
				}
			}
			set1.addAll(allRefVariants.get(r1));
			
			
			for (int r2 = 0; r2 < referencepanels.length; r2++) {
				ArrayList<KgVariant> list = allRefVariants.get(r2);
				for (KgVariant v : list) {
					if (set1.contains(v)) {
						if (v.f.isMultiAllelic()) {
							refComparisonsMultiAllelic[r1][r2]++;
						} else if (v.f.isIndel()) {
							refComparisonsIndels[r1][r2]++;
						} else {
							refComparisonsSNPs[r1][r2]++;
						}
					}
					
				}
			}
		}
		
		
		// compare to full list
		for (int r1 = 0; r1 < referencepanels.length + 1; r1++) {
			ArrayList<KgVariant> list = null;
			if (r1 < refpanelnames.length) {
				list = allRefVariants.get(r1);
			} else {
				list = new ArrayList<>();
				list.addAll(allRefVariantsSet);
			}
			int r2 = referencepanels.length;
			for (KgVariant v : list) {
				if (allRefVariantsSet.contains(v)) {
					if (v.f.isMultiAllelic()) {
						refComparisonsMultiAllelic[r1][r2]++;
						refComparisonsMultiAllelic[r2][r1]++;
					} else if (v.f.isIndel()) {
						refComparisonsIndels[r1][r2]++;
						refComparisonsIndels[r2][r1]++;
					} else {
						refComparisonsSNPs[r1][r2]++;
						refComparisonsSNPs[r2][r1]++;
					}
				}
			}
		}
		
		
		// determine per dataset which variants ones are missing, and what type they are
		int[][] snpsmissing = new int[imputedVariants.size()][referencepanels.length + 1];
		int[][] indelsmissing = new int[imputedVariants.size()][referencepanels.length + 1];
		int[][] multiallelicmissing = new int[imputedVariants.size()][referencepanels.length + 1];
		int[][] snpsoverlapping = new int[imputedVariants.size()][referencepanels.length + 1];
		int[][] indelsoverlapping = new int[imputedVariants.size()][referencepanels.length + 1];
		int[][] multiallelicoverlapping = new int[imputedVariants.size()][referencepanels.length + 1];
		
		
		for (int r = 0; r < referencepanels.length + 1; r++) {
			
			ArrayList<KgVariant> reflist = null;
			if (r == referencepanels.length) {
				reflist = new ArrayList<>();
				reflist.addAll(allRefVariantsSet);
			} else {
				reflist = allRefVariants.get(r);
			}
			
			for (int d = 0; d < imputedVariants.size(); d++) {
				ArrayList<KgVariant> list = imputedVariants.get(d);
				HashSet<KgVariant> presentVariants = new HashSet<>();
				for (KgVariant v : list) {
					presentVariants.add(v);
				}
				
				// now test agains the reference set
				for (KgVariant rv : reflist) {
					if (presentVariants.contains(rv)) {
						if (rv.f.isMultiAllelic()) {
							multiallelicoverlapping[d][r]++;
						} else if (rv.f.isIndel()) {
							indelsoverlapping[d][r]++;
						} else {
							snpsoverlapping[d][r]++;
						}
					} else {
						if (rv.f.isMultiAllelic()) {
							multiallelicmissing[d][r]++;
						} else if (rv.f.isIndel()) {
							indelsmissing[d][r]++;
						} else {
							snpsmissing[d][r]++;
						}
					}
				}
			}
		}
		
		// write to disk
		TextFile tf = new TextFile(outfilename + "refcomps.txt", TextFile.W);
		
		String preheader = "-\tSNPs\tUniqueSNPs\tIndels\tUniqueIndels\tMultiAllelic\tUniqueMultiAllelic";
		tf.writeln(preheader);
		for (int r = 0; r < referencepanels.length; r++) {
			tf.writeln(refpanelnames[r] + "\t" + Strings.concat(nrSNPs[r], Strings.tab)
					+ "\t" + Strings.concat(nrIndels[r], Strings.tab)
					+ "\t" + Strings.concat(nrMultiAllelic[r], Strings.tab));
		}
		tf.writeln("All\t" + Strings.concat(nrSNPs[referencepanels.length], Strings.tab)
				+ "\t" + Strings.concat(nrIndels[referencepanels.length], Strings.tab)
				+ "\t" + Strings.concat(nrMultiAllelic[referencepanels.length], Strings.tab));
		tf.writeln();
		
		TextFile dscomps = new TextFile(outfilename + "dscomps.txt", TextFile.W);
		for (int i = 0; i < 3; i++) {
			String header = "SNPs";
			int[][] outputarrRef = refComparisonsSNPs;
			
			int[][] outcomparrdsmissing = snpsmissing;
			int[][] outcomparrdsoverlap = snpsoverlapping;
			
			if (i == 1) {
				header = "Indels";
				outputarrRef = refComparisonsIndels;
				outcomparrdsmissing = indelsmissing;
				outcomparrdsoverlap = indelsoverlapping;
			}
			if (i == 2) {
				header = "MultiAllelic";
				outputarrRef = refComparisonsMultiAllelic;
				outcomparrdsmissing = multiallelicmissing;
				outcomparrdsoverlap = multiallelicoverlapping;
			}
			
			for (int d = 0; d < refpanelnames.length; d++) {
				header += "\t" + refpanelnames[d];
			}
			header += "\tAll";
			
			tf.writeln(header);
			dscomps.writeln(header);
			
			for (int d = 0; d < referencepanels.length + 1; d++) {
				String ln = null;
				int totalVarsInCat = 0;
				if (d == referencepanels.length) {
					ln = "All\t" + Strings.concat(outputarrRef[d], Strings.tab);
				} else {
					ln = refpanelnames[d] + "\t" + Strings.concat(outputarrRef[d], Strings.tab);
				}
				tf.writeln(ln);
			}
			tf.writeln();
			
			for (int d = 0; d < diseasenames.length; d++) {
				String ln = diseasenames[d] + "\t" + Strings.concat(outcomparrdsoverlap[d], Strings.tab) + "\t" + Strings.concat(outcomparrdsmissing[d], Strings.tab);
				dscomps.writeln(ln);
			}
			dscomps.writeln();
		}
		
		tf.close();
		dscomps.close();
		
		
	}
	
	
	public void determineIfMissingVariantsCluster(String ldFile1kg,
												  String[] imputedVCFFiles,
												  String[] diseaseAssocFiles,
												  String[] diseaseout,
												  int ldthreshold,
												  double mafthresholdref,
												  double mafthresholdds,
												  double infoscorethreshold,
												  int nriter,
												  String regionsFile,
												  boolean sampleGenomeWide,
												  boolean considerImputedButNotTestedAsMissing) throws IOException {
		
		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(regionsFile);
		
		// determine which variants are included
		ArrayList<ArrayList<VCFVariant>> imputedVariants = new ArrayList<>();
		for (int d = 0; d < imputedVCFFiles.length; d++) {
			ArrayList<VCFVariant> vars = new ArrayList<>();
			TextFile tf = new TextFile(imputedVCFFiles[d], TextFile.R);
			System.out.println("parsing: " + imputedVCFFiles[d]);
			String ln = tf.readLine();
			
			while (ln != null) {
				if (!ln.startsWith("#")) {
					VCFVariant v = new VCFVariant(ln, VCFVariant.PARSE.HEADER);
					String af = v.getInfo().get("AF");
					if (af != null && v.asFeature().overlaps(regions)) {
						String[] split = af.split(",");
						if (split.length == 1) {
							Double maf = Double.parseDouble(af);
							if (maf > 0.5) {
								maf = 1 - maf;
							}
							double infoscore = v.getImputationQualityScore();
							if (maf > mafthresholdds && infoscore > infoscorethreshold && v.asFeature().overlaps(regions)) {
								KgVariant vs = new KgVariant();
								vs.f = v.asSNPFeature();
								vs.maf = maf;
								String hweps = v.getInfo().get("HWEP");
								if (hweps != null) {
									vs.hwep = Double.parseDouble(hweps);
								}
								vars.add(v);
								
							}
						} else {
							// System.out.println("skipping multi-allelic variant " + v.getId() + " \t " + Strings.concat(v.getAlleles(), Strings.comma) + "\t" + af);
						}
						
					}
				}
				ln = tf.readLine();
			}
			tf.close();
			System.out.println(vars.size() + " variants included after imputation for disease " + d);
			imputedVariants.add(vars);
		}
		
		if (considerImputedButNotTestedAsMissing) {
			ArrayList<ArrayList<VCFVariant>> tmpimputedVariants = new ArrayList<>();
			for (int d = 0; d < diseaseAssocFiles.length; d++) {
				AssociationFile f = new AssociationFile();
				ArrayList<VCFVariant> variants = imputedVariants.get(d);
				ArrayList<AssociationResult> results = f.read(diseaseAssocFiles[d]);
				HashSet<String> assocvars = new HashSet<>();
				for (AssociationResult r : results) {
					String var = r.getSnp().getChromosome().toString() + "_" + r.getSnp().getStart();
					assocvars.add(var);
				}
				
				ArrayList<VCFVariant> vars = new ArrayList<>();
				for (VCFVariant v : variants) {
					String var = v.asFeature().getChromosome().toString() + "_" + v.getPos();
					if (assocvars.contains(var)) {
						vars.add(v);
					}
				}
				System.out.println(variants.size() + " variants before comparing to assoc file; " + vars.size() + " after for disease " + d);
				tmpimputedVariants.add(vars);
			}
			imputedVariants = tmpimputedVariants;
		}
		
		// determine which variants are in the reference panel
		ArrayList<KgVariant> kgVariantsInRegions = new ArrayList<>();
		ArrayList<KgVariant> kgVariantsNotInRegions = new ArrayList<>();
		
		TextFile tf = new TextFile(ldFile1kg, TextFile.R);
		tf.readLine();
		String ln = tf.readLine();
		while (ln != null) {
			String[] elems = Strings.tab.split(ln);
			// Variant MAF     HWEP    0       1       2       3       4       5       6       7       8       9
			String variant = new String(elems[0]);
			SNPFeature feature = SNPFeature.parseSNPFeature(variant);
			Double maf = Double.parseDouble(elems[1]);
			if (maf > mafthresholdref) {
				Double hwep = Double.parseDouble(elems[2]);
				int nproxies = 0;
				for (int b = 3 + ldthreshold; b < elems.length; b++) {
					nproxies += Integer.parseInt(elems[b]);
				}
				KgVariant v = new KgVariant();
				v.f = feature;
				v.nproxies = nproxies;
				v.maf = maf;
				v.hwep = hwep;
				
				if (feature.overlaps(regions)) {
					kgVariantsInRegions.add(v);
				} else {
					kgVariantsNotInRegions.add(v);
				}
			}
			
			ln = tf.readLine();
		}
		tf.close();
		
		System.out.println(kgVariantsInRegions.size() + " variants in regions");
		System.out.println(kgVariantsNotInRegions.size() + " variants not in regions");
		
		// compare imputed variants with reference variantset
		for (int d = 0; d < imputedVariants.size(); d++) {
			System.out.println("Assessing dataset " + d);
			// make a list of variants that we actually included in the datasets,
			// so we can link them to the reference variants
			ArrayList<VCFVariant> diseasevariants = imputedVariants.get(d);
			HashSet<String> variantIds = new HashSet<String>();
			for (VCFVariant v : diseasevariants) {
				variantIds.add(v.asFeature().getChromosome().toString() + "_" + v.getPos());
			}
			
			// 2. determine which variants are missing after imputation
			ArrayList<KgVariant> allMissedVariants = new ArrayList<>();
			ArrayList<KgVariant> allIncludedVariants = new ArrayList<>();
			for (KgVariant v : kgVariantsInRegions) {
				String id = v.f.getChromosome().toString() + "_" + v.f.getStart();
				if (!variantIds.contains(id)) {
					allMissedVariants.add(v);
				} else {
					allIncludedVariants.add(v);
				}
			}
			System.out.println(variantIds.size() + " in imputed set. " + kgVariantsInRegions.size() + " variants in ref. " + allMissedVariants.size() + " missing.");
			
			// sort variants because sorting is awesome
			Collections.sort(allMissedVariants);
			Collections.sort(allIncludedVariants);
			
			// 3. measure distance between variants in regions that are missing
			
			HashSet<KgVariantPair> missingVariantPairs = new HashSet<>();
			
			// there are more than 1 missing variants, so write some output
			String header = "region" +
					"\tRegionDistanceMean" +
					"\tRegionDistanceVar" +
					"\tRegionDistanceN" +
					"\tMeanMatchedNullDistance" +
					"\tMeanMatchedNullVariance" +
					"\t%SignificantDistanceDiff" +
					"\t%SmallerDistance" +
					"\t%SmallerDistanceAndSignificant";
			TextFile outtf = new TextFile(diseaseout[d] + ".txt", TextFile.W);
			if (allMissedVariants.size() < 2) {
				outtf.writeln("Not enough missing variants");
			} else {
				outtf.writeln(header);
				ArrayList<Feature> regionsToTest = new ArrayList<>();
				for (Feature f : regions) {
					regionsToTest.add(f);
				}
				
				for (Feature region : regionsToTest) {
					ArrayList<Integer> distances = new ArrayList<>();
					ArrayList<KgVariant> missingKgVariantsInRegion = filter(allMissedVariants, region);
					ArrayList<KgVariant> variantsToSampleFrom = null;
					if (sampleGenomeWide) {
						variantsToSampleFrom = kgVariantsInRegions;
					} else {
						variantsToSampleFrom = filter(kgVariantsInRegions, region);
					}
					
					if (missingKgVariantsInRegion.size() < 2) {
						String out = region.toString()
								+ "\t" + 0
								+ "\t" + 0
								+ "\t" + 0
								+ "\t" + 0
								+ "\t" + 0
								+ "\t" + 0
								+ "\t" + 0
								+ "\t" + 0;
						outtf.writeln(out);
					} else {
						
						
						for (int v = 0; v < missingKgVariantsInRegion.size(); v++) {
							KgVariant var = missingKgVariantsInRegion.get(v);
							KgVariant neighbor = null;
							Integer distance = null;
							if (v == 0) {
								// nearest is the next
								neighbor = missingKgVariantsInRegion.get(v + 1);
								if (neighbor.f.getChromosome().equals(var.f.getChromosome())) {
									// measure distance
									distance = Math.abs(neighbor.f.getStart() - var.f.getStart());
								}
							} else if (v == missingKgVariantsInRegion.size() - 1) {
								// previous one is nearest
								neighbor = missingKgVariantsInRegion.get(v - 1);
								if (neighbor.f.getChromosome().equals(var.f.getChromosome())) {
									// measure distance
									distance = Math.abs(neighbor.f.getStart() - var.f.getStart());
								}
							} else {
								
								// either previous or next is nearest
								KgVariant neighbor1 = missingKgVariantsInRegion.get(v - 1);
								Integer d1 = null;
								Integer d2 = null;
								if (neighbor1.f.getChromosome().equals(var.f.getChromosome())) {
									// measure distance
									d1 = Math.abs(neighbor1.f.getStart() - var.f.getStart());
								}
								KgVariant neighbor2 = missingKgVariantsInRegion.get(v + 1);
								if (neighbor2.f.getChromosome().equals(var.f.getChromosome())) {
									// measure distance
									d2 = Math.abs(neighbor2.f.getStart() - var.f.getStart());
								}
								
								if (d1 != null && d2 != null) {
									distance = Math.min(d1, d2);
									if (distance.equals(d1)) {
										neighbor = neighbor1;
									} else {
										neighbor = neighbor2;
									}
								} else if (d1 != null) {
									distance = d1;
									neighbor = neighbor1;
								} else if (d2 != null) {
									distance = d2;
									neighbor = neighbor2;
								}
							}
							
							if (distance != null) {
								// got ourselves a missing variant with some friends
								KgVariantPair p = new KgVariantPair(var, neighbor);
								if (!missingVariantPairs.contains(p)) {
									distances.add(distance);
									missingVariantPairs.add(p); // make sure we're only counting pairs once
								}
							}
						}
						
						// determine difference somehow
						double[] a = new double[distances.size()];
						for (int q = 0; q < a.length; q++) {
							a[q] = distances.get(q);
						}
						double regionmean = Descriptives.mean(a);
						double regionvariance = Descriptives.variance(a);
						int nriterssignificant = 0;
						int nritersdistancesmaller = 0;
						int nritersdistancesmallerandsignificant = 0;
						ArrayList<Double> nullmeans = new ArrayList<>();
						
						int i = 0;
						while (i < nriter) {
							// now measure distance between matched snps
							ArrayList<Integer> nullDistance = new ArrayList<>();
							System.out.println(region.toString() + " - " + missingKgVariantsInRegion.size() + " missing variants " + variantsToSampleFrom.size() + " variants to sample from..");
							nullDistance = matchVariants(variantsToSampleFrom, missingKgVariantsInRegion);
							
							double[] b = new double[nullDistance.size()];
							for (int q = 0; q < b.length; q++) {
								b[q] = nullDistance.get(q);
							}
							
							double nullmean = Descriptives.mean(b);
							nullmeans.add(nullmean);
							double nullvariance = Descriptives.variance(b);
							WilcoxonMannWhitney mwm = new WilcoxonMannWhitney();
							double pwilcoxon = mwm.returnWilcoxonMannWhitneyPValue(a, b);
							double pstudent = TTest.test(a, b);
							
							if (nullmean < regionmean) {
								nritersdistancesmaller++;
							}
							if (pwilcoxon < 0.05 / regions.size()) {
								nriterssignificant++;
								if (nullmean < regionmean) {
									nritersdistancesmallerandsignificant++;
								}
							}
							i++;
						}
						
						
						double d2 = (double) nritersdistancesmaller / nriter;
						double d1 = (double) nriterssignificant / nriter;
						double d3 = (double) nritersdistancesmallerandsignificant / nriter;
						
						String out = region.toString()
								+ "\t" + regionmean
								+ "\t" + regionvariance
								+ "\t" + distances.size()
								+ "\t" + Descriptives.mean(Primitives.toPrimitiveArr(nullmeans.toArray(new Double[0])))
								+ "\t" + Descriptives.variance(Primitives.toPrimitiveArr(nullmeans.toArray(new Double[0])))
								+ "\t" + d1
								+ "\t" + d2
								+ "\t" + d3;
						
						outtf.writeln(out);
					}
				}
			}
			outtf.close();
		}
	}
	
	private ArrayList<KgVariant> filter(ArrayList<KgVariant> kgVariantsInRegions, Feature region) {
		ArrayList<KgVariant> output = new ArrayList<>();
		for (KgVariant v : kgVariantsInRegions) {
			if (v.f.overlaps(region)) {
				output.add(v);
			}
		}
		return output;
	}
	
	public ArrayList<Integer> matchVariants(ArrayList<KgVariant> refVariantsToSampleFrom, ArrayList<KgVariant> variantsToMatch) {
		ArrayList<Integer> nullDistance = new ArrayList<>();
		
		// bin the included variants
		Bins b = new Bins(100);
		b.construct(refVariantsToSampleFrom);
		
		HashMap<KgVariant, Integer> ktoi = new HashMap<KgVariant, Integer>();
		for (int q = 0; q < refVariantsToSampleFrom.size(); q++) {
			ktoi.put(refVariantsToSampleFrom.get(q), q);
		}
		
		// iterate the variants we need to match
		for (KgVariant var1 : variantsToMatch) {
			// search variants that were included
			
			System.out.println("INPUT: " + var1.f.toString() + " - maf: " + var1.maf + "\tld: " + var1.nproxies);
			KgVariant matched = b.match(var1);
			if (matched != null) {
				System.out.println("MATCH " + matched.f.toString() + " - maf: " + matched.maf + "\tld: " + matched.nproxies);
			} else {
				System.out.println("MATCH NOT FOUND");
			}
			
			
			if (matched != null) {
				// get the variant index for variant matched
				Integer id = ktoi.get(matched);
				KgVariant neighbor = null;
				
				// get this variant's neighbor
				if (id == refVariantsToSampleFrom.size() - 1) {
					// can only go down from here
					neighbor = refVariantsToSampleFrom.get(id - 1);
				} else if (id
						== 0) {
					// reach the stars! there's only up from the bottom
					neighbor = refVariantsToSampleFrom.get(1);
				} else {
					// could go up or down.. select the closest I guess?
					KgVariant neigh1 = refVariantsToSampleFrom.get(id - 1);
					int d1 = Math.abs(matched.f.getStart() - neigh1.f.getStart());
					KgVariant neigh2 = refVariantsToSampleFrom.get(id + 1);
					int d2 = Math.abs(matched.f.getStart() - neigh2.f.getStart());
					if (neigh1.f.getChromosome().equals(matched.f.getChromosome()) && neigh2.f.getChromosome().equals(matched.f.getChromosome())) {
						// both neighbors are on same chr as our matched candidate
						if (d1 > d2) {
							neighbor = neigh2;
						} else {
							neighbor = neigh1;
						}
					} else if (neigh1.f.getChromosome().equals(matched.f.getChromosome())) {
						neighbor = neigh1;
					} else if (neigh2.f.getChromosome().equals(matched.f.getChromosome())) {
						neighbor = neigh2;
					}
				}
				
				// check whether the chromosome is identical
				if (neighbor != null && neighbor.f.getChromosome().equals(matched.f.getChromosome())) {
					// measure distance
					nullDistance.add(Math.abs(neighbor.f.getStart() - matched.f.getStart()));
				}
			}
		}
		
		return nullDistance;
	}
	
	class Bins {
		// [ld][maf][variantlist]
		ArrayList<ArrayList<ArrayList<KgVariant>>> bins = new ArrayList<>();
		int nrbins = 0;
		
		public Bins(int nrBins) {
			nrbins = nrBins;
			for (int i = 0; i < nrBins; i++) {
				ArrayList<ArrayList<KgVariant>> q = new ArrayList<>();
				for (int j = 0; j < nrBins; j++) {
					q.add(new ArrayList<>());
				}
				bins.add(q);
			}
		}
		
		public void construct(ArrayList<KgVariant> variants) {
			for (KgVariant v : variants) {
				int mafbin = (int) Math.floor(v.maf * 2 * nrbins);
				if (mafbin >= nrbins) {
					mafbin = nrbins - 1;
				}
				int ldbin = v.nproxies;
				if (ldbin >= nrbins) {
					ldbin = nrbins - 1;
				}
				bins.get(ldbin).get(mafbin).add(v);
			}
		}
		
		public KgVariant match(KgVariant v) {
			int mafbin = (int) Math.floor(v.maf * 2 * nrbins);
			if (mafbin >= nrbins) {
				mafbin = nrbins - 1;
			}
			int ldbin = v.nproxies;
			if (ldbin >= nrbins) {
				ldbin = nrbins - 1;
			}
			
			
			ArrayList<KgVariant> matched = bins.get(ldbin).get(mafbin);
			
			// if the matched set is empty,
			// we want to search in similar bins
			int tmpmafbin = mafbin;
			int tmpldbin = ldbin;
			
			int maxDeltaMaf = 10;
			int maxDeltaLD = 10;
			int ldIteration = 0;
			int nrtimeldup = 0;
			int nrtimelddown = 0;
			boolean reachedLDBottom = false;
			boolean reachedLDTop = false;
			while (ldIteration < maxDeltaLD) {
				
				// try another LD bin, if we didn't find anything in this ldbins' mafbins
				if (ldIteration > 0) {
					// go up or down, depending on the iteration,
					// unless we've already been at the top or bottom of the binlist
					tmpldbin = ldbin;
					if ((ldIteration % 2 == 0 || reachedLDBottom) && !reachedLDTop) {
						tmpldbin += nrtimeldup;
					} else if (!reachedLDBottom) {
						tmpldbin -= nrtimelddown;
					} else if (reachedLDBottom && reachedLDTop) {
						// this should never happen
						break;
					}
					
					if (tmpldbin <= 0) {
						reachedLDBottom = true;
						tmpldbin = 0;
					}
					if (tmpldbin >= nrbins - 1) {
						reachedLDTop = true;
						tmpldbin = nrbins - 1;
					}
				}
				
				
				// go up or down one mafbin, try maxdelta steps
				int nrtimesup = 1;
				int nrtimesdown = 1;
				boolean reachedTop = false;
				boolean reachedBottom = false;
				for (int mafIteration = 1; mafIteration < maxDeltaMaf + 1; mafIteration++) {
					tmpmafbin = mafbin; // reset the mafbin
					if ((mafIteration % 2 == 0 || reachedBottom) && !reachedTop) {
						// go up
						tmpmafbin += nrtimesup;
						nrtimesup++;
					} else if (!reachedBottom) {
						// go down
						tmpmafbin -= nrtimesdown;
						nrtimesdown++;
					} else if (reachedTop && reachedBottom) {
						matched = new ArrayList<>(); // return empty handed
						break;
					}
					
					if (tmpmafbin <= 0) {
						reachedBottom = true;
						tmpmafbin = 0;
					}
					if (tmpmafbin >= nrbins - 1) {
						reachedTop = true;
						tmpmafbin = nrbins - 1;
					}
					
					matched = bins.get(tmpldbin).get(tmpmafbin);
					if (!matched.isEmpty()) {
						break;
					}
				}
				
				ldIteration++;
			}
			
			// now pick a random snp from the list
			if (!matched.isEmpty()) {
				int index = (int) Math.floor(Math.random() * matched.size());
				if (index == matched.size()) {
					index = matched.size() - 1;
				}
				return matched.get(index);
			} else {
				return null;
			}
			
		}
		
	}
	
	
}
