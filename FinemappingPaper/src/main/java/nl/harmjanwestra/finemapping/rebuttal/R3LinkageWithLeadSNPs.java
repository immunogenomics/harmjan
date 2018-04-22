package nl.harmjanwestra.finemapping.rebuttal;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Triple;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFTabix;
import nl.harmjanwestra.utilities.vcf.VCFVariant;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

public class R3LinkageWithLeadSNPs {
	
	public static void main(String[] args) {
		String disk = "d:";
		String[] diseasefiles = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
				"D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz"
		};
		
		String[] gwasfiles = new String[]{
				disk + "\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Okada\\RA.OKADA.gz",
				disk + "/Sync/SyncThing/Postdoc/2016-03-RAT1D-Finemapping/Data/ImmunoBase/hg19_gwas_ic_t1d_onengut_meta_4_19_1.tab.gz"
		};
//		String[] leadsnpfiles = new String[]{
//				"D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\leadsnps\\annot_2013_ORs.txt",
//				"D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\leadsnps\\t1d.txt"
//		};
		
		String[] leadsnpfiles = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\leadsnps\\2018-04-15-ComparisonsToMakeRA.txt",
				"D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\leadsnps\\2018-04-15-ComparisonsToMakeT1D.txt"
		};
		
		String[] regionfiles = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-significantregions-75e7.bed",
				"D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-significantregions-75e7.bed"
		};
		String tabix = disk + "\\Sync\\SyncThing\\Data\\Ref\\1kg\\ALL.chrCHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";
		String samplefile = disk + "\\Sync\\SyncThing\\Data\\Ref\\1kg-europeanpopulations.txt.gz";
//		String[] outfiles = new String[]{
//				"D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\comparisontoimmunobase2\\ra.txt",
//				"D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\comparisontoimmunobase2\\t1d.txt"
//		};
		String[] outfiles = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\comparisontoimmunobase2\\2018-04-15-ra.txt",
				"D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\comparisontoimmunobase2\\2018-04-15-t1d.txt"
		};
		
		R3LinkageWithLeadSNPs r = new R3LinkageWithLeadSNPs();
		try {
			r.runNew(diseasefiles, gwasfiles, leadsnpfiles, regionfiles, tabix, samplefile, outfiles);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void runOld(String[] diseasefiles, String[] gwasfiles, String[] leadsnpfiles, String[] regionfiles, String tabixPrefix, String samplefilterfile, String[] outfiles) throws IOException {
		
		
		QTLOverlap q = new QTLOverlap();
		AssociationFile af = new AssociationFile();
		BedFileReader bf = new BedFileReader();
		
		for (int d = 0; d < diseasefiles.length; d++) {
			ArrayList<Feature> regions = bf.readAsList(regionfiles[d]);
			ArrayList<AssociationResult> finemapresults = af.readRegions(diseasefiles[d], regions);
			System.out.println(finemapresults.size() + " finemap results");
			ArrayList<AssociationResult> gwasresults = af.readRegions(gwasfiles[d], regions);
			System.out.println(gwasresults.size() + " gwas results");
			Collections.sort(regions, new FeatureComparator(true));
			
			// load lead snps
			ArrayList<Feature> leadsnps = new ArrayList<>();
			String leadsnpfile = leadsnpfiles[d];
			TextFile tf = new TextFile(leadsnpfile, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				Chromosome chr = Chromosome.parseChr(elems[0]);
				Integer pos = Integer.parseInt(elems[1]);
				Feature f = new SNPFeature(chr, pos, pos);
				if (elems.length > 2) {
					f.setName(elems[2]);
				}
				leadsnps.add(f);
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			
			// determine whether the top variants are in
			TextFile outf = new TextFile(outfiles[d], TextFile.W);
			String header = "Region"
					+ "\tPreviouslyReportedLeadSNP"
					+ "\tFinemap-TopAssociatedVariant"
					+ "\tAFControls"
					+ "\tAfCases"
					+ "\tImpQual"
					+ "\tAlleles"
					+ "\tMinorAllele"
					+ "\tOR"
					+ "\tPval"
					+ "\tleadSNPinGWASData-snp"
					+ "\talleles"
					+ "\tminor"
					+ "\tor"
					+ "\tpval"
					+ "\tdistance"
					+ "\tdprime"
					+ "\trsq"
					+ "\tfinemaptopsnpingwasdata-snp"
					+ "\talleles"
					+ "\tminor"
					+ "\tor"
					+ "\tpval"
					+ "\tleadsnpinfmdata-snp"
					+ "\tAfControls"
					+ "\tAfCases"
					+ "\tImpQual"
					+ "\talleles"
					+ "\tminor"
					+ "\tor"
					+ "\tpval";
			outf.writeln(header);
			for (int r = 0; r < regions.size(); r++) {
				Feature region = regions.get(r);
				String tabixFile = tabixPrefix.replaceAll("CHR", "" + region.getChromosome().getNumber());
				VCFTabix t = new VCFTabix(tabixFile);
				
				boolean[] filter = t.getSampleFilter(samplefilterfile);
				ArrayList<VCFVariant> allvariants = t.getAllVariants(region, filter);
				System.out.println(allvariants.size() + " variants in region..");
				
				ArrayList<AssociationResult> regionFinemapResults = getAssoc(finemapresults, region);
				ArrayList<AssociationResult> regionGWASResults = getAssoc(gwasresults, region);
				
				// determine lead snps in original analysis
				ArrayList<Feature> regionleadsnps = getLeadSNPs(leadsnps, region);
				
				// find the top assocaiated snp in the region
				AssociationResult topassoc = null;
				AssociationResult topassocInGWAS = null;
				
				double maxp = 0;
				for (AssociationResult a : regionFinemapResults) {
					if (a.getLog10Pval() > maxp) {
						topassoc = a;
						maxp = a.getLog10Pval();
					}
				}
				
				for (AssociationResult a : regionGWASResults) {
					if (a.getSnp().overlaps(topassoc.getSnp())) {
						topassocInGWAS = a;
					}
				}
				
				for (Feature regionLeadSNP : regionleadsnps) {
					// find the matching assoc
					AssociationResult leadSNPInGWASdat = null;
					AssociationResult leadSNPInFinemapData = null;
					
					
					for (AssociationResult a : regionFinemapResults) {
						if (a.getSnp().overlaps(regionLeadSNP)) {
							// found it!
							leadSNPInFinemapData = a;
						}
					}
					
					for (AssociationResult a : regionGWASResults) {
						if (a.getSnp().overlaps(regionLeadSNP)) {
							// found it!
							leadSNPInGWASdat = a;
						}
					}
					
					/*
					Region
					PreviouslyReportedLeadSNP
					Finemap-TopAssociatedVariant
					AFControls
					AfCases
					ImpQual
					Alleles
					MinorAllele
					OR
					Pval
					 */
					
					
					String out = region.toString()
							+ "\t" + regionLeadSNP.toString()
							+ "\t" + topassoc.getSnp().toString()
							+ "\t" + topassoc.getSnp().getAFControls()
							+ "\t" + topassoc.getSnp().getAFCases()
							+ "\t" + topassoc.getSnp().getImputationQualityScore()
							+ "\t" + Strings.concat(topassoc.getSnp().getAlleles(), Strings.comma)
							+ "\t" + topassoc.getSnp().getMinorAllele()
							+ "\t" + topassoc.getORs()[0][0]
							+ "\t" + topassoc.getLog10Pval();
					
					
					/*
					leadSNPinGWASData-snp
					alleles
					minor
					or
					pval
					distance
					dprime
					rsq
					 */
					if (leadSNPInGWASdat != null) {
						// determine LD with topassoc
						double ld = 0;
						VCFVariant v1 = q.getVariant(topassoc.getSnp(), allvariants);
						VCFVariant v2 = q.getVariant(leadSNPInGWASdat.getSnp(), allvariants);
						DetermineLD lc = new DetermineLD();
						Pair<Double, Double> ldv = lc.getLD(v1, v2);
						int distance = topassoc.getSnp().getStart() - leadSNPInGWASdat.getSnp().getStart();
						out += "\t" + leadSNPInGWASdat.getSnp().toString()
//								+ "\t" + leadSNPInGWASdat.getSnp().getAFControls()
//								+ "\t" + leadSNPInGWASdat.getSnp().getAFCases()
//								+ "\t" + leadSNPInGWASdat.getSnp().getImputationQualityScore()
								+ "\t" + Strings.concat(leadSNPInGWASdat.getSnp().getAlleles(), Strings.comma)
								+ "\t" + leadSNPInGWASdat.getSnp().getMinorAllele()
								+ "\t" + leadSNPInGWASdat.getORs()[0][0]
								+ "\t" + leadSNPInGWASdat.getLog10Pval()
								+ "\t" + distance
								+ "\t" + ldv.getLeft() // dprime
								+ "\t" + ldv.getRight(); // rsq
					} else {
						out += "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-";
					}
					/*
					snp
					alleles
					minor
					or
					pval
					 */
					if (topassocInGWAS != null) {
						out += "\t" + topassocInGWAS.getSnp().toString()
//								+ "\t" + topassocInGWAS.getSnp().getAFControls()
//								+ "\t" + topassocInGWAS.getSnp().getAFCases()
//								+ "\t" + topassocInGWAS.getSnp().getImputationQualityScore()
								+ "\t" + Strings.concat(topassocInGWAS.getSnp().getAlleles(), Strings.comma)
								+ "\t" + topassocInGWAS.getSnp().getMinorAllele()
								+ "\t" + topassocInGWAS.getORs()[0][0]
								+ "\t" + topassocInGWAS.getLog10Pval();
					} else {
						out += "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
						;
					}
					
					/*
					leadsnpinfmdata-snp
					AfControls
					AfCases
					ImpQual
					alleles
					minor
					or
					pval
					 */
					
					if (leadSNPInFinemapData != null) {
						out += "\t" + leadSNPInFinemapData.getSnp().toString()
								+ "\t" + leadSNPInFinemapData.getSnp().getAFControls()
								+ "\t" + leadSNPInFinemapData.getSnp().getAFCases()
								+ "\t" + leadSNPInFinemapData.getSnp().getImputationQualityScore()
								+ "\t" + Strings.concat(leadSNPInFinemapData.getSnp().getAlleles(), Strings.comma)
								+ "\t" + leadSNPInFinemapData.getSnp().getMinorAllele()
								+ "\t" + leadSNPInFinemapData.getORs()[0][0]
								+ "\t" + leadSNPInFinemapData.getLog10Pval();
					} else {
						out += "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
						;
					}
					outf.writeln(out);
				}
				
				t.close();
			}
			outf.close();
			
		}
	}
	
	
	public void runNew(String[] diseasefiles, String[] gwasfiles, String[] leadsnpfiles, String[] regionfiles, String tabixPrefix, String samplefilterfile, String[] outfiles) throws IOException {
		
		
		QTLOverlap q = new QTLOverlap();
		AssociationFile af = new AssociationFile();
		BedFileReader bf = new BedFileReader();
		
		for (int d = 0; d < diseasefiles.length; d++) {
			ArrayList<Feature> regions = bf.readAsList(regionfiles[d]);
			ArrayList<AssociationResult> finemapresults = af.readRegions(diseasefiles[d], regions);
			System.out.println(finemapresults.size() + " finemap results");
			ArrayList<AssociationResult> gwasresults = af.readRegions(gwasfiles[d], regions);
			System.out.println(gwasresults.size() + " gwas results");
			Collections.sort(regions, new FeatureComparator(true));
			
			
			ArrayList<Triple<Feature, SNPFeature, SNPFeature>> regionsAndTopSNPs = new ArrayList<Triple<Feature, SNPFeature, SNPFeature>>();
			
			// load lead snps
			String leadsnpfile = leadsnpfiles[d];
			TextFile tf = new TextFile(leadsnpfile, TextFile.R);
			tf.readLine();
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				Feature region = Feature.parseFeature(elems[0]);
				String[] snpstr1 = elems[1].split("_");
				String[] snpstr2 = elems[2].split("_");
				
				SNPFeature snp1 = null;
				SNPFeature snp2 = null;
				
				if (snpstr1.length >= 3) {
					snp1 = SNPFeature.parseSNPFeature(elems[1]);
				}
				if (snpstr2.length >= 3) {
					snp2 = SNPFeature.parseSNPFeature(elems[2]);
				}
				Triple<Feature, SNPFeature, SNPFeature> f = new Triple<>(region, snp1, snp2);
				regionsAndTopSNPs.add(f);
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			
			// determine whether the top variants are in
			TextFile outf = new TextFile(outfiles[d], TextFile.W);
			String header = "Region"
					+ "\tPreviouslyReportedLeadSNP"
					+ "\tFinemap-TopAssociatedVariant"
					+ "\tAFControls"
					+ "\tAfCases"
					+ "\tImpQual"
					+ "\tAlleles"
					+ "\tMinorAllele"
					+ "\tOR"
					+ "\tPval"
					+ "\tleadSNPinGWASData-snp"
					+ "\talleles"
					+ "\tminor"
					+ "\tor"
					+ "\tpval"
					+ "\tdistance"
					+ "\tdprime"
					+ "\trsq"
					+ "\tfinemaptopsnpingwasdata-snp"
					+ "\talleles"
					+ "\tminor"
					+ "\tor"
					+ "\tpval"
					+ "\tleadsnpinfmdata-snp"
					+ "\tAfControls"
					+ "\tAfCases"
					+ "\tImpQual"
					+ "\talleles"
					+ "\tminor"
					+ "\tor"
					+ "\tpval";
			outf.writeln(header);
			for (int r = 0; r < regions.size(); r++) {
				Feature region = regions.get(r);
				String tabixFile = tabixPrefix.replaceAll("CHR", "" + region.getChromosome().getNumber());
				VCFTabix t = new VCFTabix(tabixFile);
				
				boolean[] filter = t.getSampleFilter(samplefilterfile);
				ArrayList<VCFVariant> allvariants = t.getAllVariants(region, filter);
				System.out.println(allvariants.size() + " variants in region..");
				
				ArrayList<AssociationResult> regionFinemapResults = getAssoc(finemapresults, region);
				ArrayList<AssociationResult> regionGWASResults = getAssoc(gwasresults, region);
				
				// determine lead snps in original analysis
				ArrayList<Triple<Feature, SNPFeature, SNPFeature>> regionleadsnps = getLeadSNPsForRegion(regionsAndTopSNPs, region);
				
				for (Triple<Feature, SNPFeature, SNPFeature> combo : regionleadsnps) {
					Feature gwasLeadVar = combo.getMiddle();
					Feature finemapLeadVar = combo.getRight();
					
					// find the matching assoc
					AssociationResult gwasLeadInGWASData = null;
					AssociationResult gwasLeadInFinemapData = null;
					
					
					AssociationResult finemapLeadInFinemapData = null;
					AssociationResult finemapLeadInGWASData = null;
					
					// get top variant in region
					if (finemapLeadVar == null) {
						finemapLeadVar = getTopVar(regionFinemapResults);
					}
					
					
					if (gwasLeadVar == null) {
						gwasLeadVar = getTopVar(regionGWASResults);
					}
					
					
					gwasLeadInGWASData = getVar(regionGWASResults, gwasLeadVar);
					gwasLeadInFinemapData = getVar(regionFinemapResults, gwasLeadVar);
					finemapLeadInFinemapData = getVar(regionFinemapResults, finemapLeadVar);
					finemapLeadInGWASData = getVar(regionGWASResults, finemapLeadVar);
					
					if (true) {
						// bench this!
					}
					/*
					Region
					PreviouslyReportedLeadSNP
					Finemap-TopAssociatedVariant
					AFControls
					AfCases
					ImpQual
					Alleles
					MinorAllele
					OR
					Pval
					 */
					
					
					String out = region.toString()
							+ "\t" + gwasLeadVar.toString()
							+ "\t" + finemapLeadInFinemapData.getSnp().toString()
							+ "\t" + finemapLeadInFinemapData.getSnp().getAFControls()
							+ "\t" + finemapLeadInFinemapData.getSnp().getAFCases()
							+ "\t" + finemapLeadInFinemapData.getSnp().getImputationQualityScore()
							+ "\t" + Strings.concat(finemapLeadInFinemapData.getSnp().getAlleles(), Strings.comma)
							+ "\t" + finemapLeadInFinemapData.getSnp().getMinorAllele()
							+ "\t" + finemapLeadInFinemapData.getORs()[0][0]
							+ "\t" + finemapLeadInFinemapData.getLog10Pval();
					
					
					/*
					leadSNPinGWASData-snp
					alleles
					minor
					or
					pval
					distance
					dprime
					rsq
					 */
					if (gwasLeadInGWASData != null) {
						// determine LD with topassoc
						double ld = 0;
						VCFVariant v1 = q.getVariant(finemapLeadInFinemapData.getSnp(), allvariants);
						VCFVariant v2 = q.getVariant(gwasLeadInGWASData.getSnp(), allvariants);
						DetermineLD lc = new DetermineLD();
						Pair<Double, Double> ldv = lc.getLD(v1, v2);
						int distance = finemapLeadInFinemapData.getSnp().getStart() - gwasLeadInGWASData.getSnp().getStart();
						out += "\t" + gwasLeadInGWASData.getSnp().toString()
//								+ "\t" + leadSNPInGWASdat.getSnp().getAFControls()
//								+ "\t" + leadSNPInGWASdat.getSnp().getAFCases()
//								+ "\t" + leadSNPInGWASdat.getSnp().getImputationQualityScore()
								+ "\t" + Strings.concat(gwasLeadInGWASData.getSnp().getAlleles(), Strings.comma)
								+ "\t" + gwasLeadInGWASData.getSnp().getMinorAllele()
								+ "\t" + gwasLeadInGWASData.getORs()[0][0]
								+ "\t" + gwasLeadInGWASData.getLog10Pval()
								+ "\t" + distance
								+ "\t" + ldv.getLeft() // dprime
								+ "\t" + ldv.getRight(); // rsq
					} else {
						out += "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-";
					}
					/*
					snp
					alleles
					minor
					or
					pval
					 */
					if (finemapLeadInGWASData != null) {
						out += "\t" + finemapLeadInGWASData.getSnp().toString()
//								+ "\t" + topassocInGWAS.getSnp().getAFControls()
//								+ "\t" + topassocInGWAS.getSnp().getAFCases()
//								+ "\t" + topassocInGWAS.getSnp().getImputationQualityScore()
								+ "\t" + Strings.concat(finemapLeadInGWASData.getSnp().getAlleles(), Strings.comma)
								+ "\t" + finemapLeadInGWASData.getSnp().getMinorAllele()
								+ "\t" + finemapLeadInGWASData.getORs()[0][0]
								+ "\t" + finemapLeadInGWASData.getLog10Pval();
					} else {
						out += "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
						;
					}
					
					/*
					leadsnpinfmdata-snp
					AfControls
					AfCases
					ImpQual
					alleles
					minor
					or
					pval
					 */
					
					if (gwasLeadInFinemapData != null) {
						out += "\t" + gwasLeadInFinemapData.getSnp().toString()
								+ "\t" + gwasLeadInFinemapData.getSnp().getAFControls()
								+ "\t" + gwasLeadInFinemapData.getSnp().getAFCases()
								+ "\t" + gwasLeadInFinemapData.getSnp().getImputationQualityScore()
								+ "\t" + Strings.concat(gwasLeadInFinemapData.getSnp().getAlleles(), Strings.comma)
								+ "\t" + gwasLeadInFinemapData.getSnp().getMinorAllele()
								+ "\t" + gwasLeadInFinemapData.getORs()[0][0]
								+ "\t" + gwasLeadInFinemapData.getLog10Pval();
					} else {
						out += "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
								+ "\t-"
						;
					}
					outf.writeln(out);
					
				}
				
				
				t.close();
			}
			outf.close();
			
		}
	}
	
	private AssociationResult getVar(ArrayList<AssociationResult> regionGWASResults, Feature finemapLeadVar) {
		for (AssociationResult a : regionGWASResults) {
			if (a.getSnp().overlaps(finemapLeadVar)) {
				return a;
			}
		}
		return null;
	}
	
	private SNPFeature getTopVar(ArrayList<AssociationResult> regionFinemapResults) {
		double maxp = 0;
		AssociationResult r = null;
		
		for (AssociationResult a : regionFinemapResults) {
			if (a.getLog10Pval() > maxp) {
				maxp = a.getLog10Pval();
				r = a;
			}
		}
		return r.getSnp();
	}
	
	private ArrayList<Feature> getLeadSNPs(ArrayList<Feature> leadsnps, Feature region) {
		ArrayList<Feature> output = new ArrayList<>();
		for (Feature f : leadsnps) {
			if (f.overlaps(region)) {
				output.add(f);
			}
		}
		return output;
	}
	
	private ArrayList<Triple<Feature, SNPFeature, SNPFeature>> getLeadSNPsForRegion(ArrayList<Triple<Feature, SNPFeature, SNPFeature>> leadsnps, Feature region) {
		ArrayList<Triple<Feature, SNPFeature, SNPFeature>> output = new ArrayList<>();
		for (Triple<Feature, SNPFeature, SNPFeature> f : leadsnps) {
			if (f.getLeft().overlaps(region)) {
				output.add(f);
			}
		}
		return output;
	}
	
	private ArrayList<AssociationResult> getAssoc(ArrayList<AssociationResult> results, Feature region) {
		ArrayList<AssociationResult> output = new ArrayList<>();
		for (AssociationResult r : results) {
			if (r.getSnp().overlaps(region)) {
				output.add(r);
			}
		}
		return output;
	}
	
	
}
