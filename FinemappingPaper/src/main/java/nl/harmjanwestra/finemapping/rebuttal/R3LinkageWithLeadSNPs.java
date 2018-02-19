package nl.harmjanwestra.finemapping.rebuttal;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
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
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz"
		};
		String[] gwasfiles = new String[]{
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Okada\\RA.OKADA.gz",
				disk + "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/ImmunoBase/hg19_gwas_ic_t1d_onengut_meta_4_19_1.tab.gz"
		};
		String[] leadsnpfiles = new String[]{
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\leadsnps\\annot_2013_ORs.txt",
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\leadsnps\\t1d.txt"
		};
		String[] regionfiles = new String[]{
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-significantregions-75e7.bed",
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-significantregions-75e7.bed"
		};
		String tabix = disk + "\\Data\\Ref\\1kg\\ALL.chrCHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";
		String samplefile = disk + "\\Data\\Ref\\1kg-europeanpopulations.txt.gz";
		String[] outfiles = new String[]{
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\comparisontoimmunobase2\\ra.txt",
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\comparisontoimmunobase2\\t1d.txt"
		};
		
		R3LinkageWithLeadSNPs r = new R3LinkageWithLeadSNPs();
		try {
			r.run(diseasefiles, gwasfiles, leadsnpfiles, regionfiles, tabix, samplefile, outfiles);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void run(String[] diseasefiles, String[] gwasfiles, String[] leadsnpfiles, String[] regionfiles, String tabixPrefix, String samplefilterfile, String[] outfiles) throws IOException {
		
		
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
	
	private ArrayList<Feature> getLeadSNPs(ArrayList<Feature> leadsnps, Feature region) {
		ArrayList<Feature> output = new ArrayList<>();
		for (Feature f : leadsnps) {
			if (f.overlaps(region)) {
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
