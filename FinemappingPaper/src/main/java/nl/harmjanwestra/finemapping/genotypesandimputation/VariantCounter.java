package nl.harmjanwestra.finemapping.genotypesandimputation;

import nl.harmjanwestra.finemapping.rebuttal.KgVariant;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Triple;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 5/28/16.
 */
public class VariantCounter {
	
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
	
	public static void main(String[] args) {
		try {
			VariantCounter c = new VariantCounter();
			c.countAccuracy();
//			c.countINFO();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public void determineMissingVariants() throws IOException {
		
		String[] files = new String[]{
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-07-10-Accuracy/T1D-COSMO.txt",
				"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/2016-07-10-Accuracy/RA-COSMO.txt"
		};
		
		String variantsOnIC = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-06-21-ImputationQuality/RAAndT1D-recode-maf0005-ICRegionsW100kb-samplenamefix.vcf.gz-updatedRSId-stats.vcf.gz";
		String seqpanelvcf = "/Data/tmp/2016-05-28/seqpanelfiltered-maf0005-cr0950-rd10-gq30-runNamesFixed-RASampleNamesFixed-badSamplesRemoved-mixupsFixed.vcf.gz";
		String bedregions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		String samplelist = "";
		
		double mafthreshold = 0.01;
		double upperthreshold = 1;
		double infothreshold = 0.5;
		
		boolean includeICVariants = false;
		boolean includeId = true;
		boolean includeIndels = true;
		
		
		HashSet<String> variantsOnICHash = loadVariantHash(variantsOnIC, includeId);
		System.out.println(variantsOnICHash.size() + " total on IC");
		
		
		Triple<ArrayList<VCFVariant>, ArrayList<VCFVariant>, ArrayList<VCFVariant>> seqpanelvariants = loadSequencedVariants(
				seqpanelvcf, bedregions, mafthreshold, upperthreshold, variantsOnICHash, includeId, includeIndels, samplelist
		);
		
		ArrayList<VCFVariant> seqpanel = seqpanelvariants.getLeft();
		ArrayList<VCFVariant> variantsOnImmunoChip = seqpanelvariants.getMiddle();
		ArrayList<VCFVariant> variantsNotOnImmunoChip = seqpanelvariants.getRight();
		
		System.out.println(seqpanel.size() + " variants in VCF");
		System.out.println(variantsNotOnImmunoChip.size() + " not on IC");
		System.out.println(variantsOnImmunoChip.size() + " on IC");
		
		
		// get a list of imputed variants for each of the sequencing panels
		
		
		System.out.println("MAF> " + mafthreshold);
		System.out.println("INFO> " + infothreshold);
		
		
		HashSet<String> sequencedVariantsHash = null;
		if (includeICVariants) {
			ArrayList<VCFVariant> allvars = new ArrayList<>();
			allvars.addAll(variantsNotOnImmunoChip);
			allvars.addAll(variantsOnImmunoChip);
			sequencedVariantsHash = hashit(allvars, includeId);
		} else {
			sequencedVariantsHash = hashit(variantsNotOnImmunoChip, includeId);
		}
		
		System.out.println(sequencedVariantsHash.size() + " after hashing");
		System.out.println(files.length + " files");
		HashSet<String> variantsFound = new HashSet<String>();
		
		for (int f = 0; f < files.length; f++) {
			// get the imputation accuracies for these variants
			TextFile tf2 = new TextFile(files[f], TextFile.R);
			
			String[] elems = tf2.readLineElems(TextFile.tab);
			int nrSequenced = 0;
			int nrSequencedPassingRSQ = 0;
			int nrSequencedPassingMaf = 0;
			int nrSequencdPassingMafAndRSQ = 0;
			while (elems != null) {
				if (!elems[rsqlcol].equals("null")) {
					double val = Double.parseDouble(elems[rsqlcol]);
					double maf = Double.parseDouble(elems[maf2col]);
					
					String[] varElems = elems[0].split("_");
					
					boolean sequenced = isVariantInHash(varElems, sequencedVariantsHash, includeId);
					
					if (sequenced) {
						String variant = Chromosome.parseChr(varElems[0]).toString() + "_" + varElems[1] + "_" + varElems[2];
//						System.out.println(variant);
						variantsFound.add(variant);
					}
					
					if (sequenced) {
						nrSequenced++;
						if (maf > mafthreshold) {
							nrSequencedPassingMaf++;
							if (val > infothreshold) {
								nrSequencdPassingMafAndRSQ++;
							}
						}
						if (val > infothreshold) {
							nrSequencedPassingRSQ++;
						}
					}
				}
				
				elems = tf2.readLineElems(TextFile.tab);
			}
			tf2.close();

//			System.out.println("----");
//			System.out.println();

//			System.out.println("nrSequenced\t" + nrSequenced + "\t" + ((double) nrSequenced / sequencedVariantsHash.size()));
//			System.out.println("nrSequencedPassingRSQ\t" + nrSequencedPassingRSQ + "\t" + ((double) nrSequencedPassingRSQ / sequencedVariantsHash.size()) + "\t" + ((double) nrSequencedPassingRSQ / nrSequenced));
//			System.out.println("nrSequencedPassingMaf\t" + nrSequencedPassingMaf + "\t" + ((double) nrSequencedPassingMaf / sequencedVariantsHash.size()));
//			System.out.println("nrSequencdPassingMafAndRSQ\t" + nrSequencdPassingMafAndRSQ + "\t" + ((double) nrSequencdPassingMafAndRSQ / sequencedVariantsHash.size()));
//			System.out.println();
			
			System.out.println(nrSequenced + "\t" + nrSequencedPassingMaf + "\t" + nrSequencedPassingRSQ + "\t" + nrSequencdPassingMafAndRSQ);
		}
		
		System.out.println(variantsFound.size() + " total variants found.");
		int ctr = 0;
		for (String var : sequencedVariantsHash) {
			if (!variantsFound.contains(var)) {
				System.out.println(var);
				ctr++;
			}
		}
		
		System.out.println(ctr + " variants missing");
		
	}
	
	enum SNPTTYPE {
		INDEL, MULTI, SNP
	}
	
	public void countAccuracy() throws IOException {
		String disk = "d:";
		String[] variantsOnIC = new String[]{
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2016-06-21-ImputationQuality\\RAAndT1D-recode-maf0005-ICRegionsW100kb-samplenamefix.vcf.gz-updatedRSId-stats.vcf.gz"
		};
		
		String[] seqpanelvcfs = new String[]{
//				"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\panels\\unifiedgenotyper-maf0005-cr0950-rd10-gq30.vcf.gz-samplenamefix-mixupfix-nonmatchingremoved.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\panels\\hapcaller-maf0005-cr0950-rd10-gq30.vcf.gz-samplenamefix-mixupfix-nonmatchingremoved.vcf.gz"
		};
		
		String[] seqpanelnames = new String[]{
//				"UnifiedGenotyper"
				"HaplotypeCaller"
		};
		
		
		String bedregions = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		bedregions = disk + "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		String[] files = new String[]{
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
		String[] labels = new String[]{
				"RA / COSMO / BEAGLE",
				"RA / COSMO / EAGLE / PBWT",
				"RA / EUR / BEAGLE",
				"RA / HRC / EAGLE / PBWT",
				"T1D / COSMO / BEAGLE",
				"T1D / COSMO / EAGLE / PBWT",
				"T1D / EUR / BEAGLE",
				"T1D / HRC / EAGLE / PBWT",
				"T1D / HRC / EAGLE / SHAPEIT / PBWT",
				"T1D / HRC / EAGLE / MACH"
		};
		
		String outfile = disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2-VariantCounter\\output";
		String samplelist = null;
		
		// get a list of maf > 0.005 variants that are on the sequencingpanel
		double mafthreshold = 0.01;
		double upperthreshold = 1;
		double infothreshold = 0.5;
		boolean includeICVariants = true;
		boolean includeAlleles = false;
		boolean includeId = false;
		boolean includeIndels = true;
		
		
		for (int r = 0; r < seqpanelvcfs.length; r++) {
			HashSet<String> variantsOnICHash = loadVariantHash(variantsOnIC[0], includeId);
			System.out.println(variantsOnICHash.size() + " total on IC");
			
			Triple<ArrayList<VCFVariant>, ArrayList<VCFVariant>, ArrayList<VCFVariant>> seqpanelvariants = loadSequencedVariants(
					seqpanelvcfs[r], bedregions, mafthreshold, upperthreshold, variantsOnICHash, includeId, includeIndels, samplelist
			);
			
			ArrayList<VCFVariant> seqpanel = seqpanelvariants.getLeft();
			TextFile out = new TextFile(outfile, TextFile.W);
			out.writeln("Chromosome\tPosition\tRsid\tRef\tAlt\tAF\tCallRate\tHWE-P");
			for (VCFVariant var : seqpanel) {
				out.writeln(var.getChr().toString()
						+ "\t" + var.getPos()
						+ "\t" + var.getId()
						+ "\t" + var.getAlleles()[0]
						+ "\t" + Strings.concat(var.getAlleles(), Strings.comma, 1, var.getAlleles().length)
						+ "\t" + Strings.concat(var.getAlleleFrequencies(), Strings.semicolon)
						+ "\t" + var.getCallrate()
						+ "\t" + var.getHwep()
				);
			}
			out.close();
//				System.exit(-1);
			
			ArrayList<VCFVariant> variantsOnImmunoChip = seqpanelvariants.getMiddle();
			ArrayList<VCFVariant> variantsNotOnImmunoChip = seqpanelvariants.getRight();
			
			System.out.println(seqpanel.size() + " variants in VCF");
			System.out.println(variantsNotOnImmunoChip.size() + " not on IC");
			System.out.println(variantsOnImmunoChip.size() + " on IC");
			
			// get a list of imputed variants for each of the sequencing panels
			System.out.println("MAF> " + mafthreshold);
			System.out.println("INFO> " + infothreshold);
			
			HashSet<String> sequencedVariantsHash = null;
			HashMap<String, SNPTTYPE> types = null;
			
			if (includeICVariants) {
				ArrayList<VCFVariant> allvars = new ArrayList<>();
				allvars.addAll(variantsNotOnImmunoChip);
				allvars.addAll(variantsOnImmunoChip);
				
				sequencedVariantsHash = hashit(allvars, includeId);
			} else {
				sequencedVariantsHash = hashit(variantsNotOnImmunoChip, includeId);
			}


//			HashSet<KgVariant> allRefVariantsSet = new HashSet<>();
//			ArrayList<KgVariant> refVariants = new ArrayList<>();
//			for (int i = 0; i < seqpanel.size(); i++) {
//				VCFVariant v = seqpanel.get(i);
//				KgVariant var = new KgVariant();
//				var.f = v.asSNPFeature();
//
//				var.maf = v.getMAF();
//				var.hwep = v.getHwep();
//
//				var.f.useNameForComparison(false);
//				if (var.f.isMultiAllelic()) {
//					var.f.useAllelesForComparison(false);
//				}
//				if (!includeAlleles) {
//					var.f.useAllelesForComparison(false);
//				}
//				if (!includeId) {
//					var.f.setName(null);
//				}
//
//				if (v.isIndel() && includeIndels || !v.isIndel()) {
//					refVariants.add(var);
//					allRefVariantsSet.add(var);
//				}
//			}
			
			System.out.println(sequencedVariantsHash.size() + " after hashing");
			System.out.println(files.length + " files");
			for (int f = 0; f < files.length; f++) {
				// get the imputation accuracies for these variants
				TextFile tf2 = new TextFile(files[f], TextFile.R);
				
				String[] elems = tf2.readLineElems(TextFile.tab);
				int nrSequenced = 0;
				int nrSequencedPassingRSQ = 0;
				int nrSequencedPassingMaf = 0;
				int nrSequencdPassingMafAndRSQ = 0;
				
				// determine which variants are missing and what type they are...
				int multiAllelicOverlap = 0;
				int indeloverlap = 0;
				int snpoverlap = 0;
				
				int multiAllelicMissing = 0;
				int indeloverlapMissing = 0;
				int snpoverlapMissing = 0;
				
				TextFile outf = new TextFile("D:\\Data\\tmp\\23018\\overlap" + f + ".txt", TextFile.W);
				while (elems != null) {
					if (elems.length > 1) {
						String[] varElems = elems[0].split("_");
						boolean sequenced = isVariantInHash(varElems, sequencedVariantsHash, includeId);
						if (sequenced) {
							nrSequenced++;
							outf.writeln(elems[0]);
						}
						
						if (!elems[rsqlcol].equals("null")) {
							double val = Double.parseDouble(elems[rsqlcol]);
							double maf = Double.parseDouble(elems[maf2col]);
							
							if (maf > mafthreshold) {
								nrSequencedPassingMaf++;
								if (val > infothreshold) {
									nrSequencdPassingMafAndRSQ++;
								}
							}
							if (val > infothreshold) {
								nrSequencedPassingRSQ++;
							}
						}
					}
					elems = tf2.readLineElems(TextFile.tab);
				}
				outf.close();
				tf2.close();
				
				System.out.println(seqpanelnames[r] + "\t" + labels[f]
						+ "\t" + nrSequenced
						+ "\t" + nrSequencedPassingMaf
						+ "\t" + nrSequencedPassingRSQ
						+ "\t" + nrSequencdPassingMafAndRSQ);
			}
		}
		
		
	}
	
	public void countINFO() throws IOException {
		String disk = "d:";
		String[] variantsOnIC = new String[]{
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2016-06-21-ImputationQuality\\RAAndT1D-recode-maf0005-ICRegionsW100kb-samplenamefix.vcf.gz-updatedRSId-stats.vcf.gz"
		};
		
		String[] seqpanelvcfs = new String[]{
//				"c:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\panels\\unifiedgenotyper-maf0005-cr0950-rd10-gq30.vcf.gz-samplenamefix-mixupfix-nonmatchingremoved.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2\\panels\\hapcaller-maf0005-cr0950-rd10-gq30.vcf.gz-samplenamefix-mixupfix-nonmatchingremoved.vcf.gz"
		};
		
		String[] seqpanelnames = new String[]{
//				"UnifiedGenotyper"
				"HaplotypeCaller"
		};
		
		
		String bedregions = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		bedregions = disk + "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		
		
		String[] files = new String[]{
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
		String[] labels = new String[]{
				"RA / COSMO / BEAGLE",
				"RA / COSMO / EAGLE / PBWT",
				"RA / EUR / BEAGLE",
				"RA / HRC / EAGLE / PBWT",
				"T1D / COSMO / BEAGLE",
				"T1D / COSMO / EAGLE / PBWT",
				"T1D / EUR / BEAGLE",
				"T1D / HRC / EAGLE / PBWT",
				"T1D / HRC / SHAPEIT / PBWT",
				"T1D / HRC / EAGLE / MACH"
		};
		
		String[] rawinputs = new String[]{
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\RA-RAW.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\RA-RAW.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\RA-RAW.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\RA-RAW.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\T1D-RAW.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\T1D-RAW.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\T1D-RAW.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\T1D-RAW.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\T1D-RAW.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\T1D-RAW.vcf.gz"
			
		};
		String[] matchedinputs = new String[]{
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\RA-COSMO.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\RA-COSMO.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\RA-EUR.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\RA-HRC.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\T1D-COSMO.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\T1D-COSMO.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\T1D-EUR.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\T1D-HRC.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\T1D-HRC.vcf.gz",
				disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\T1D-HRC.vcf.gz",
			
			
		};
		
		
		String outfile = disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\accuracy\\R2-VariantCounter\\outputINFO";
		String samplelist = null;
		
		// get a list of maf > 0.005 variants that are on the sequencingpanel
		double mafthreshold = 0.01;
		double upperthreshold = 1;
		double infothreshold = 0.3;
		boolean includeICVariants = true;
		boolean includeAlleles = false;
		boolean includeId = true;
		boolean includeIndels = true;
		
		
		BedFileReader r = new BedFileReader();
		ArrayList<Feature> regions = r.readAsList(bedregions);
		
		
		System.out.println(files.length + " files");
		for (int f = 0; f < files.length; f++) {
			
			
			String rawinput = rawinputs[f];
			String matchedinput = matchedinputs[f];
			
			// count input
			int nrinput = 0;
			TextFile tf = new TextFile(rawinput, TextFile.R);
			String ln1 = tf.readLine();
			while (ln1 != null) {
				if (!ln1.startsWith("#")) {
					VCFVariant v = new VCFVariant(ln1);
					if (v.asFeature().overlaps(regions)) {
						nrinput++;
						
					}
				}
				ln1 = tf.readLine();
			}
			tf.close();
			
			
			int nrinputmaf005 = 0;
			int nrinputmaf010 = 0;
			TextFile tf3 = new TextFile(matchedinput, TextFile.R);
			String ln3 = tf3.readLine();
			while (ln3 != null) {
				if (!ln3.startsWith("#")) {
					VCFVariant v = new VCFVariant(ln3);
					if (v.asFeature().overlaps(regions)) {
						double af = Double.parseDouble(v.getInfo().get("AF"));
						if (af > 0.5) {
							af = 1 - af;
						}
						if (af > 0.005) {
							nrinputmaf005++;
						}
						if (af > 0.01) {
							nrinputmaf010++;
						}
					}
				}
				ln3 = tf3.readLine();
			}
			tf3.close();
			
			
			// get the imputation accuracies for these variants
			TextFile tf2 = new TextFile(files[f], TextFile.R);
			
			String ln2 = tf2.readLine();
			int nrSequenced = 0;
			int nrSequencedPassingRSQ = 0;
			int nrSequencedPassingMaf = 0;
			int nrSequencdPassingMafAndRSQ = 0;
			int nrSequencedNotPassingRSQ = 0;
			int nrSequencdPassingMafAndRSQNonIndel = 0;
			while (ln2 != null) {
				if (!ln2.startsWith("#")) {
					VCFVariant v = new VCFVariant(ln2);
					if (v.asFeature().overlaps(regions)) {
						String afstr = v.getInfo().get("AF");
						String[] split = afstr.split(",");
						double lowest = 1;
						for (int d = 0; d < split.length; d++) {
							Double af = Double.parseDouble(split[d]);
							if (af < lowest) {
								lowest = af;
							}
						}
						double af = lowest;
						if (af > 0.5) {
							af = 1 - af;
						}
						nrSequenced++;
						double info = v.getImputationQualityScore();
						if (af > mafthreshold) {
							nrSequencedPassingMaf++;
						}
						if (info >= infothreshold) {
							nrSequencedPassingRSQ++;
						} else {
							nrSequencedNotPassingRSQ++;
						}
						
						
						if (info > infothreshold) {
							if (af > mafthreshold) {
								nrSequencdPassingMafAndRSQ++;
								if (!v.isIndel()) {
									nrSequencdPassingMafAndRSQNonIndel++;
								}
							}
							
						}
					}
					
					
				}
				ln2 = tf2.readLine();
			}
			tf2.close();
			
			System.out.println(labels[f]
					+ "\t" + nrinput
					+ "\t" + nrinputmaf005
					+ "\t" + nrinputmaf010
					+ "\t" + nrSequenced
					+ "\t" + nrSequencedPassingRSQ
					+ "\t" + nrSequencedPassingMaf
					+ "\t" + nrSequencdPassingMafAndRSQ
					+ "\t" + nrSequencdPassingMafAndRSQNonIndel);
//			System.out.println(nrSequencedNotPassingRSQ);
		}
		
		
	}
	
	
	public HashSet<String> loadVariantHash(String variantsOnIC, boolean includeId) throws IOException {
		TextFile tf = new TextFile(variantsOnIC, TextFile.R);
		HashSet<String> variantIds = new HashSet<String>();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (!elems[0].startsWith("#")) {
				if (elems.length > 2) {
					if (includeId) {
						String id = Chromosome.parseChr(elems[0]).toString() + "_" + elems[1] + "_" + elems[2];
						variantIds.add(id);
					} else {
						String id = Chromosome.parseChr(elems[0]).toString() + "_" + elems[1];
						variantIds.add(id);
					}
				}
				
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		return variantIds;
	}
	
	public HashSet<String> hashit(ArrayList<VCFVariant> variantsNotOnImmunoChip, boolean includeId) {
		HashSet<String> variantHash = new HashSet<String>();
		for (VCFVariant var : variantsNotOnImmunoChip) {
			String variant = "";
			if (includeId) {
				variant = var.getChrObj().toString() + "_" + var.getPos() + "_" + var.getId();
			} else {
				variant = var.getChrObj().toString() + "_" + var.getPos();
			}
			variantHash.add(variant);
		}
		
		
		return variantHash;
	}
	
	public boolean isIndel(String[] elems) {
		
		String alleles = elems[3] + "," + elems[4];
		String[] alleleElems = alleles.split(",");
		
		for (String s : alleleElems) {
			if (s.length() > 1) {
				return true;
			}
		}
		return false;
	}
	
	public Triple<ArrayList<VCFVariant>, ArrayList<VCFVariant>, ArrayList<VCFVariant>> loadSequencedVariants(String seqpanelvcf,
																											 String bedregionsFile,
																											 double mafthreshold,
																											 double upperthreshold,
																											 HashSet<String> variantsOnICHash,
																											 boolean includeId,
																											 boolean includeIndels,
																											 String sampleList) throws IOException {
		
		
		boolean[] includesamples = null;
		if (sampleList != null) {
			VCFGenotypeData d = new VCFGenotypeData(seqpanelvcf);
			ArrayList<String> allsamples = d.getSamples();
			d.close();
			
			TextFile tf = new TextFile(sampleList, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			HashSet<String> sampleSet = new HashSet<String>();
			while (elems != null) {
				sampleSet.add(elems[0]);
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			
			includesamples = new boolean[allsamples.size()];
			int ctr = 0;
			for (int i = 0; i < allsamples.size(); i++) {
				if (sampleSet.contains(allsamples.get(i))) {
					includesamples[i] = true;
					ctr++;
				}
			}
			System.out.println(ctr + " samples overlapping with: " + sampleList);
		}
		
		
		BedFileReader bfr = new BedFileReader();
		ArrayList<Feature> regions = bfr.readAsList(bedregionsFile);
		
		TextFile tf = new TextFile(seqpanelvcf, TextFile.R);
		String ln = tf.readLine();
		
		ArrayList<VCFVariant> variantsNotOnImmunoChip = new ArrayList<>();
		ArrayList<VCFVariant> variantsOnImmunoChip = new ArrayList<>();
		ArrayList<VCFVariant> seqpanel = new ArrayList<>();
		
		while (ln != null) {
			if (!ln.startsWith("#")) {
				String[] elems = ln.split("\t");
//				Chromosome chr = Chromosome.parseChr(elems[0]);
//				if (chr.isAutosome()) {
				
				VCFVariant variant = new VCFVariant(ln, VCFVariant.PARSE.ALL);
//					if (variant.getMAF() > mafthreshold && variant.getMAF() < upperthreshold) {
				if (variant.asFeature().getChromosome().isAutosome() && variant.getMAF() > mafthreshold && variant.asFeature().overlaps(regions)) {
					
					boolean varOnIc = isVariantInHash(elems, variantsOnICHash, includeId);
					boolean iswithinregion = isWithinRegion(regions, elems);
					boolean indel = isIndel(elems);
//						if (iswithinregion) {
					seqpanel.add(variant);
					if (includeIndels || (!includeIndels && !indel)) {
						if (!varOnIc) {
							variantsNotOnImmunoChip.add(variant);
						} else {
							variantsOnImmunoChip.add(variant);
						}
					}
//						}
				}
//				}
			}
			ln = tf.readLine();
		}
		tf.close();
		return new Triple<>(seqpanel, variantsOnImmunoChip, variantsNotOnImmunoChip);
	}
	
	public boolean isVariantInHash(String[] elems, HashSet<String> variantHash, boolean includeId) {
		if (includeId) {
			String variant = Chromosome.parseChr(elems[0]).toString() + "_" + elems[1] + "_" + elems[2];
			return variantHash.contains(variant);
		} else {
			String variant = Chromosome.parseChr(elems[0]).toString() + "_" + elems[1];
			return variantHash.contains(variant);
		}
		
	}
	
	public boolean mafbelowfthreshold(String[] elems, double t) {
		double maf = getMaf(elems);
		if (maf < t) {
			return true;
		}
		return false;
	}
	
	public boolean isWithinRegion(ArrayList<Feature> list, String[] elems) {
		Feature varfeat = new Feature();
		
		int pos = Integer.parseInt(elems[1]);
		varfeat.setChromosome(Chromosome.parseChr(elems[0]));
		varfeat.setStart(pos);
		varfeat.setStop(pos);
		for (Feature f : list) {
			if (f.overlaps(varfeat)) {
				return true;
			}
		}
		return false;
	}
	
	public double getMaf(String[] elems) {
		String[] infoElems = elems[7].split(";");
		double maf = 1;
		for (String info : infoElems) {
			if (info.startsWith("AF")) {
				String[] elems2 = info.split("=");
				String[] elems3 = elems2[1].split(",");
				for (String afStr : elems3) {
					double af = Double.parseDouble(afStr);
					if (af > 0.5) {
						af = 1 - af;
					}
					if (af < maf) {
						maf = af;
					}
					
				}
			}
		}
		return maf;
	}
	
	public double getInfo(String[] elems) {
		String[] infoElems = elems[7].split(";");
		double infoscore = 0;
		for (String info : infoElems) {
			if (info.startsWith("INFO")) {
				String[] elems2 = info.split("=");
				infoscore = Double.parseDouble(elems2[1]);
			}
		}
		return infoscore;
	}
	
	public boolean isIndel1(String[] elems) {
		String alleles = elems[minorallele1col];
		String[] alleleElems = alleles.split(",");
		boolean b = false;
		for (String s : alleleElems) {
			if (s.length() > 1) {
				b = true;
			}
		}
		return b;
	}
	
	public boolean isIndel2(String[] elems) {
		String alleles = elems[minorallele2col];
		String[] alleleElems = alleles.split(",");
		boolean b = false;
		for (String s : alleleElems) {
			if (s.length() > 1) {
				b = true;
			}
		}
		return b;
	}
	
	public boolean maf1belowfthreshold(String[] elems, double t) {
		Double d = Double.parseDouble(elems[maf1col]);
		if (d < t) {
			return true;
		} else {
			return false;
		}
	}
	
	public boolean maf2belowfthreshold(String[] elems, double t) {
		Double d = Double.parseDouble(elems[maf2col]);
		if (d < t) {
			return true;
		} else {
			return false;
		}
	}
	
	
}
