package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import nl.harmjanwestra.utilities.plink.PlinkFamFile;
import nl.harmjanwestra.utilities.vcf.SampleAnnotation;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.vcf.filter.genotypefilters.AllelicDepthFilter;
import nl.harmjanwestra.utilities.vcf.filter.genotypefilters.GenotypeQualityFilter;
import nl.harmjanwestra.utilities.vcf.filter.genotypefilters.ReadDepthFilter;
import nl.harmjanwestra.utilities.vcf.filter.genotypefilters.VCFGenotypeFilter;
import nl.harmjanwestra.utilities.vcf.filter.variantfilters.*;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by Harm-Jan on 02/04/16.
 */
public class VCFFilter {
	
	
	public void filter(String in,
					   String out,
					   String fam,
					   double mafthreshold,
					   double callratethreshold,
					   double hwepthreshold,
					   Double missingnessthreshold,
					   Integer readdepth,
					   Integer gqual,
					   Double allelicBalance,
					   boolean onlyAutosomes) throws IOException {
		
		TextFile tf1 = new TextFile(in, TextFile.R);
		TextFile tf2 = new TextFile(out, TextFile.W);
		String ln = tf1.readLine();
		
		
		System.out.println("Filtering: " + in);
		System.out.println("Out: " + out);
		System.out.println("FAM: " + fam);
		System.out.println("mafthreshold > " + mafthreshold);
		System.out.println("readdepth > " + readdepth);
		System.out.println("gqual > " + gqual);
		System.out.println("hwep > " + hwepthreshold);
		System.out.println("callrate > " + callratethreshold);
		System.out.println("missingness > " + missingnessthreshold);
		
		TextFile filterlog = new TextFile(out + "_log.txt", TextFile.W);
		
		VCFVariantFilters variantFilters = new VCFVariantFilters();
		variantFilters.add(new VCFVariantCallRateFilter(callratethreshold));

		if (missingnessthreshold > 0) {
			variantFilters.add(new VCFVariantMissingnessPFilter(missingnessthreshold));
		}
		
		if (fam != null) {
			variantFilters.add(new VCFVariantHWEPFilter(hwepthreshold, VCFVariantHWEPFilter.MODE.CONTROLS));
		} else if (hwepthreshold > 0) {
			variantFilters.add(new VCFVariantHWEPFilter(hwepthreshold, VCFVariantHWEPFilter.MODE.OVERALL));
		}
		
		if (fam != null) {
			variantFilters.add(new VCFVariantMAFFilter(mafthreshold, VCFVariantMAFFilter.MODE.CONTROLS));
		} else {
			variantFilters.add(new VCFVariantMAFFilter(mafthreshold, VCFVariantMAFFilter.MODE.OVERALL));
		}
		
		System.out.println(variantFilters.toString());
		if (!Gpio.exists(in)) {
			System.out.println("Could not find IN file");
			System.exit(-1);
		}
		
		SampleAnnotation sampleAnnotation = null;
		if (fam != null) {
			if (!Gpio.exists(fam)) {
				System.out.println("Could not find FAM file");
				System.exit(-1);
			}
			PlinkFamFile pf = new PlinkFamFile(fam);
			sampleAnnotation = pf.getSampleAnnotation();
		}
		
		
		ArrayList<VCFGenotypeFilter> genotypeFilters = new ArrayList<>();
		if (gqual != null) {
			genotypeFilters.add(new GenotypeQualityFilter(gqual));
			System.out.println("Adding genotyping filter: " + gqual);
		}
		if (allelicBalance != null) {
			if (readdepth == null) {
				genotypeFilters.add(new AllelicDepthFilter(allelicBalance, 0));
				System.out.println("Adding allele balance filter: ab" + allelicBalance + "\tread depth: " + 0);
			} else {
				genotypeFilters.add(new AllelicDepthFilter(allelicBalance, readdepth));
				System.out.println("Adding allele balance filter: ab" + allelicBalance + "\tread depth: " + readdepth);
			}
		} else if (readdepth != null) {
			System.out.println("Adding read depth filter: " + readdepth);
			genotypeFilters.add(new ReadDepthFilter(readdepth));
		}
		
		String logheader = "Variant"
				+ "\tAlleles"
				+ "\tMinorAllele"
				+ "\tMAFCases"
				+ "\tMAFControls"
				+ "\tMAF"
				+ "\tHWEPCases"
				+ "\tHWEPControls"
				+ "\tHWEP"
				+ "\tCallRateCases"
				+ "\tCallRateControls"
				+ "\tCallRate"
				+ "\tCallRate-P"
				+ "\tPassesThresholds"
				+ "\tFailedFilter";
		
		System.out.println("Starting filter...");
		filterlog.writeln(logheader);
		int lnctr = 0;
		int kept = 0;
		while (ln != null) {
			if (ln.startsWith("##")) {
				tf2.writeln(ln);
			} else if (ln.startsWith("#CHROM")) {
				tf2.writeln("##VCFGenotypeFilter=maf>" + mafthreshold + ",cr>" + callratethreshold + ",dp>" + readdepth + ",gq>" + gqual + ",ab>" + allelicBalance + ",missingnessp>" + missingnessthreshold);
				
				// reorder sample annotation
				ArrayList<String> samples = new ArrayList<>();
				String[] elems = ln.split("\t");
				for (int i = 9; i < elems.length; i++) {
					samples.add(elems[i]);
				}
				
				System.out.println(samples.size() + " individuals in VCF");
				if (sampleAnnotation != null) {
					sampleAnnotation.reorder(samples, true);
					int nrSamplesWithAnnotation = sampleAnnotation.getSampleName().length;
					if (nrSamplesWithAnnotation != samples.size()) {
						System.err.println("Error: nr of samples with annotation: " + nrSamplesWithAnnotation + ", while number of samples in VCF: " + samples.size());
						System.exit(-1);
					}
				}
				tf2.writeln(ln);
			} else {
				boolean passesThresholds = false;
				VCFVariant var = new VCFVariant(ln, genotypeFilters, true, sampleAnnotation);
				
				boolean autosomal = var.getChrObj().isAutosome();
				if (autosomal) {
					if (variantFilters.passesFilters(var)) {
						tf2.writeln(ln);
						kept++;
						passesThresholds = true;
					}
				}

//				if (var.getId().equals("rs3826110")) {
//					VCFVariantMissingnessPFilter v = new VCFVariantMissingnessPFilter(missingnessthreshold);
//					System.out.println(var.getId() + "\t" + var.getDiffMissingnessP() + "\t" + v.passesThreshold(var) + "\t" + passesThresholds);
//				}
				
				String logout = var.toString()
						+ "\t" + Strings.concat(var.getAlleles(), Strings.comma)
						+ "\t" + var.getMinorAllele()
						+ "\t" + var.getMAFCases()
						+ "\t" + var.getMAFControls()
						+ "\t" + var.getMAF()
						+ "\t" + var.getHwepCases()
						+ "\t" + var.getHwepControls()
						+ "\t" + var.getHwep()
						+ "\t" + var.getCallrateCases()
						+ "\t" + var.getCallrateControls()
						+ "\t" + var.getCallrate()
						+ "\t" + var.getDiffMissingnessP()
						+ "\t" + passesThresholds;
				if (variantFilters.getFailedFilter() != null) {
					logout += "\t" + variantFilters.getFailedFilter().toString();
				} else {
					logout += "\tNA";
				}
				filterlog.writeln(logout);
				
			}
			ln = tf1.readLine();
			if (lnctr % 1000 == 0) {
				System.out.print("\r"+lnctr + " lines parsed. " + kept + " kept.");
			}
			lnctr++;
		}
		System.out.println();
		
		tf1.close();
		tf2.close();
		
		filterlog.close();
		System.out.println("Done");
		System.out.println(lnctr + " lines parsed. " + kept + " kept.");
	}

	/*

		// Maybe this filtering should be moved to a separate class or
		// something?
		if (minimalReadDepth > 0 && dpCol != -1) {

		}
		if (minimalGenotypeQual > 0 && gqCol != -1) {

		}
	 */
	
	public void filterNonACTGVariants(String vcfin, String vcfout) throws IOException {
		System.out.println("in: " + vcfin);
		System.out.println("out: " + vcfout);
		TextFile out = new TextFile(vcfout, TextFile.W);
		TextFile tf = new TextFile(vcfin, TextFile.R);
		
		int written = 0;
		int total = 0;
		String ln = tf.readLine();
		while (ln != null) {
			if (ln.startsWith("#")) {
				out.writeln(ln);
			} else {
				VCFVariant var = new VCFVariant(ln, VCFVariant.PARSE.HEADER);
				String[] alleles = var.getAlleles();
				boolean allelesok = true;
				for (String allele : alleles) {
					if (allele.toUpperCase().contains("A") ||
							allele.toUpperCase().contains("C") ||
							allele.toUpperCase().contains("T") ||
							allele.toUpperCase().contains("G")) {
					} else {
						allelesok = false;
					}
				}
				if (allelesok) {
					out.writeln(ln);
					written++;
				} else {
					System.out.println("Weird alleles: " + var.getId() + "\t" + Strings.concat(alleles, Strings.forwardslash));
				}
				total++;
			}
			ln = tf.readLine();
		}
		
		System.out.println("weird alleles: " + written + "/" + total + " written..");
		tf.close();
		out.close();
	}
}
