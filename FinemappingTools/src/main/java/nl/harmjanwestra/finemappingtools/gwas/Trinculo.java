package nl.harmjanwestra.finemappingtools.gwas;

import cern.colt.matrix.tdouble.DoubleMatrix2D;

import nl.harmjanwestra.finemappingtools.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.DiseaseStatus;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.console.ProgressBar;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.vcf.VCFVariantComparator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.*;

public class Trinculo extends LRTest {
	
	
	public Trinculo(LRTestOptions options) throws IOException {
		super(options);
		run();
	}
	
	private void run() throws IOException {
		
		// trinculo multinom --dosage genotypes --pheno phenos.txt --phenoname Pheno --basepheno Control --covar pcs.txt --normalize --defaultprior --out example5
		
		
		// phenos.txt - space sep
		/*
		FID IID Pheno	Order
		IND1 IND1 Disease1	1
		IND2 IND2 Control	0
		 */
		
		/* pcs.txt
		FID IID PC1 PC2
		IND1 IND1 0.08411336 -0.752277
		IND2 IND2 0.7129942 -0.9827497
		IND3 IND3 -0.2415766 0.8311323
		 */
		System.out.println("Running TRINCULO converter...");
		options.setSplitMultiAllelic(true);
		
		DoubleMatrix2D covariates = sampleAnnotation.getCovariates();
		
		String[] samplenames = sampleAnnotation.getSampleName();
		TextFile out = new TextFile(options.getOutputdir() + "pcs.txt", TextFile.W);
		String header = "FID IID";
		for (int c = 0; c < covariates.columns(); c++) {
			header += " cov" + c;
		}
		out.writeln(header);
		
		for (int s = 0; s < samplenames.length; s++) {
			String lnout = samplenames[s] + " " + samplenames[s];
			for (int c = 0; c < covariates.columns(); c++) {
				lnout += " " + covariates.get(s, c);
			}
			out.writeln(lnout);
		}
		out.close();
		
		DiseaseStatus[][] statuses = sampleAnnotation.getSampleDiseaseStatus();
		System.out.println(statuses[0].length + " diseases loaded.");
		TextFile out2 = new TextFile(options.getOutputdir() + "pheno.txt", TextFile.W);
		String header2 = "FID IID Pheno\tOrder";
		out2.writeln(header2);
		for (int s = 0; s < samplenames.length; s++) {
			String lnout = samplenames[s] + " " + samplenames[s];
			
			boolean isCase = false;
			int disease = 0;
			
			for (int d = 0; d < statuses[s].length; d++) {
				DiseaseStatus stat = statuses[s][d];
				if (stat.equals(DiseaseStatus.CASE)) {
					disease = d;
					isCase = true;
				}
			}
			if (!isCase) {
				lnout += " ControlPheno\t0";
			} else {
				lnout += " Disease" + (disease + 1) + "\t" + (disease + 1);
			}
			out2.writeln(lnout);
		}
		out2.close();
		
		// read the association file
		if (options.getAssocFile() == null) {
			System.err.println("Please supply assoc file.");
			System.exit(-1);
		}
		
		AssociationFile f = new AssociationFile();
		String assocfile = options.getAssocFile();
		ArrayList<AssociationResult> assocResults = f.read(assocfile);
		
		System.out.println(assocResults.size() + " association results loaded from: " + assocfile);
		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(options.getBedfile());
		ArrayList<VCFVariant> variants = readVariants(options.getOutputdir() + "variantlog.txt", regions);
		Collections.sort(variants, new VCFVariantComparator());
		
		if (variants.isEmpty()) {
			System.out.println("Sorry no variants found.");
			System.exit(-1);
		}
		
		
		if (exService == null) {
			exService = Executors.newWorkStealingPool(options.getNrThreads());
		}
		
		for (Feature region : regions) {
			// make output directory
			
			// filter variants
			ArrayList<VCFVariant> regionVariants = getRegionVariants(variants, region);
			ArrayList<AssociationResult> regionAssocResults = getRegionAssocResults(assocResults, region);
			if (!regionVariants.isEmpty() && !regionAssocResults.isEmpty()) {
				// get the intersect between variants and associationresults
				System.out.println("Processing region: " + region.toString());
				HashSet<String> variantIDSetAssocResults = new HashSet<String>();
				for (int i = 0; i < regionAssocResults.size(); i++) {
					AssociationResult r = regionAssocResults.get(i);
					SNPFeature snp = r.getSnp();
					String id = snp.getChromosome().toString() + "_" + snp.getStart() + "_" + snp.getName() + "_" + snp.getAlleles()[0] + "_" + snp.getAlleles()[1];
					variantIDSetAssocResults.add(id);
				}
				
				HashMap<String, Integer> variantIndex = new HashMap<>();
				int ctr = 0;
				ArrayList<VCFVariant> sharedVariants = new ArrayList<>();
				for (int i = 0; i < regionVariants.size(); i++) {
					VCFVariant v = regionVariants.get(i);
					String id = v.getChrObj().toString() + "_" + v.getPos() + "_" + v.getId() + "_" + v.getAlleles()[0] + "_" + v.getAlleles()[1];
					if (variantIDSetAssocResults.contains(id)) {
						if (variantIndex.containsKey(id)) {
							System.err.println("Error. Index already contains variant: " + id);
						} else {
							variantIndex.put(id, ctr);
							ctr++;
							sharedVariants.add(v);
						}
					}
				}
				
				if (sharedVariants.isEmpty()) {
					System.err.println("No shared variants for region: " + region.toString());
				} else {
					System.out.println(sharedVariants.size() + " variants to export..");
					// dose file - space sep
					//sampleID SNP1 SNP2 SNP3 SNP4 SNP5
					//IND1 0 0 0 0 1
					
					
					ArrayList<VCFVariant> tmp = new ArrayList<>();
					for (VCFVariant v : sharedVariants) {
						if (v.isBiallelic()) {
							tmp.add(v);
						}
					}
					sharedVariants = tmp;
					
					// write the file
					
					
					ProgressBar pb2 = new ProgressBar(samplenames.length, "Prepping:");
					float[][] output = new float[samplenames.length][sharedVariants.size()];
					for (int v = 0; v < sharedVariants.size(); v++) {
						VCFVariant var = sharedVariants.get(v);
						DoubleMatrix2D dosages = var.getDosagesAsMatrix2D();
						for (int s = 0; s < samplenames.length; s++) {
							double dosage = dosages.getQuick(s, 0);
							if (!Double.isNaN(dosage) && dosage > -1) {
								output[s][v] = (float) dosage;
							} else {
								output[s][v] = -9;
							}
						}
						pb2.iterate();
						
					}
					pb2.close();
					
					TextFile doseout = new TextFile(options.getOutputdir() + "region-" + region.toString() + ".dosage", TextFile.W);
					String headerdosage = "sampleID";
					for (int v = 0; v < sharedVariants.size(); v++) {
						VCFVariant var = sharedVariants.get(v);
						headerdosage += " " + var.asSNPFeature().toString();
					}
					doseout.writeln(headerdosage);
					ProgressBar pb = new ProgressBar(samplenames.length, "Writing: " + options.getOutputdir() + "region-" + region.toString() + ".dosage");
					for (int s = 0; s < samplenames.length; s++) {
						StringBuilder doseln = new StringBuilder(samplenames[s]);
						for (int v = 0; v < sharedVariants.size(); v++) {
							double dosage = output[s][v];
							doseln.append(" ").append(dosage);
						}
						pb.iterate();
						doseout.writeln(doseln.toString());
					}
					pb.close();
					doseout.close();
					
					TextFile mapout = new TextFile(options.getOutputdir() + "region-" + region.toString() + ".map", TextFile.W);
					for (int v = 0; v < sharedVariants.size(); v++) {
						VCFVariant var = sharedVariants.get(v);
						
						String lnout = var.asSNPFeature().getChromosome().getNumber() + " " + var.asSNPFeature().toString() + " 0 " + var.getPos();
						mapout.writeln(lnout);
					}
					mapout.close();
					
					TextFile pedout = new TextFile(options.getOutputdir() + "region-" + region.toString() + ".ped", TextFile.W);
					ProgressBar pb3 = new ProgressBar(samplenames.length, "Writing: " + options.getOutputdir() + "region-" + region.toString() + ".ped");
					for (int s = 0; s < samplenames.length; s++) {
						StringBuilder doseln = new StringBuilder(samplenames[s])
								.append(" ").append(samplenames[s])
								.append(" 0 0 0 0");
						
						for (int v = 0; v < sharedVariants.size(); v++) {
							VCFVariant var = sharedVariants.get(v);
							String[] alleles = var.getAlleles();
							double dosage = output[s][v];
							if (dosage < 0) {
								doseln.append(" ").append("0").append(" ").append("0");
							} else if (dosage < 0.5) {
								doseln.append(" ").append(alleles[0]).append(" ").append(alleles[0]);
							} else if (dosage > 1.5) {
								doseln.append(" ").append(alleles[1]).append(" ").append(alleles[1]);
							} else {
								doseln.append(" ").append(alleles[0]).append(" ").append(alleles[1]);
							}
						}
						pb3.iterate();
						pedout.writeln(doseln.toString());
					}
					pb3.close();
					pedout.close();
				}
			}
		}
		
		if (exService != null) {
			exService.shutdown();
		}
		
	}
	
	private ArrayList<AssociationResult> getRegionAssocResults(ArrayList<AssociationResult> assocResults, Feature region) {
		ArrayList<AssociationResult> output = new ArrayList<>();
		for (int r = 0; r < assocResults.size(); r++) {
			AssociationResult result = assocResults.get(r);
			if (result.getSnp().overlaps(region)) {
				output.add(result);
			}
		}
		return output;
	}
	
	public ArrayList<VCFVariant> getRegionVariants(ArrayList<VCFVariant> variants, Feature region) {
		ArrayList<VCFVariant> output = new ArrayList<VCFVariant>();
		for (VCFVariant v : variants) {
			if (!v.isMultiallelic() && v.asFeature().overlaps(region)) {
				output.add(v);
			}
		}
		return output;
	}
	
}
