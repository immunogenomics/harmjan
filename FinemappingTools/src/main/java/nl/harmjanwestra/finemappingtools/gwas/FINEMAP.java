package nl.harmjanwestra.finemappingtools.gwas;

import nl.harmjanwestra.finemappingtools.gwas.CLI.LRTestOptions;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.vcf.VCFVariantComparator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 7/12/17.
 */
public class FINEMAP extends LRTest {

	public FINEMAP(LRTestOptions options) throws IOException {
		super(options);
	}

	public void run() throws IOException {


		// region1.k --> prior prob thresholds for k causal variants
		// e.g.: 0.6 0.3 0.1


		System.out.println("Running GUESS converter...");
		options.setSplitMultiAllelic(true);

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(options.getBedfile());
		ArrayList<VCFVariant> variants = readVariants(options.getOutputdir() + "variantlog.txt", regions);
		Collections.sort(variants, new VCFVariantComparator());

		if (variants.isEmpty()) {
			System.out.println("Sorry no variants found.");
			System.exit(-1);
		}

		// read the association file
		AssociationFile f = new AssociationFile();
		String assocfile = options.getAssocFile();
		ArrayList<AssociationResult> assocResults = f.read(assocfile);

		// open master file
		// master file:
		// region1.snp and region1.config, and region1.log are output files...
		// z; 	ld; 	snp; 	config; 	k; 	log; 	n-ind
		// region1.z; 	region1.ld; 	region1.snp; 	region1.config; 	region1.k; 	region1.log; 	5363

		boolean append = false;
		TextFile masterOut = null;
		if (Gpio.exists(options.getOutputdir() + "master")) {
			masterOut = new TextFile(new java.io.File(options.getOutputdir() + "master"), TextFile.MODE.APPEND);
		} else {
			masterOut = new TextFile(options.getOutputdir() + "master", TextFile.W);
			masterOut.writeln("z;\tld;\tsnp;\tconfig;\t	k;\tlog;\tn-ind");
		}

		for (Feature region : regions) {
			// make output directory

			// filter variants
			ArrayList<VCFVariant> regionVariants = getRegionVariants(variants, region);
			ArrayList<AssociationResult> associationResults = getRegionAssocResults(assocResults, region);
			if (!regionVariants.isEmpty() && !associationResults.isEmpty()) {
				// get the intersect between variants and associationresults

				HashSet<String> variantIDSetAssocResults = new HashSet<String>();
				for (int i = 0; i < associationResults.size(); i++) {
					AssociationResult r = associationResults.get(i);
					SNPFeature snp = r.getSnp();
					String id = snp.getChromosome().toString() + "_" + snp.getStart() + "_" + snp.getName() + "_" + snp.getAlleles()[0] + "_" + snp.getAlleles()[1];
					variantIDSetAssocResults.add(id);
				}

				HashMap<String, Integer> variantIndex = new HashMap<String, Integer>();
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
					// index association results, reorder
					AssociationResult[] sharedResults = new AssociationResult[sharedVariants.size()];
					for (int i = 0; i < associationResults.size(); i++) {
						AssociationResult r = associationResults.get(i);
						SNPFeature snp = r.getSnp();
						String id = snp.getChromosome().toString() + "_" + snp.getStart() + "_" + snp.getName() + "_" + snp.getAlleles()[0] + "_" + snp.getAlleles()[1];
						Integer index = variantIndex.get(id);
						if (index != null) {
							sharedResults[index] = r;
						}
					}

					// check whether there are any nulls (you never know)
					boolean run = true;
					for (int i = 0; i < sharedResults.length; i++) {
						if (sharedResults[i] == null) {
							System.err.println("One of the shared association results hasn't been set?");
							run = false;
						}
					}


					if (run) {
						// write Z-score file
						// region.z --> association Z-scores
		/*
		rs1 	3.58
rs2 	4.36
rs3 	3.71
		 */
						TextFile zOut = new TextFile(region.toString() + ".z", TextFile.W);

						zOut.close();

						// calculate correlation matrix
						// region.ld --> same order as region1.z; pearson correlations
		/*
		1.00 	0.95 	0.98
0.95 	1.00 	0.96
0.97 	0.96 	1.00
		 */
					}
				}
			}
		}
		masterOut.close();


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
