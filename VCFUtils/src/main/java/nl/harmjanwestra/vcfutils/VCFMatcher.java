package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.genotypes.GenotypeTools;
import nl.harmjanwestra.utilities.vcf.VCFFunctions;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 3/10/16.
 */
public class VCFMatcher {


	public void run(String refVCF,
	                String testVCF,
	                String outprefix,
	                boolean linux,
	                String vcfsort) throws IOException {


		String logoutfile = outprefix + "-log.txt";

		String vcfout = outprefix + "-matched.vcf.gz";
		String uniquevarsout = outprefix + "-uniqueVars.txt";
		System.out.println("Merging: ");
		System.out.println("ref: " + refVCF);
		System.out.println("test: " + testVCF);
		System.out.println("out: " + vcfout);


		// pos --> Feature / Alleles / minor allele
		HashMap<Feature, ArrayList<Triple<String[], String, Double>>> refData = new HashMap<>(1000000);

		// reload the variants, get the positions with multiple variants at that position
		VCFGenotypeData iterator = new VCFGenotypeData(refVCF);


		int nrLoaded = 0;
		System.out.println("Inventorizing variants in: " + refVCF);
		while (iterator.hasNext()) {
			VCFVariant f = iterator.next();
			Triple<String[], String, Double> data = new Triple<String[], String, Double>(f.getAlleles(), f.getMinorAllele(), f.getMAF());
			ArrayList<Triple<String[], String, Double>> siteData = refData.get(f.asFeature());
			if (siteData == null) {
				siteData = new ArrayList<>();
			}
			siteData.add(data);
			refData.put(f.asFeature(), siteData);
			nrLoaded++;
			if (nrLoaded % 10000 == 0) {
				System.out.println(nrLoaded + " loaded so far");
			}
		}

		System.out.println(refData.size() + " features after loading ref vcf");
		iterator.close();


		TextFile logOut = new TextFile(logoutfile, TextFile.W);
		logOut.writeln("chr\t" +
				"pos\t" +
				"refVariant#atPos\t" +
				"refVariantId\t" +
				"refVariantRefAllele\t" +
				"refVariantAltAlleles\t" +
				"refVariantMinorAllele\t" +
				"refVariantMAF\t" +
				"testVariantId\t" +
				"testVariantRefAllele\t" +
				"testVariantAltAlleles\t" +
				"testVariantMinorAllele\t" +
				"testVariantMAF\t" +
				"testVariantRefAlleleAM\t" +
				"testVariantAltAllelesAM\t" +
				"testVariantMinorAlleleAM\t" +
				"Reason");


		TextFile uniqueTest = new TextFile(uniquevarsout, TextFile.W);
		VCFMerger merger = new VCFMerger();

		TextFile outvcf = new TextFile(vcfout, TextFile.W);
		TextFile headerin = new TextFile(testVCF, TextFile.R);
		String ln = headerin.readLine();
		while (ln != null) {
			if (ln.startsWith("#")) {
				outvcf.writeln(ln);
			}
			ln = headerin.readLine();
		}
		headerin.close();

		VCFGenotypeData iterator2 = new VCFGenotypeData(testVCF);

		int nrWritten = 0;
		int nrNotWritten = 0;
		while (iterator2.hasNext()) {
			boolean variantAlreadyWritten = false;
			VCFVariant testVar = iterator2.next();
			ArrayList<Triple<String[], String, Double>> referenceData = refData.get(testVar.asFeature());
			if (referenceData == null) {
				// write to unique file
				uniqueTest.writeln(testVar.getChr() + "\t" + testVar.getPos() + "\t" + testVar.getId());
			} else {

				for (Triple<String[], String, Double> refVariant : referenceData) {

					String logln = testVar.getChr()
							+ "\t" + testVar.getPos()
							+ "\t" + referenceData.size()
							+ "\t" + testVar.getId()
							+ "\t" + refVariant.getLeft()[0]
							+ "\t" + Strings.concat(refVariant.getLeft(), Strings.comma, 1, refVariant.getLeft().length)
							+ "\t" + refVariant.getMiddle()
							+ "\t" + refVariant.getRight();

					String logoutputln = logln;


					logoutputln += "\t" + testVar.getId()
							+ "\t" + testVar.getAlleles()[0]
							+ "\t" + Strings.concat(testVar.getAlleles(), Strings.comma, 1, testVar.getAlleles().length)
							+ "\t" + testVar.getMinorAllele()
							+ "\t" + testVar.getMAF();


					String[] refAlleles = refVariant.getLeft();
					String refMinorAllele = refVariant.getMiddle();
					String[] testVariantAlleles = testVar.getAlleles();
					String testVariantMinorAllele = testVar.getMinorAllele();


					int nridenticalalleles = merger.countIdenticalAlleles(refAlleles, testVariantAlleles);
					GenotypeTools gtools = new GenotypeTools();

					boolean complement = false;
					if (nridenticalalleles == 0) {
						// try complement
						complement = true;
						String[] complementAlleles2 = gtools.convertToComplement(testVariantAlleles);
						testVariantMinorAllele = gtools.getComplement(testVariantMinorAllele);
						if (testVariantMinorAllele == null) {
							nridenticalalleles = 0;
						} else {
							nridenticalalleles = merger.countIdenticalAlleles(refAlleles, complementAlleles2);
						}
					}

					boolean flipped = false;
					boolean writeVariant = false;
					if (refVariant.getLeft().length == 2 && testVar.getAlleles().length == 2) {
						if (nridenticalalleles == 2) {
							if (testVariantMinorAllele.equals(refMinorAllele) || (testVar.getMAF() > 0.45 && refVariant.getRight() > 0.45)) {
								// check whether the reference allele is equal
								String[] tmpAlleles = testVariantAlleles;
								if (complement) {
									testVar.convertAllelesToComplement();
									tmpAlleles = testVar.getAlleles();
								}

								if (!refAlleles[0].equals(tmpAlleles[0])) {
									testVar.flipReferenceAlelele();
									flipped = true;
								}

								logoutputln += "\t" + testVar.getAlleles()[0] + "\t" + Strings.concat(testVar.getAlleles(), Strings.comma, 1, testVar.getAlleles().length);

								writeVariant = true;
								if (complement) {
									logoutputln += "\tOK-Complement";
								} else {
									logoutputln += "\tOK";
								}
								if (flipped) {
									logoutputln += "-flippedAlleles";
								}


							} else {
								// write to log?
								logoutputln += "\t-\t-\tNotOK-DiffMinor";

							}
						} else {
							// write to log?
							logoutputln += "\t-\t-\tNotOK-IncompatibleAlleles";

						}
					} else if (nridenticalalleles > 1) {
						// recode the genotypes towards the joint set of alleles
						// get a list of all alleles at this locus...
						HashSet<String> uniqueAlleles = new HashSet<String>();
						uniqueAlleles.addAll(Arrays.asList(refAlleles));
						uniqueAlleles.addAll(Arrays.asList(testVariantAlleles));

						HashMap<String, Integer> alleleMap = new HashMap<String, Integer>();
						ArrayList<String> newAlleles = new ArrayList<String>();
						for (int i = 0; i < refAlleles.length; i++) {
							alleleMap.put(refAlleles[i], i);
							newAlleles.add(refAlleles[i]);
						}

						for (int i = 0; i < testVariantAlleles.length; i++) {
							String alleleStr = testVariantAlleles[i];
							if (!alleleMap.containsKey(alleleStr)) {
								alleleMap.put(alleleStr, alleleMap.size());
								newAlleles.add(alleleStr);
							}
						}

						// recode testVariant
						testVar.recodeAlleles(alleleMap, newAlleles.toArray(new String[0]));
						logoutputln += "\t" + testVar.getAlleles()[0] + "\t" + Strings.concat(testVar.getAlleles(), Strings.comma, 1, testVar.getAlleles().length);

						// mergecheese
						logoutputln += "\tOK-MultiAllelic-AllelesRecoded";
					}

					if (writeVariant && !variantAlreadyWritten) {
						outvcf.writeln(testVar.toVCFString());
						variantAlreadyWritten = true;
					}
					logOut.writeln(logoutputln);
				}
			}
			if (variantAlreadyWritten) {
				nrWritten++;
			} else {
				nrNotWritten++;
			}
		}

		System.out.println(nrWritten + " written (" + ((double) nrWritten / (nrWritten + nrNotWritten)) + ") \t" + nrNotWritten + " not written\t" + (nrWritten + nrNotWritten) + " total");
		uniqueTest.close();
		outvcf.close();
		logOut.close();

		VCFFunctions t = new VCFFunctions();

		t.sortVCF(linux, vcfsort, vcfout, outprefix + "-matched-sorted.vcf.gz", outprefix + "-matched-sort.sh");
		Gpio.delete(outprefix + "-matched-sort.sh");
		Gpio.delete(vcfout);
	}


	public void runReferenceOnlyHasVariantPositions(String refVCF,
	                                                String testVCF,
	                                                String outprefix,
	                                                boolean linux,
	                                                String vcfsort) throws IOException {


		String logoutfile = outprefix + "-log.txt";

		String vcfout = outprefix + "-matched.vcf.gz";
		String uniquevarsout = outprefix + "-uniqueVars.txt";
		System.out.println("Merging: ");
		System.out.println("ref: " + refVCF);
		System.out.println("test: " + testVCF);
		System.out.println("out: " + vcfout);


		// pos --> Feature / Alleles / minor allele
		HashMap<Feature, ArrayList<VCFVariant>> refData = new HashMap<>(1000000);

		// reload the variants, get the positions with multiple variants at that position
		VCFGenotypeData iterator = new VCFGenotypeData(refVCF);


		int nrLoaded = 0;
		System.out.println("Inventorizing variants in: " + refVCF);
		while (iterator.hasNext()) {
			VCFVariant f = iterator.nextLoadHeader();

			String minorAllele = f.getMinorAlleleFromInfoField();
			if (minorAllele != null) {
				ArrayList<VCFVariant> siteData = refData.get(f.asFeature());
				if (siteData == null) {
					siteData = new ArrayList<>();
				}
				siteData.add(f);
				refData.put(f.asFeature(), siteData);
				nrLoaded++;
				if (nrLoaded % 10000 == 0) {
					System.out.println(nrLoaded + " loaded so far");
				}
			}
		}

		System.out.println(refData.size() + " features after loading ref vcf");
		iterator.close();


		TextFile logOut = new TextFile(logoutfile, TextFile.W);
		logOut.writeln("chr\t" +
				"pos\t" +
				"refVariant#atPos\t" +
				"refVariantId\t" +
				"refVariantRefAllele\t" +
				"refVariantAltAlleles\t" +
				"refVariantMinorAllele\t" +
				"refVariantMAF\t" +
				"testVariantId\t" +
				"testVariantRefAllele\t" +
				"testVariantAltAlleles\t" +
				"testVariantMinorAllele\t" +
				"testVariantMAF\t" +
				"testVariantRefAlleleAM\t" +
				"testVariantAltAllelesAM\t" +
				"testVariantMinorAlleleAM\t" +
				"Reason");


		TextFile uniqueTest = new TextFile(uniquevarsout, TextFile.W);
		VCFMerger merger = new VCFMerger();

		TextFile outvcf = new TextFile(vcfout, TextFile.W);
		TextFile headerin = new TextFile(testVCF, TextFile.R);
		String ln = headerin.readLine();
		while (ln != null) {
			if (ln.startsWith("#")) {
				outvcf.writeln(ln);
			}
			ln = headerin.readLine();
		}
		headerin.close();

		VCFGenotypeData iterator2 = new VCFGenotypeData(testVCF);

		int nrWritten = 0;
		int nrNotWritten = 0;
		while (iterator2.hasNext()) {
			boolean variantAlreadyWritten = false;
			VCFVariant testVar = iterator2.next();
			ArrayList<VCFVariant> referenceData = refData.get(testVar.asFeature());
			if (referenceData == null) {
				// write to unique file
				uniqueTest.writeln(testVar.getChr() + "\t" + testVar.getPos() + "\t" + testVar.getId());
			} else {
				for (VCFVariant refVariant : referenceData) {

					String[] refAlleles = refVariant.getAlleles();
					double refVariantMinorAlleleFreq = 0;
					String refMinorAllele = refVariant.getMinorAllele();


					String logln = testVar.getChr()
							+ "\t" + testVar.getPos()
							+ "\t" + referenceData.size()
							+ "\t" + testVar.getId()
							+ "\t" + refAlleles[0]
							+ "\t" + Strings.concat(refAlleles, Strings.comma, 1, refAlleles.length)
							+ "\t" + refMinorAllele
							+ "\t" + refVariantMinorAlleleFreq;

					String logoutputln = logln;


					logoutputln += "\t" + testVar.getId()
							+ "\t" + testVar.getAlleles()[0]
							+ "\t" + Strings.concat(testVar.getAlleles(), Strings.comma, 1, testVar.getAlleles().length)
							+ "\t" + testVar.getMinorAllele()
							+ "\t" + testVar.getMAF();


					String[] testVariantAlleles = testVar.getAlleles();
					String testVariantMinorAllele = testVar.getMinorAllele();


					int nridenticalalleles = merger.countIdenticalAlleles(refAlleles, testVariantAlleles);
					GenotypeTools gtools = new GenotypeTools();

					boolean complement = false;
					if (nridenticalalleles == 0) {
						// try complement
						complement = true;
						String[] complementAlleles2 = gtools.convertToComplement(testVariantAlleles);
						testVariantMinorAllele = gtools.getComplement(testVariantMinorAllele);
						if (testVariantMinorAllele == null) {
							nridenticalalleles = 0;
						} else {
							nridenticalalleles = merger.countIdenticalAlleles(refAlleles, complementAlleles2);
						}
					}

					boolean flipped = false;
					boolean writeVariant = false;
					if (refAlleles.length == 2 && testVar.getAlleles().length == 2) {
						if (nridenticalalleles == 2) {
							if (testVariantMinorAllele.equals(refMinorAllele) || (testVar.getMAF() > 0.45 && refVariant.getMAF() > 0.45)) {
								// check whether the reference allele is equal
								String[] tmpAlleles = testVariantAlleles;
								if (complement) {
									testVar.convertAllelesToComplement();
									tmpAlleles = testVar.getAlleles();
								}

								if (!refAlleles[0].equals(tmpAlleles[0])) {
									testVar.flipReferenceAlelele();
									flipped = true;
								}

								logoutputln += "\t" + testVar.getAlleles()[0] + "\t" + Strings.concat(testVar.getAlleles(), Strings.comma, 1, testVar.getAlleles().length);

								writeVariant = true;
								if (complement) {
									logoutputln += "\tOK-Complement";
								} else {
									logoutputln += "\tOK";
								}
								if (flipped) {
									logoutputln += "-flippedAlleles";
								}


							} else {
								// write to log?
								logoutputln += "\t-\t-\tNotOK-DiffMinor";

							}
						} else {
							// write to log?
							logoutputln += "\t-\t-\tNotOK-IncompatibleAlleles";

						}
					} else if (nridenticalalleles > 1) {
						// recode the genotypes towards the joint set of alleles
						// get a list of all alleles at this locus...
						HashSet<String> uniqueAlleles = new HashSet<String>();
						uniqueAlleles.addAll(Arrays.asList(refAlleles));
						uniqueAlleles.addAll(Arrays.asList(testVariantAlleles));

						HashMap<String, Integer> alleleMap = new HashMap<String, Integer>();
						ArrayList<String> newAlleles = new ArrayList<String>();
						for (int i = 0; i < refAlleles.length; i++) {
							alleleMap.put(refAlleles[i], i);
							newAlleles.add(refAlleles[i]);
						}

						for (int i = 0; i < testVariantAlleles.length; i++) {
							String alleleStr = testVariantAlleles[i];
							if (!alleleMap.containsKey(alleleStr)) {
								alleleMap.put(alleleStr, alleleMap.size());
								newAlleles.add(alleleStr);
							}
						}

						// recode testVariant
						testVar.recodeAlleles(alleleMap, newAlleles.toArray(new String[0]));
						logoutputln += "\t" + testVar.getAlleles()[0] + "\t" + Strings.concat(testVar.getAlleles(), Strings.comma, 1, testVar.getAlleles().length);

						// mergecheese
						logoutputln += "\tOK-MultiAllelic-AllelesRecoded";
					}

					if (writeVariant && !variantAlreadyWritten) {
						outvcf.writeln(testVar.toVCFString());
						variantAlreadyWritten = true;
					}
					logOut.writeln(logoutputln);
				}
			}
			if (variantAlreadyWritten) {
				nrWritten++;
			} else {
				nrNotWritten++;
			}
		}

		System.out.println(nrWritten + " written (" + ((double) nrWritten / (nrWritten + nrNotWritten)) + ") \t" + nrNotWritten + " not written\t" + (nrWritten + nrNotWritten) + " total");
		uniqueTest.close();
		outvcf.close();
		logOut.close();

		VCFFunctions t = new VCFFunctions();

		t.sortVCF(linux, vcfsort, vcfout, outprefix + "-matched-sorted.vcf.gz", outprefix + "-matched-sort.sh");
		Gpio.delete(outprefix + "-matched-sort.sh");
		Gpio.delete(vcfout);
	}


}
