package nl.harmjanwestra.vcfutils;

import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureTree;
import nl.harmjanwestra.utilities.genotypes.GenotypeTools;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.NavigableSet;

/**
 * Created by hwestra on 5/2/16.
 */
public class VCFVariantStats {

	public static void main(String[] args) {
		Main m = new Main();
		args = new String[]{
				"--summarize3compare",
				"-i", "/Data/tmp/2016-05-03/cosmo.stats.vcf.gz",
				"-i2", "/Data/tmp/2016-05-03/cosmo.stats.vcf.gz",
				"-o", "/Data/tmp/2016-05-03/CosmoVsHRC",
				"-b", "/Data/tmp/2016-05-03/AllICLoci.bed"
		};
		Main.main(args);

	}

	public void run(String vcf1, String vcf2) throws IOException {

		TextFile out = new TextFile(vcf2, TextFile.W);
		System.out.println("Out: " + vcf2);
		String[] vcfFiles = vcf1.split(",");
		for (int i = 0; i < vcfFiles.length; i++) {
			String file = vcfFiles[i];

			if (Gpio.exists(file)) {
				TextFile tf = new TextFile(file, TextFile.R);
				String ln = tf.readLine();
				System.out.println("Parsing: " + file);
				int lnnum = 0;
				while (ln != null) {
					if (ln.startsWith("##")) {
						if (i == 0) {
							out.writeln(ln);
						}
					} else if (ln.startsWith("#CHROM")) {

						if (i == 0) {
							out.writeln("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
						}

					} else {

						String chrStr = ln.substring(0, 100);
						String[] chrStrElems = chrStr.split("\t");
						Chromosome chr = Chromosome.parseChr(chrStrElems[0]);
						if (!chr.equals(Chromosome.X)) {

							VCFVariant variant = new VCFVariant(ln, VCFVariant.PARSE.ALL);

							// AC / AN / AF
							String AN = "AN=" + variant.getTotalAlleleCount();
							String AF = "AF=" + Strings.concat(variant.getAlleleFrequencies(), Strings.comma, 1, variant.getAlleles().length);
							String AC = "AC=" + Strings.concat(variant.getNrAllelesObserved(), Strings.comma, 1, variant.getAlleles().length);
							String INFO = "INFO=" + variant.getImputationQualityScore();

							String infoString = AC + ";" + AF + ";" + AN + ";" + INFO;

							StringBuilder outbuilder = new StringBuilder(1000);
							outbuilder.append(variant.getChr())
									.append("\t").append(variant.getPos())
									.append("\t").append(variant.getId())
									.append("\t").append(variant.getAlleles()[0])
									.append("\t").append(Strings.concat(variant.getAlleles(), Strings.comma, 1, variant.getAlleles().length))
									.append("\t.").append("\t.")
									.append("\t").append(infoString);

							out.writeln(outbuilder.toString());
						}

					}
					ln = tf.readLine();
					lnnum++;
					if (lnnum % 10000 == 0) {
						System.out.print(lnnum + " lines parsed.\r");
					}
				}
				tf.close();
			}
		}
		System.out.println();
		System.out.println("Done");
		out.close();


	}


	public void compare(String vcf1, String vcf2, String out, String regionsFile, String filterfile, boolean filterExcludes) throws IOException {
		HashSet<Feature> filter = null;
		if (filterfile != null) {
			filter = new HashSet<Feature>();
			TextFile tf = new TextFile(filterfile, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				if (elems.length >= 3 && !elems[0].startsWith("#")) {
					Feature f = new Feature();
					Chromosome chr = Chromosome.parseChr(elems[0]);
					int pos = Integer.parseInt(elems[1]);
					String id = elems[2];
					f.setName(id);
					f.setChromosome(chr);
					f.setStart(pos);
					f.setStop(pos);
					filter.add(f);
				}
				elems = tf.readLineElems(TextFile.tab);
			}

			tf.close();
			System.out.println(filter.size() + " variants loaded from " + filterfile);
		}

		BedFileReader bedFileReader = new BedFileReader();
		ArrayList<Feature> bedFeatures = bedFileReader.readAsList(regionsFile);

		System.out.println(bedFeatures.size() + " regions in " + regionsFile);

		ArrayList<VCFVariant> variants1 = loadVariants(vcf1, bedFeatures, filter, filterExcludes);
		System.out.println(variants1.size() + " variants in " + vcf1);
		ArrayList<VCFVariant> variants2 = loadVariants(vcf2, bedFeatures, filter, filterExcludes);
		System.out.println(variants1.size() + " variants in " + vcf2);

		GenotypeTools tools = new GenotypeTools();
		VCFMerger merger = new VCFMerger();

		TextFile outf = new TextFile(out + "shared.txt", TextFile.W);
		TextFile outfailed = new TextFile(out + "failedcomp.txt", TextFile.W);
		TextFile outmultiallelic = new TextFile(out + "multiallelic.txt", TextFile.W);
		TextFile outunique1 = new TextFile(out + "uniqueTo1.txt", TextFile.W);
		TextFile outunique2 = new TextFile(out + "uniqueTo2.txt", TextFile.W);

		outunique1.writeln("#input:" + vcf1);
		outunique1.writeln("Chr\tPos\tSNP\tRef\tAlt\tNrAlleles1\tAN\tAF\tAC\tINFO\tMAF");

		outunique2.writeln("#input:" + vcf2);
		outunique2.writeln("Chr\tPos\tSNP\tRef\tAlt\tNrAlleles\tAN\tAF\tAC\tINFO\tMAF");

		outf.writeln("Chr\tPos\tSNP1\tRef1\tAlt1\tNrAlleles1\tAN1\tAF1\tAC1\tINFO1\tMAF1\t" +
				"SNP2\tRef2\tAlt2\tNrAlleles2\tAN2\tAF2\tAC2\tINFO2\tMAF2\tEqualID\tNrEqualAlleles\tComplement\tReferenceEqual");

		int shared = 0;
		int uniquein1 = 0;
		int uniquein2 = 0;

		int multiAllelic = 0;

		for (
				int chr = 1;
				chr < 23; chr++)

		{
			// hash the variants in list 1
			HashMap<Feature, ArrayList<VCFVariant>> set = filter(variants1, chr);

			// get overlapping variants
			for (int j = 0; j < variants2.size(); j++) {
				VCFVariant var2 = variants2.get(j);
				if (var2.getChrObj().getNumber() == chr) {

					if (var2.getAlleles().length > 2) {
						multiAllelic++;
					}

					Feature f = var2.asFeature();
					f.setName(null);
					ArrayList<VCFVariant> variantsIn1 = set.get(f);
					if (variantsIn1 == null) {

						if (var2.getAlleles().length == 2) {

							String af = var2.getInfo().get("AF");
							Double maf = null;
							if (af != null) {
								maf = Double.parseDouble(af);
								if (maf > 0.5) {
									maf = 1 - maf;
								}

							}

							String outstr = var2.getChr() + "\t" + var2.getPos()
									+ "\t" + var2.getId() + "\t" + var2.getAlleles()[0] + "\t" + Strings.concat(var2.getAlleles(), Strings.comma, 1, var2.getAlleles().length) + "\t" + var2.getAlleles().length
									+ "\t" + var2.getInfo().get("AN") + "\t" + var2.getInfo().get("AF") + "\t" + var2.getInfo().get("AC") + "\t" + var2.getInfo().get("INFO") + "\t" + maf;
							outunique2.writeln(outstr);
							uniquein2++;
						}

					} else {
						shared++;
						for (int k = 0; k < variantsIn1.size(); k++) {
							VCFVariant var1 = variantsIn1.get(k);
							// identical
							String[] alleles1 = var1.getAlleles();
							String[] alleles2 = var2.getAlleles();
							if (alleles1.length == 2 && alleles2.length == 2) {

								int nrIdenticalAlleles = merger.countIdenticalAlleles(alleles1, alleles2);
								boolean complement = false;
								if (nrIdenticalAlleles < 2) {
									// try complement
									var1.convertAllelesToComplement();
									alleles1 = var1.getAlleles();
									nrIdenticalAlleles = merger.countIdenticalAlleles(alleles1, alleles2);
									complement = true;
								}

								if (nrIdenticalAlleles == 2) {
									// consider as a unique variant in stead
									boolean referenceEqual = true;
									// check whether the reference alleles are identical
									if (!var1.getAlleles()[0].equals(var2.getAlleles()[0])) {
										// check which allele corresponds to the reference allele
										referenceEqual = false;

										Double AN = Double.parseDouble(var2.getInfo().get("AN"));
										Double AF = Double.parseDouble(var2.getInfo().get("AF"));
										Double AC = Double.parseDouble(var2.getInfo().get("AC"));

										var2.getInfo().put("AF", "" + (1d - AF));
										var2.getInfo().put("AC", "" + (AN - AC));
										var2.flipReferenceAlelele();
									}

									String af2 = var2.getInfo().get("AF");
									Double maf2 = null;
									if (af2 != null) {
										maf2 = Double.parseDouble(af2);
										if (maf2 > 0.5) {
											maf2 = 1 - maf2;
										}
									}

									String af1 = var1.getInfo().get("AF");
									Double maf1 = null;
									if (af1 != null) {
										maf1 = Double.parseDouble(af1);
										if (maf1 > 0.5) {
											maf1 = 1 - maf1;
										}
									}


									String outStr = var1.getChr() + "\t" + var1.getPos()
											+ "\t" + var1.getId() + "\t" + var1.getAlleles()[0] + "\t" + Strings.concat(var1.getAlleles(), Strings.comma, 1, var1.getAlleles().length) + "\t" + var1.getAlleles().length
											+ "\t" + var1.getInfo().get("AN") + "\t" + var1.getInfo().get("AF") + "\t" + var1.getInfo().get("AC") + "\t" + var1.getInfo().get("INFO") + "\t" + maf1
											+ "\t" + var2.getId() + "\t" + var2.getAlleles()[0] + "\t" + Strings.concat(var2.getAlleles(), Strings.comma, 1, var2.getAlleles().length) + "\t" + var2.getAlleles().length
											+ "\t" + var2.getInfo().get("AN") + "\t" + var2.getInfo().get("AF") + "\t" + var2.getInfo().get("AC") + "\t" + var2.getInfo().get("INFO") + "\t" + maf2
											+ "\t" + var1.getId().equals(var2.getId())
											+ "\t" + nrIdenticalAlleles
											+ "\t" + complement
											+ "\t" + referenceEqual;
									outf.writeln(outStr);
								} else {
									String af2 = var2.getInfo().get("AF");
									Double maf2 = null;
									if (af2 != null) {
										maf2 = Double.parseDouble(af2);
										if (maf2 > 0.5) {
											maf2 = 1 - maf2;
										}
									}

									String af1 = var1.getInfo().get("AF");
									Double maf1 = null;
									if (af1 != null) {
										maf1 = Double.parseDouble(af1);
										if (maf1 > 0.5) {
											maf1 = 1 - maf1;
										}
									}
									String outStr = var1.getChr() + "\t" + var1.getPos()
											+ "\t" + var1.getId() + "\t" + var1.getAlleles()[0] + "\t" + Strings.concat(var1.getAlleles(), Strings.comma, 1, var1.getAlleles().length) + "\t" + var1.getAlleles().length
											+ "\t" + var1.getInfo().get("AN") + "\t" + var1.getInfo().get("AF") + "\t" + var1.getInfo().get("AC") + "\t" + var1.getInfo().get("INFO") + "\t" + maf1
											+ "\t" + var2.getId() + "\t" + var2.getAlleles()[0] + "\t" + Strings.concat(var2.getAlleles(), Strings.comma, 1, var2.getAlleles().length) + "\t" + var2.getAlleles().length
											+ "\t" + var2.getInfo().get("AN") + "\t" + var2.getInfo().get("AF") + "\t" + var2.getInfo().get("AC") + "\t" + var2.getInfo().get("INFO") + "\t" + maf2
											+ "\t" + var1.getId().equals(var2.getId())
											+ "\t" + nrIdenticalAlleles
											+ "\t" + complement
											+ "\t" + null;
									outfailed.writeln(outStr);
								}


								if (complement) {
									// change back to original orientation
									var1.convertAllelesToComplement();
								}
							} else {
								outmultiallelic.writeln();
							}
						}
					}
				}

			}

			set = filter(variants2, chr);
			for (int j = 0; j < variants1.size(); j++) {
				VCFVariant var1 = variants1.get(j);
				if (var1.getChrObj().getNumber() == chr) {
					Feature f = var1.asFeature();
					f.setName(null);
					if (var1.getAlleles().length == 2) {
						ArrayList<VCFVariant> variantsIn2 = set.get(f);
						if (variantsIn2 == null) {
							uniquein1++;
							String af1 = var1.getInfo().get("AF");
							Double maf1 = null;
							if (af1 != null) {
								maf1 = Double.parseDouble(af1);
								if (maf1 > 0.5) {
									maf1 = 1 - maf1;
								}
							}
							String outstr = var1.getChr() + "\t" + var1.getPos()
									+ "\t" + var1.getId() + "\t" + var1.getAlleles()[0] + "\t" + Strings.concat(var1.getAlleles(), Strings.comma, 1, var1.getAlleles().length) + "\t" + var1.getAlleles().length
									+ "\t" + var1.getInfo().get("AN") + "\t" + var1.getInfo().get("AF") + "\t" + var1.getInfo().get("AC") + "\t" + var1.getInfo().get("INFO") + "\t" + maf1;
							outunique1.writeln(outstr);
						}
					}
				}
			}
		}

		outf.close();
		outfailed.close();
		outmultiallelic.close();
		outunique1.close();
		outunique2.close();

		System.out.println(shared + " variants shared.");
		System.out.println(uniquein1 + " variants unique in " + vcf1);
		System.out.println(uniquein2 + " variants unique in " + vcf2);
		System.out.println(multiAllelic + " multi allelic");

	}

	private HashMap<Feature, ArrayList<VCFVariant>> filter(ArrayList<VCFVariant> variants1, int i) {
		HashMap<Feature, ArrayList<VCFVariant>> output = new HashMap<>();
		for (int q = 0; q < variants1.size(); q++) {
			VCFVariant var = variants1.get(q);
			if (var.getChrObj().getNumber() == i) {
				Feature f = var.asFeature();
				f.setName(null); // just in case.
				ArrayList<VCFVariant> list = output.get(f);
				if (list == null) {
					list = new ArrayList<VCFVariant>();
				}
				list.add(var);
				output.put(f, list);
			}
		}
		return output;
	}

	public ArrayList<VCFVariant> loadVariants(String vcf, ArrayList<Feature> bedFeatures, HashSet<Feature> filter, boolean filterExcludes) throws IOException {
		ArrayList<VCFVariant> variants1 = new ArrayList<VCFVariant>();


		FeatureTree tree = new FeatureTree(bedFeatures, true);

		System.out.println(tree.getTreeSize() + " regions in tree, " + tree.getListSize() + " regions in list");


		int loaded = 0;
		int parsed = 0;
		int filtered = 0;
		VCFGenotypeData data1 = new VCFGenotypeData(vcf);
		VCFVariant variant1 = data1.nextLoadHeader();
		while (data1.hasNext()) {
			Feature f1 = variant1.asFeature();
			NavigableSet<Feature> overlappingSet = tree.getFeatureSet(f1);
			if (!overlappingSet.isEmpty()) {
				if (filter == null || (filterExcludes && !filter.contains(f1)) || (!filterExcludes && filter.contains(f1))) {
					variants1.add(variant1);
					loaded++;
				} else {
					filtered++;
				}

			} else {
				overlappingSet = tree.getFeatureSet(f1.getChromosome(), 0, f1.getChromosome().getLength());

				for (Feature f : overlappingSet) {
					if (f.overlaps(f1)) {
						System.out.println("excluded: " + f1.toBedString() + " but overlaps: " + f.toBedString());
					}
				}
			}

			variant1 = data1.nextLoadHeader();
			parsed++;
			if (parsed % 10000 == 0) {
				System.out.print(parsed + " variants parsed. " + loaded + " in memory.\r");
			}
		}
		System.out.println();
		System.out.println(parsed + " variants parsed. " + loaded + " in memory. " + filtered + " filtered");
		System.out.println("Done loading: " + vcf);
		return variants1;
	}

}
