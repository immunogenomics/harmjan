package nl.harmjanwestra.finemapping.annotation;

import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFTabix;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Triple;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 10/30/16.
 */
public class BedOverlap {


	public static void main(String[] args) {

		BedOverlap overlap = new BedOverlap();
		String variantfile = "/Data/Projects/Chikashi/TAKnonHLA7SNPs.txt";
		String tabixprefix = "/Data/Ref/beagle_1kg/1kg.phase3.v5a.chr";
		String tabixSampleFile = "/Data/Ref/1kg-europeanpopulations.txt.gz";
		String bedfiles = "/Data/Enhancers/TrynkaEpigenetics/hg19/h3k4me3files.txt";
		bedfiles = "/Data/Enhancers/TrynkaEpigenetics/hg19/alldnase.txt";
		int bparoundpeakcenter = 150;
		double ld = 0.8;
		String output = "/Data/Projects/Chikashi/ChikashiH3K4me3.txt";
		output = "/Data/Projects/Chikashi/ChikashiDNASE.txt";

		try {
			overlap.run(variantfile, tabixprefix, tabixSampleFile, bedfiles, bparoundpeakcenter, ld, output);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String variantfile,
					String tabixprefix,
					String limitsamples,
					String bedfiles,
					int bpAroundCenter,
					double ldthreshold,
					String output) throws IOException {

		// load variants and define regions
		ArrayList<SNPFeature> variants = new ArrayList<SNPFeature>();
		TextFile tf = new TextFile(variantfile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length > 1) {
				String regionStr = elems[0];
				String var = elems[1];
				String[] varelems = var.split("_");
				String name = varelems[1];
				String pos = varelems[0];
				String[] poselems = pos.split(":");
				Chromosome chr = Chromosome.parseChr(poselems[0]);
				Integer posInt = Integer.parseInt(poselems[1]);
				SNPFeature feat = new SNPFeature();
				feat.setChromosome(chr);
				feat.setStart(posInt);
				feat.setStop(posInt);
				feat.setName(name);

				variants.add(feat);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		System.out.println(variants.size() + " variants found");

		// define regions based on LD
		ArrayList<Triple<SNPFeature, ArrayList<VCFVariant>, ArrayList<Double>>> variantsAndRegions = getVariants(tabixprefix, limitsamples, variants, ldthreshold);

		ArrayList<Pair<String, ArrayList<Feature>>> bedFilesAndNames = new ArrayList<>();
		TextFile in2 = new TextFile(bedfiles, TextFile.R);
		String[] elems2 = in2.readLineElems(TextFile.tab);
		while (elems2 != null) {
			String name = null;
			String file = null;
			if (elems2.length > 1) {
				name = elems2[0];
				file = elems2[1];
			} else if (elems2.length == 1) {
				name = elems2[0];
				file = elems2[0];
			}

			if (file != null) {
				BedFileReader reader = new BedFileReader();
				ArrayList<Feature> feats = reader.readAsList(file);

				if (bpAroundCenter > 0) {
					for (Feature f : feats) {
						int start = f.getStart();
						int stop = f.getStop();
						int mid = (stop - start) / 2;
						int midpos = start + mid;
						f.setStart(midpos - bpAroundCenter);
						f.setStop(midpos + bpAroundCenter);
					}
				}
				bedFilesAndNames.add(new Pair<String, ArrayList<Feature>>(name, feats));
			}


			elems2 = in2.readLineElems(TextFile.tab);
		}
		in2.close();


		// iterate variants
		TextFile outf = new TextFile(output, TextFile.W);
		outf.writeln("Variant\tProxy\tBedName\tBedFeature");
		for (Triple<SNPFeature, ArrayList<VCFVariant>, ArrayList<Double>> pair : variantsAndRegions) {

			SNPFeature origVariant = pair.getLeft();
			ArrayList<VCFVariant> proxies = pair.getMiddle();
			ArrayList<Double> ld = pair.getRight();
			if (proxies == null) {
				for (Pair<String, ArrayList<Feature>> p : bedFilesAndNames) {
					ArrayList<Feature> bedregions = p.getRight();
					for (Feature f : bedregions) {
						if (origVariant.overlaps(f)) {
							outf.writeln(origVariant.toString()
									+ "\t" + origVariant.toString()
									+ "\t" + 1
									+ "\t" + p.getLeft()
									+ "\t" + f.toString());
						}
					}
				}
			} else {
				// get region bounds
				int max = 0;
				int min = Integer.MAX_VALUE;
				for (VCFVariant v : proxies) {
					if (v.getPos() > max) {
						max = v.getPos();
					}
					if (v.getPos() < min) {
						min = v.getPos();
					}
				}

				Feature region = new Feature(pair.getLeft().getChromosome(), min, max);
				for (Pair<String, ArrayList<Feature>> p : bedFilesAndNames) {

					ArrayList<Feature> bedregions = p.getRight();
					ArrayList<Feature> filtered = filterBed(bedregions, region);

					int ctr = 0;
					for (VCFVariant variant : proxies) {
						for (Feature f : filtered) {
							if (variant.asFeature().overlaps(f)) {

								outf.writeln(origVariant.toString()
										+ "\t" + variant.asSNPFeature().toString()
										+ "\t" + ld.get(ctr)
										+ "\t" + p.getLeft()
										+ "\t" + f.toString());


							}
						}
						ctr++;
					}

				}
			}


		}
		outf.close();
	}

	private ArrayList<Feature> filterBed(ArrayList<Feature> bedregions, Feature region) {

		ArrayList<Feature> output = new ArrayList<>();
		for (Feature f : bedregions) {
			if (f.overlaps(region)) {
				output.add(f);
			}
		}
		return output;
	}

	private ArrayList<Triple<SNPFeature, ArrayList<VCFVariant>, ArrayList<Double>>> getVariants(String tabixprefix, String sampleFilter, ArrayList<SNPFeature> variants, double ldthreshold) throws IOException {

		ArrayList<Triple<SNPFeature, ArrayList<VCFVariant>, ArrayList<Double>>> output = new ArrayList<>();
		for (int i = 0; i < variants.size(); i++) {
			ArrayList<VCFVariant> proxies = new ArrayList<>();
			ArrayList<Double> ldvals = new ArrayList<>();
			SNPFeature variant = variants.get(i);

			int regionstart = variant.getStart() - 1000000;
			int regionstop = variant.getStop() + 1000000;
			if (regionstart < 0) {
				regionstart = 1;
			}
			if (regionstop > variant.getChromosome().getLengthB37()) {
				regionstop = variant.getChromosome().getLengthB37();
			}

			Feature region = new Feature(variant.getChromosome(), regionstart, regionstop);
			VCFTabix tabix = new VCFTabix(tabixprefix + variant.getChromosome().getNumber() + ".vcf.gz");
			boolean[] filter = tabix.getSampleFilter(sampleFilter);

			TabixReader.Iterator iterator = tabix.query(region);

			String ln = iterator.next();
			ArrayList<VCFVariant> tmpvariants = new ArrayList<>();
			while (ln != null) {
				tmpvariants.add(new VCFVariant(ln, VCFVariant.PARSE.ALL, filter));
				ln = iterator.next();
			}
			tabix.close();
			System.out.println(tmpvariants.size() + " variants in region " + region.toString());
			VCFVariant queryvariant = null;
			for (VCFVariant var : tmpvariants) {
				if (var.getChrObj().equals(variant.getChromosome())) {
					if (var.getPos() == variant.getStart()) {
						queryvariant = var;
					}
				}
			}
			if (queryvariant != null) {
				System.out.println(variant.toString() + " present in tabix.");
			}

			if (queryvariant != null) {

				if (!Double.isNaN(ldthreshold)) {
					DetermineLD calc = new DetermineLD();
					for (VCFVariant var : tmpvariants) {
						Pair<Double, Double> ld = calc.getLD(queryvariant, var);
						double rsq = ld.getRight();
						if (rsq > ldthreshold) {
							proxies.add(var);
							ldvals.add(rsq);
						}
					}
				} else {
					proxies.add(queryvariant);
				}
				System.out.println(proxies.size() + " proxyfinder in tabix");
			} else {
				System.out.println(variant.toString() + " not found in reference");
				proxies = null;
			}

			output.add(new Triple<SNPFeature, ArrayList<VCFVariant>, ArrayList<Double>>(variant, proxies, ldvals));
		}
		return output;
	}

}
