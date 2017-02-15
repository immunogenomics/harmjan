package nl.harmjanwestra.finemapping;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFTabix;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by hwestra on 10/28/16.
 */
public class eQTLOverlap {


	public void run(String regionvariantfile,
					String tabixprefix,
					String limitsamples,
					String[] eQTLfiles,
					double ldthreshold,
					String output) throws IOException {
		// load variants and regions

		ArrayList<Pair<Feature, SNPFeature>> variantsPerRegion = new ArrayList<Pair<Feature, SNPFeature>>();
		TextFile tf = new TextFile(regionvariantfile, TextFile.R);
		HashSet<Feature> uniqueRegions = new HashSet<Feature>();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
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


			Feature f = Feature.parseFeature(regionStr);
			variantsPerRegion.add(new Pair<Feature, SNPFeature>(f, feat));

			uniqueRegions.add(f);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();


		ArrayList<ArrayList<EQTL>> alleQTL = loadEQTLFiles(eQTLfiles, uniqueRegions);

		VCFTabix tabix = new VCFTabix(tabixprefix);
		boolean[] samplefilter = tabix.getSampleFilter(limitsamples);
		DetermineLD ldcalc = new DetermineLD();
		// for each region
		TextFile outtf = new TextFile(output, TextFile.W);
		outtf.writeln("Region\tVariant\teSNP\teQTLGene\tLD");
		for (Feature region : uniqueRegions) {
			// load all 1kg variants in region
			ArrayList<VCFVariant> vcfVariants = tabix.getAllVariants(region, samplefilter);

			ArrayList<ArrayList<EQTL>> eQTLInRegion = filterEQTL(alleQTL, region);

			for (Pair<Feature, SNPFeature> queryPair : variantsPerRegion) {

				if (queryPair.getLeft().overlaps(region)) {
					SNPFeature variantFeat = queryPair.getRight();

					VCFVariant variant = getVariant(vcfVariants, variantFeat);

					// load all eQTLs that have a variant in the region
					for (ArrayList<EQTL> eqtlList : eQTLInRegion) {


						// for each variant in the query set
						for (EQTL e : eqtlList) {
							VCFVariant evariant = getVariant(vcfVariants, e.getSnp());
							// check LD between variants
							Pair<Double, Double> ld = ldcalc.getLD(variant, evariant);
							double rsq = ld.getRight();
							// output
							if (rsq > ldthreshold) {
								outtf.writeln(region.toString() + "\t" + variant.toString() + "\t" + evariant.toString() + "\t" + e.getGenename() + "\t" + rsq);
							}
						}
					}
				}
			}
		}
		outtf.close();
	}

	private VCFVariant getVariant(ArrayList<VCFVariant> vcfVariants, SNPFeature variantFeat) {

		for (VCFVariant v : vcfVariants) {
			if (v.getChrObj().equals(variantFeat.getChromosome())) {
				if (v.getPos() == variantFeat.getStart()) {
					return v;
				}
			}
		}
		return null;
	}

	private ArrayList<ArrayList<EQTL>> filterEQTL(ArrayList<ArrayList<EQTL>> alleQTL, Feature region) {

		ArrayList<ArrayList<EQTL>> output = new ArrayList<>();
		for (ArrayList<EQTL> a : alleQTL) {
			ArrayList<EQTL> out = new ArrayList<>();
			for (EQTL e : a) {
				if (e.overlaps(region)) {
					out.add(e);
				}
			}
			output.add(out);
		}

		return output;
	}

	private ArrayList<ArrayList<EQTL>> loadEQTLFiles(String[] eQTLfiles, HashSet<Feature> uniqueRegions) throws IOException {


		ArrayList<ArrayList<EQTL>> eQTLs = new ArrayList<>();
		for (int e = 0; e < eQTLfiles.length; e++) {
			String f = eQTLfiles[e];
			TextFile tf = new TextFile(f, TextFile.R);

			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {

				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();


		}
		return eQTLs;
	}

}
