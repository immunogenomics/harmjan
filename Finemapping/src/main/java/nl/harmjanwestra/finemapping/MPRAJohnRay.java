package nl.harmjanwestra.finemapping;

import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 8/8/16.
 */
public class MPRAJohnRay {


	public static void main(String[] args) {
		MPRAJohnRay r = new MPRAJohnRay();
		try {
//			r.filterResults();
			r.calculateProxies();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void calculateProxies() throws IOException {

		String tabixrefprefix = "/Data/Ref/beagle_1kg/1kg.phase3.v5a.chr1.vcf.gz";
		String outfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-08-08-JohnRayMPRA/filterout/ldfile.txt";
		String raybedfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-08-08-JohnRayMPRA/20160617_both_sig_MPRA_SNPs.bed";

		//
		Feature region = new Feature();
		region.setStart(137882875);
		region.setStop(138275085);
		region.setChromosome(Chromosome.SIX);

		// 6_137999562_rs397886367
		Feature qsnp1 = new Feature();
		qsnp1.setChromosome(Chromosome.SIX);
		qsnp1.setStart(137999562);
		qsnp1.setStop(137999563);
		qsnp1.setName("rs397886367");

		// 6_138243739_rs58721818
		Feature qsnp2 = new Feature();
		qsnp2.setChromosome(Chromosome.SIX);
		qsnp2.setStart(138243739);
		qsnp2.setStop(138243739);
		qsnp2.setName("rs58721818");

		ArrayList<Feature> bedregions = readBed(raybedfile);
		System.out.println(bedregions.size() + " ray bed regions");
		ArrayList<VCFVariant> variants = new ArrayList<>();
		for (Feature f : bedregions) {
			String tabixfile = tabixrefprefix + f.getChromosome().getNumber() + ".vcf.gz";
			TabixReader reader = new TabixReader(tabixfile);
			TabixReader.Iterator window = reader.query(f.getChromosome().getNumber() + ":" + (f.getStart() - 10) + "-" + (f.getStop() + 10));
			String next = window.next();
			VCFVariant variant1 = null;

			int nr = 0;
			boolean found = false;
			while (next != null) {
				VCFVariant variant = new VCFVariant(next, VCFVariant.PARSE.ALL);
				if (variant.asFeature().overlaps(f)) {
					if (variant.getId().equals(f.getName())) {
						variant1 = new VCFVariant(next, VCFVariant.PARSE.ALL);
						variants.add(variant1);
						found = true;
					}
				}
				next = window.next();
				nr++;
			}

			reader.close();

			if (!found) {
				System.out.println("Could not find variant: " + f.getName());
			}
		}
		System.out.println(variants.size() + " variants loaded");


		VCFVariant qsnp1obj = null;
		VCFVariant qsnp2obj = null;

		String tabixfile = tabixrefprefix + qsnp1.getChromosome().getNumber() + ".vcf.gz";
		TabixReader reader = new TabixReader(tabixfile);
		TabixReader.Iterator window = reader.query(qsnp1.getChromosome().getNumber() + ":" + (region.getStart() - 10) + "-" + (region.getStop() + 10));


		String next = window.next();
		int nrTotal = 0;
		while (next != null) {
			VCFVariant variant = new VCFVariant(next, VCFVariant.PARSE.ALL);

			if (variant.getPos() == qsnp1.getStart()) {
				qsnp1obj = variant;
			}
			if (variant.getPos() == qsnp2.getStart()) {
				qsnp2obj = variant;
			}

			next = window.next();
			nrTotal++;
		}
		reader.close();

		System.out.println(nrTotal + " variants in region.");


		DetermineLD ld = new DetermineLD();


		if (qsnp1obj == null || qsnp2obj == null) {
			if (qsnp1obj == null) {
				System.out.println("Could not find obj1: " + qsnp1.getName());
				System.out.println(qsnp1.getChromosome().toString() + "_" + qsnp1.getStart() + "_" + qsnp1.getName());
			}
			if (qsnp2obj == null) {
				System.out.println("Could not find obj2: " + qsnp2.getName());
				System.out.println(qsnp2.getChromosome().toString() + "_" + qsnp2.getStart() + "_" + qsnp2.getName());
			}
			System.exit(-1);
		} else {
			System.out.println("Found both snps");
		}

		TextFile out = new TextFile(outfile, TextFile.W);
		out.writeln("Chromosome\tPosition\tId\tDistanceToSignal1\tRSquaredSignal1\tDistanceToSignal2\tRSquaredSignal2");
		for (int i = 0; i < variants.size(); i++) {

			VCFVariant var = variants.get(i);

			Pair<Double, Double> ld1 = ld.getLD(qsnp1obj, var);
			Pair<Double, Double> ld2 = ld.getLD(qsnp2obj, var);

			out.writeln(var.getChrObj().toString()
					+ "\t" + var.getPos()
					+ "\t" + var.getId()
					+ "\t" + (var.getPos() - qsnp1obj.getPos())
					+ "\t" + ld1.getRight()
					+ "\t" + (var.getPos() - qsnp2obj.getPos())
					+ "\t" + ld2.getRight()
			);

		}
		out.close();

	}

	public void filterResults() throws IOException {
		String raresults = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/RA-assoc0.3-COSMO-merged-posterior.txt.gz";
		raresults = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ImmunoBase/hg19_gwas_ra_okada_4_19_1.tab.gz";
		raresults = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ImmunoBase/hg19_gwas_ic_ra_eyre_4_19_1.tab.gz";
		raresults = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-08-08-JohnRayMPRA/RA-assoc0.3-COSMO-chr6-gwas-1.txt";
		String raybedfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-08-08-JohnRayMPRA/20160617_both_sig_MPRA_SNPs.bed";
		String tabixprefix = "/Data/Ref/1kg/cosmo.1kg.phase3.v5.chr";
		String outSharedhits = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-08-08-JohnRayMPRA/filterout/shared-cond1.txt";
		String outNonShared = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-08-08-JohnRayMPRA/filterout/nonshared-cond1.txt";
		String outProxies = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-08-08-JohnRayMPRA/filterout/proxiesshared-cond1.txt";

		int windowsize = 250000;
		double ldthreshold = 0.8;

		Feature region = new Feature();
		region.setStart(137882875);
		region.setStop(138275085);
		region.setChromosome(Chromosome.SIX);


		AssociationFile f = new AssociationFile();
		ArrayList<AssociationResult> results = f.read(raresults, region);

		System.out.println(results.size() + " assoc results");

		ArrayList<Feature> bedregions = readBed(raybedfile);

		System.out.println(bedregions.size() + " ray bed regions");


		TextFile out1 = new TextFile(outSharedhits, TextFile.W);
		TextFile out2 = new TextFile(outNonShared, TextFile.W);
		TextFile out3 = new TextFile(outProxies, TextFile.W);
		out1.writeln(f.getHeader());
		out2.writeln(f.getHeader());
		out3.writeln(f.getHeader());
		int shared = 0;
		for (AssociationResult r : results) {
			boolean overlap = false;
			for (Feature raysnp : bedregions) {
				if (raysnp.overlaps(r.getSnp())) {
					overlap = true;
				}
			}

			if (overlap) {
				out1.writeln(r.toString());
				shared++;
			} else {
				out2.writeln(r.toString());
			}
		}

		out1.close();
		out2.close();
		System.out.println(shared + " shared results");

		// get LD for all variants in Ray's list
		ArrayList<ArrayList<Feature>> proxies = getProxyVariants(tabixprefix, bedregions, windowsize, ldthreshold);
		shared = 0;
		for (AssociationResult r : results) {
			boolean overlap = false;
			for (int q = 0; q < proxies.size(); q++) {
				ArrayList<Feature> prox = proxies.get(q);
				if (prox != null) {
					for (Feature raysnp : prox) {
						if (raysnp.overlaps(r.getSnp())) {
							overlap = true;
							shared++;
						}
					}
				}
			}

			if (overlap) {
				out3.writeln(r.toString());
			}
		}
		out3.close();
		System.out.println(shared + " shared results");

	}

	private ArrayList<Feature> readBed(String raybedfile) throws IOException {
		ArrayList<Feature> output = new ArrayList<>();
		TextFile tf = new TextFile(raybedfile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length > 3) {
				Feature f = new Feature();
				f.setChromosome(Chromosome.parseChr(elems[0]));
				f.setStart(Integer.parseInt(elems[1]));
				f.setStop(Integer.parseInt(elems[2]));
				f.setName(elems[3]);
				output.add(f);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		return output;
	}

	private ArrayList<ArrayList<Feature>> getProxyVariants(String tabixrefprefix, ArrayList<Feature> variants, int windowSize, double ldthreshold) throws IOException {

		ArrayList<ArrayList<Feature>> output = new ArrayList<>();
		DetermineLD ld = new DetermineLD();
		for (int i = 0; i < variants.size(); i++) {
			Feature snp1 = variants.get(i);
			System.out.println("Getting proxyfinder for: " + i + "/" + variants.size() + " - " + snp1.getName());
			String tabixfile = tabixrefprefix + snp1.getChromosome().getNumber() + ".vcf.gz";
			TabixReader reader = new TabixReader(tabixfile);
			TabixReader.Iterator window = reader.query(snp1.getChromosome().getNumber() + ":" + (snp1.getStart() - 10) + "-" + (snp1.getStop() + 10));
			String next = window.next();
			VCFVariant variant1 = null;


			int nr = 0;
			while (next != null) {
				VCFVariant variant = new VCFVariant(next, VCFVariant.PARSE.ALL);
				if (variant.asFeature().overlaps(snp1)) {
					if (variant.getId().equals(snp1.getName())) {
						variant1 = new VCFVariant(next, VCFVariant.PARSE.ALL);
					}
				}
				next = window.next();
				nr++;
			}
			reader.close();

			if (variant1 == null) {
				System.out.println(snp1.getName() + " not found in reference");
			}

			if (variant1 != null) {
				ArrayList<Feature> snpout = new ArrayList<>();

				reader = new TabixReader(tabixfile);
				window = reader.query(snp1.getChromosome().getNumber() + ":" + (snp1.getStart() - windowSize) + "-" + (snp1.getStop() + windowSize));
				next = window.next();

				int nrproxies = 0;
				int nrvariants = 0;
				while (next != null) {
					VCFVariant variant = new VCFVariant(next, VCFVariant.PARSE.ALL);

					// calculate LD
					Pair<Double, Double> ldvals = ld.getLD(variant1, variant);
					double rsq = ldvals.getRight();
					if (rsq > ldthreshold) {
						snpout.add(variant.asFeature());
						nrproxies++;
					}

					next = window.next();
					nrvariants++;
					nr++;
				}
				reader.close();
				System.out.println(nrproxies + "/" + nrvariants + " are proxyfinder for this snp");
				output.add(snpout);
			} else {
				output.add(null);
			}
		}

		return output;
	}

}
