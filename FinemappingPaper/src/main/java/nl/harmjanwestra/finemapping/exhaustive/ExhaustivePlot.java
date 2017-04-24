package nl.harmjanwestra.finemapping.exhaustive;

import com.itextpdf.text.DocumentException;
import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.panels.HeatmapPanel;
import nl.harmjanwestra.utilities.graphics.panels.LDPanel;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFTabix;
import nl.harmjanwestra.utilities.vcf.VCFVariant;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 8/30/16.
 */
public class ExhaustivePlot {


	public static void main(String[] args) {

		try {
			System.out.println("works");
			BedFileReader reader = new BedFileReader();
			ArrayList<Feature> regions = reader.readAsList("/Data/Projects/2016-Finemapping/Exhaustive/data/2017-03-28-RegionsExhaustive.txt");
//			String[] diseases = new String[]{"RA", "T1D"};
			String[] diseases = new String[]{"RA", "T1D", "META"};
//			diseases = new String[]{"RA"};
			int nrToPlot = 25;
			for (String d : diseases) {
				for (Feature region : regions) {
//					if (region.getChromosome().equals(Chromosome.TEN)) {
					ExhaustivePlot p = new ExhaustivePlot();
					String output = "/Data/Projects/2016-Finemapping/Exhaustive/output/" + d + "-tyk2-" + region.toString();

					String assocfile = "/Data/Projects/2016-Finemapping/Exhaustive/data/" + d + "-assoc0.3-COSMO-tyk2-chr" + region.getChromosome().getNumber() + "-pairwise-rewrite.txt.gz";
					String tabixfile = "/Data/Ref/beagle_1kg/1kg.phase3.v5a.chrCHR.vcf.gz";
					String tabixsamplelimit = "/Data/Ref/1kg-europeanpopulations.txt.gz";

					String snpcombos = "/Data/Projects/2016-Finemapping/Exhaustive/snpcombos.txt";
//					if (region.getChromosome().equals(Chromosome.SIX)) {

					if (Gpio.exists(assocfile)) {
						p.findCombination(region, assocfile, snpcombos, tabixfile, tabixsamplelimit);
//							p.plotTopNSNPs(region, assocfile, output, nrToPlot, tabixfile, tabixsamplelimit);
					}
//					}
//					}


				}
			}


		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}

	}

	private void findCombination(Feature region, String assocfile, String snpcombos, String tabixPrefix, String tabixsamplelimit) throws IOException, DocumentException {

		System.out.println("Trying to find snp combos in region: " + region.toString());
		System.out.println("Assoc: " + assocfile);

		TextFile tfc = new TextFile(snpcombos, TextFile.R);

		String[] elems = tfc.readLineElems(TextFile.tab);

		SNPFeature snpf1 = null;
		SNPFeature snpf2 = null;

		while (elems != null) {

			if (elems.length >= 3) {
				Feature region1 = Feature.parseFeature(elems[0]);
				SNPFeature feature1 = SNPFeature.parseSNPFeature(elems[1]);
				SNPFeature feature2 = SNPFeature.parseSNPFeature(elems[2]);
//				System.out.println(region1.toString() + "\t" + feature1.toString());
//			System.out.println(region1.toString() + "\tLooking for snp: " + feature1.toString() + " and " + feature2.toString());
				if (region1.overlaps(region)) {
					snpf1 = feature1;
					snpf2 = feature2;
				}
			}
			elems = tfc.readLineElems(TextFile.tab);
		}
		tfc.close();


		if (snpf1 == null || snpf2 == null) {
			System.out.println("Could not read snp pair...");
		} else {
			System.out.println(region.toString() + "\tLooking for snp: " + snpf1.toString() + " and " + snpf2.toString());

			String[] variants = new String[]{"", ""};

			TextFile tf2 = new TextFile(assocfile, TextFile.R);
			String headerln = tf2.readLine();
			elems = tf2.readLineElems(TextFile.tab);
			double maxP = 0;

			// snp index
			HashMap<SNPFeature, Integer> snpToIn = new HashMap<SNPFeature, Integer>();
			ArrayList<SNPFeature> allSNPs = new ArrayList<SNPFeature>();

			ArrayList<AssocPair> assocVals = new ArrayList<>();

			AssocPair select = null;

			while (elems != null) {
//			lineIsWhatWereLookingFor = false;
				Feature s1 = new Feature();
				Chromosome chr = Chromosome.parseChr(elems[0]);
				Integer pos1 = Integer.parseInt(elems[1]);
				String snp1Id = elems[2];
				Integer pos2 = Integer.parseInt(elems[5]);
				String snp2Id = elems[6];

				SNPFeature snp1 = new SNPFeature();
				snp1.setChromosome(chr);
				snp1.setStart(pos1);
				snp1.setStop(pos1);
				snp1.setName(snp1Id);

				if (!snpToIn.containsKey(snp1)) {
					snpToIn.put(snp1, allSNPs.size());
					allSNPs.add(snp1);
				} else {
					Integer id = snpToIn.get(snp1);
					if (id == null) {
						System.out.println("meh");
					}
					snp1 = allSNPs.get(id);
				}

				SNPFeature snp2 = new SNPFeature();
				snp2.setChromosome(chr);
				snp2.setStart(pos2);
				snp2.setStop(pos2);
				snp2.setName(snp2Id);

				if (!snpToIn.containsKey(snp2)) {
					snpToIn.put(snp2, allSNPs.size());
					allSNPs.add(snp2);
				} else {
					Integer id = snpToIn.get(snp2);
					snp2 = allSNPs.get(id);
				}

				if (region.overlaps(snp1) && region.overlaps(snp2)) {
					String snp1str = snp1.getChromosome().getNumber() + "_" + snp1.getStart() + "_" + snp1.getName();
					String snp2str = snp2.getChromosome().getNumber() + "_" + snp2.getStart() + "_" + snp2.getName();
//					System.out.println(snp1str);
//					System.exit(0);
					Double p = Double.parseDouble(elems[elems.length - 1]);
					if ((snp1str.equals(variants[0]) && snp2str.equals(variants[1]))
							|| (snp1str.equals(variants[1]) && snp2str.equals(variants[0]))) {
						System.out.println("found it.");
					}

					if (p > maxP) {
						maxP = p;
					}

					AssocPair v = new AssocPair(snp1, snp2, p);
					assocVals.add(v);

					if ((snp1.equals(snpf1) && snp2.equals(snpf2))
							|| (snp2.equals(snpf1) && snp1.equals(snpf2))
							) {
						select = v;
					}

				}
				elems = tf2.readLineElems(TextFile.tab);
			}
			tf2.close();


			if (select == null) {
				System.out.println("Could not find pair");
			} else {
				int nrsmaller = 0;
				System.out.println("P-value pair selected: " + select.getP());
				AssocPair top = null;
				double maxp = 0;
				for (AssocPair v : assocVals) {
					if (v.getP() > maxp) {
						maxp = v.getP();
						top = v;
					}
					if (v.getP() > select.getP()) {
						nrsmaller++;
					}
				}

				System.out.println("Top pair: " + top.getSnp1().toString() + "\t" + top.getSnp2().toString() + "\t" + top.getP());
				// calculate LD with selected pair...?
				// get the LD for these variants

				SNPFeature[] snplist = new SNPFeature[4];
				snplist[0] = select.getSnp1();
				snplist[1] = select.getSnp2();
				snplist[2] = top.getSnp1();
				snplist[3] = top.getSnp2();


				String tabixfile = tabixPrefix.replaceAll("CHR", "" + region.getChromosome().getNumber());
				VCFTabix tabix = new VCFTabix(tabixfile);
				boolean[] includeSamples = null;
				if (tabixsamplelimit != null) {
					includeSamples = tabix.getSampleFilter(tabixsamplelimit);
				}

				TabixReader.Iterator window = tabix.query(region);
				String next = window.next();

				VCFVariant[] vars = new VCFVariant[4];
				while (next != null) {
					VCFVariant variant = new VCFVariant(next, VCFVariant.PARSE.HEADER);

					for (int v = 0; v < snplist.length; v++) {
						if (variant.asFeature().overlaps(snplist[v])) {
							vars[v] = new VCFVariant(next, VCFVariant.PARSE.ALL, includeSamples);
						}
					}
					next = window.next();
				}
				tabix.close();

				// calculate LD
				DetermineLD ld = new DetermineLD();

				for (int v = 0; v < 2; v++) {
					for (int v2 = 2; v2 < snplist.length; v2++) {
						Pair<Double, Double> ldvals = ld.getLD(vars[v], vars[v2]);
						System.out.println(snplist[v] + " - " + snplist[v2] + "\tr2: " + ldvals.getRight() + "\tdprime: " + ldvals.getLeft());
					}
				}

				System.out.println(region.toString() + "\t" + nrsmaller + " more significant values out of " + assocVals.size());
			}
			System.out.println();
			System.out.println();
		}


	}


	public static void plotFullHalfMatrix(Feature region, String assocfile, String output) throws IOException, DocumentException {

		String[] variants = new String[]{"", ""};

		TextFile tf2 = new TextFile(assocfile, TextFile.R);
		String headerln = tf2.readLine();
		String[] elems = tf2.readLineElems(TextFile.tab);
		double maxP = 0;
		ArrayList<Pair<Integer, Integer>> positions = new ArrayList<Pair<Integer, Integer>>();
		ArrayList<Double> pvals = new ArrayList<>();
		while (elems != null) {
//			lineIsWhatWereLookingFor = false;
			Feature s1 = new Feature();
			Chromosome chr = Chromosome.parseChr(elems[0]);
			Integer pos1 = Integer.parseInt(elems[1]);
			String snp1Id = elems[2];
			Integer pos2 = Integer.parseInt(elems[5]);
			String snp2Id = elems[6];

			SNPFeature snp1 = new SNPFeature();
			snp1.setChromosome(chr);
			snp1.setStart(pos1);
			snp1.setStop(pos1);
			snp1.setName(snp1Id);

			SNPFeature snp2 = new SNPFeature();
			snp2.setChromosome(chr);
			snp2.setStart(pos2);
			snp2.setStop(pos2);
			snp2.setName(snp2Id);


			if (region.overlaps(snp1) && region.overlaps(snp2)) {


				String snp1str = snp1.getChromosome().getNumber() + "_" + snp1.getStart() + "_" + snp1.getName();
				String snp2str = snp2.getChromosome().getNumber() + "_" + snp2.getStart() + "_" + snp2.getName();
//					System.out.println(snp1str);
//					System.exit(0);
				Double p = Double.parseDouble(elems[elems.length - 1]);
				if ((snp1str.equals(variants[0]) && snp2str.equals(variants[1]))
						|| (snp1str.equals(variants[1]) && snp2str.equals(variants[0]))) {
					System.out.println("found it.");

				}

				if (p > maxP) {
					maxP = p;
				}

				positions.add(new Pair<>(snp1.getStart(), snp2.getStart()));
				pvals.add(p);


			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();

		System.out.println(pvals.size() + " pvals");


		ArrayList<Double> pvals2 = new ArrayList<>();
		for (int i = 0; i < pvals.size(); i++) {
			pvals2.add(pvals.get(i) / maxP);
		}

		// make the actual plot
		Grid grid = new Grid(400, 400, 1, 1, 100, 100);

		LDPanel panel = new LDPanel(1, 1);
		panel.setData(region, positions, pvals2);
		grid.addPanel(panel);
		panel.scaleQuadratic(true);
		grid.draw(output + ".png");
	}


	public void plotTopNSNPs(Feature region, String assocfile, String output, int topN, String tabixPrefix, String tabixsamplelimit) throws IOException, DocumentException {

		String[] variants = new String[]{"", ""};

		TextFile tf2 = new TextFile(assocfile, TextFile.R);
		String headerln = tf2.readLine();
		String[] elems = tf2.readLineElems(TextFile.tab);
		double maxP = 0;

		// snp index
		HashMap<SNPFeature, Integer> snpToIn = new HashMap<SNPFeature, Integer>();
		ArrayList<SNPFeature> allSNPs = new ArrayList<SNPFeature>();

		ArrayList<AssocPair> assocVals = new ArrayList<>();

		while (elems != null) {
//			lineIsWhatWereLookingFor = false;
			Feature s1 = new Feature();
			Chromosome chr = Chromosome.parseChr(elems[0]);
			Integer pos1 = Integer.parseInt(elems[1]);
			String snp1Id = elems[2];
			Integer pos2 = Integer.parseInt(elems[5]);
			String snp2Id = elems[6];

			SNPFeature snp1 = new SNPFeature();
			snp1.setChromosome(chr);
			snp1.setStart(pos1);
			snp1.setStop(pos1);
			snp1.setName(snp1Id);

			if (!snpToIn.containsKey(snp1)) {
				snpToIn.put(snp1, allSNPs.size());
				allSNPs.add(snp1);
			} else {
				Integer id = snpToIn.get(snp1);
				if (id == null) {
					System.out.println("meh");
				}
				snp1 = allSNPs.get(id);
			}

			SNPFeature snp2 = new SNPFeature();
			snp2.setChromosome(chr);
			snp2.setStart(pos2);
			snp2.setStop(pos2);
			snp2.setName(snp2Id);

			if (!snpToIn.containsKey(snp2)) {
				snpToIn.put(snp2, allSNPs.size());
				allSNPs.add(snp2);
			} else {
				Integer id = snpToIn.get(snp2);
				snp2 = allSNPs.get(id);
			}

			if (region.overlaps(snp1) && region.overlaps(snp2)) {
				String snp1str = snp1.getChromosome().getNumber() + "_" + snp1.getStart() + "_" + snp1.getName();
				String snp2str = snp2.getChromosome().getNumber() + "_" + snp2.getStart() + "_" + snp2.getName();
//					System.out.println(snp1str);
//					System.exit(0);
				Double p = Double.parseDouble(elems[elems.length - 1]);
				if ((snp1str.equals(variants[0]) && snp2str.equals(variants[1]))
						|| (snp1str.equals(variants[1]) && snp2str.equals(variants[0]))) {
					System.out.println("found it.");
				}

				if (p > maxP) {
					maxP = p;
				}

				assocVals.add(new AssocPair(snp1, snp2, p));

			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();

		Collections.sort(assocVals);
		System.out.println(assocVals.get(0).getSnp1() + "\t" + assocVals.get(0).getSnp2() + "\t" + assocVals.get(0).getP());

		System.out.println(assocVals.size() + " values loaded");

		HashSet<SNPFeature> topSNPs = new HashSet<SNPFeature>();
		int ctr = 0;
		while (topSNPs.size() < topN) {
			AssocPair p = assocVals.get(ctr);
			topSNPs.add(p.getSnp1());
			topSNPs.add(p.getSnp2());
			ctr++;
		}

		// sort top SNPs by position
		ArrayList<SNPFeature> allSNPFeatures = new ArrayList<SNPFeature>();
		allSNPFeatures.addAll(topSNPs);
		Collections.sort(allSNPFeatures, new FeatureComparator(true));

		double[][] vals = new double[allSNPFeatures.size()][allSNPFeatures.size()];
		// fill with NaNs
		for (int i = 0; i < vals.length; i++) {
			for (int j = i; j < vals.length; j++) {
				vals[i][j] = Double.NaN;
				vals[j][i] = Double.NaN;
			}
		}

		HashMap<SNPFeature, Integer> topSNPIndex = new HashMap<SNPFeature, Integer>();
		HashMap<Integer, Integer> topSNPPosIndex = new HashMap<>();
		String[] rowlabels = new String[allSNPFeatures.size()];
		for (int q = 0; q < allSNPFeatures.size(); q++) {
			topSNPIndex.put(allSNPFeatures.get(q), q);
			topSNPPosIndex.put(allSNPFeatures.get(q).getStart(), q);
			rowlabels[q] = allSNPFeatures.get(q).getName();
		}
		// iterate through all assoc PVals to get pairs corresponding to snp set
		double maxPval = 0;
		double thresh = -Math.log10(7.5e-7);
		double minPval = Double.MAX_VALUE;
		for (AssocPair p : assocVals) {
			if (topSNPs.contains(p.getSnp1()) && topSNPs.contains(p.getSnp2())) {
				Integer id1 = topSNPIndex.get(p.getSnp1());
				Integer id2 = topSNPIndex.get(p.getSnp2());

				if (p.getP() > thresh) {
					vals[id1][id2] = p.getP();
					vals[id2][id1] = p.getP();
					if (p.getP() > maxPval) {
						maxPval = p.getP();
					}
					if (p.getP() < minPval) {
						minPval = p.getP();
					}
				} else {
					vals[id1][id2] = Double.NaN;
					vals[id2][id1] = Double.NaN;
				}
			}
		}


		// get the LD for these variants
		String tabixfile = tabixPrefix.replaceAll("CHR", "" + region.getChromosome().getNumber());
		VCFTabix tabix = new VCFTabix(tabixfile);
		boolean[] includeSamples = null;
		if (tabixsamplelimit != null) {
			includeSamples = tabix.getSampleFilter(tabixsamplelimit);
		}

		TabixReader.Iterator window = tabix.query(region);
		String next = window.next();
		VCFVariant[] vcfVariants = new VCFVariant[topSNPs.size()];
		while (next != null) {
			VCFVariant variant = new VCFVariant(next, VCFVariant.PARSE.HEADER);
			if (topSNPPosIndex.containsKey(variant.getPos())) {
				Integer id = topSNPPosIndex.get(variant.getPos());
				vcfVariants[id] = new VCFVariant(next, VCFVariant.PARSE.ALL, includeSamples);
			}
			next = window.next();
		}
		tabix.close();
		double[][] ldmatrix = new double[topSNPs.size()][topSNPs.size()];
		DetermineLD ldcalc = new DetermineLD();
		for (int i = 0; i < vcfVariants.length; i++) {
			for (int j = i + 1; j < vcfVariants.length; j++) {
				ldmatrix[i][j] = ldcalc.getLD(vcfVariants[i], vcfVariants[j]).getRight();
				ldmatrix[j][i] = ldmatrix[i][j];
			}
		}

		if (vals.length != 0) {

			// merge LD matrix into pval matrix
			for (int i = 0; i < vals.length; i++) {
				for (int j = i; j < vals.length; j++) {
					if (i != j) {
						vals[i][j] = Double.NaN;
						vals[i][j] = ldmatrix[i][j];
					}
				}
			}


			System.out.println("Max P: " + maxPval);
			double remainder = maxPval % 1;
			maxPval += (1 - remainder);
			System.out.println("Max P: " + maxPval);
			System.out.println("Min P: " + minPval);
			remainder = minPval % 1;
			minPval -= remainder;
			System.out.println("Min P: " + minPval);

			Range pvalRange = new Range(0, minPval, 1, maxPval);
			Range ldRange = new Range(0, 0, 1, 1);

			Grid grid = new Grid(600, 600, 1, 1, 100, 100);
			HeatmapPanel panel = new HeatmapPanel(1, 1);
			panel.setData(vals, rowlabels, rowlabels);
			panel.setPlotMode(HeatmapPanel.MODE.TWODS);
			grid.addPanel(panel);
			panel.setRangeLower(pvalRange);
			panel.setRangeUpper(ldRange);
			grid.draw(output + "-pvals.pdf");

//			grid = new Grid(600, 600, 1, 1, 100, 100);
//			panel = new HeatmapPanel(1, 1);
//			panel.setData(ldmatrix, rowlabels, rowlabels);
//			panel.setPlotMode(HeatmapPanel.MODE.UPPER);
//			grid.addPanel(panel);
//			grid.draw(output + "-ld.pdf");
//			TextFile tf = new TextFile(output + "mat.txt", TextFile.W);
//			for (int i = 0; i < vals.length; i++) {
//				tf.writeln(Strings.concat(vals[i], Strings.tab));
//			}
//			tf.close();
		} else {
			System.out.println("Nothing to plot?");
		}
	}


}

class AssocPair implements Comparable<AssocPair> {

	private SNPFeature snp1;
	private SNPFeature snp2;

	double p;

	public AssocPair(SNPFeature snp1, SNPFeature snp2, Double p) {
		this.p = p;
		this.snp1 = snp1;
		this.snp2 = snp2;
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (o == null || getClass() != o.getClass()) return false;
		AssocPair assocPair = (AssocPair) o;
		return Double.compare(assocPair.p, p) == 0;
	}

	@Override
	public int hashCode() {
		long temp = Double.doubleToLongBits(p);
		return (int) (temp ^ (temp >>> 32));
	}

	@Override
	public int compareTo(AssocPair o) {
		if (o.p == p) {
			return 0;
		} else if (p > o.p) {
			return -1;
		} else {
			return 1;
		}
	}

	public SNPFeature getSnp1() {
		return snp1;
	}

	public SNPFeature getSnp2() {
		return snp2;
	}

	public double getP() {
		return p;
	}
}