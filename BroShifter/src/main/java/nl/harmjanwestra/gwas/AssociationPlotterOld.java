package nl.harmjanwestra.gwas;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.broshifter.AnnotationLoader;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.PosteriorPValFile;
import nl.harmjanwestra.utilities.association.VCFRSquares;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.*;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.*;
import nl.harmjanwestra.utilities.gtf.GTFAnnotation;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Created by hwestra on 11/23/15.
 */
public class AssociationPlotterOld {


	public static void main(String[] args) {
		AssociationPlotterOld p = new AssociationPlotterOld();
		try {

//			Feature queryVariantFeature = new Feature(Chromosome.TWO, 60462160, 62303311);
//			p.bedgraphfilter(queryVariantFeature, "/Data/Projects/2015-cd4timelinepilot/BedGraphsSeqProj/samples.txt");
//			p.makePlotSpecificRegion();

			Feature feat = new Feature(Chromosome.SIX, 32627241 - 5000, 32634466 + 5000);

//			String out = "/Data/tmp/dqb1/rna/plotBinaryTrait-atac.pdf";
//			String bg = "/Data/tmp/dqb1/atac/allbedfiles.txt";

			String out = "/Data/tmp/dqb1/rna/plotBinaryTrait.pdf";
			String bg = "/Data/tmp/dqb1/rna/allfiles.txt";
			String gtf = "/Data/Ref/Annotation/UCSC/genes.gtf";
//			String gtf = "/Data/Ref/Annotation/Gencode/gencode.v19.annotation.gtf.gz";
			p.makePlotSpecificRegionReads(gtf, bg, out, feat);

//			p.makePlotDec3();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void bedgraphfilter(Feature region, String bedGraphListFile) throws IOException {

		TextFile tf = new TextFile(bedGraphListFile, TextFile.R);
		String[] bfiles = tf.readAsArray();
		tf.close();

		for (String file : bfiles) {
			TextFile in = new TextFile(file, TextFile.R);
			TextFile out = new TextFile(file + "-filtered.bedGraph.gz", TextFile.W);

			String[] elems = in.readLineElems(TextFile.tab);
			while (elems != null) {
				Chromosome chr = Chromosome.parseChr(elems[0]);
				if (chr.equals(region.getChromosome())) {
					Integer pos = Integer.parseInt(elems[1]);
					if (pos >= region.getStart() && pos <= region.getStop()) {
						out.writeln(Strings.concat(elems, Strings.tab));
					}
				}
				elems = in.readLineElems(TextFile.tab);
			}
			in.close();
			out.close();
		}
	}


	public void makePlotSpecificRegionReads(String gtf, String bedGraphListFile, String outdir, Feature region) throws IOException {
		GTFAnnotation annot = new GTFAnnotation(gtf);

		TreeSet<Gene> genes = annot.getGeneTree();
		Gene geneStart = new Gene("", region.getChromosome(), Strand.POS, region.getStart(), region.getStart());
		Gene geneStop = new Gene("", region.getChromosome(), Strand.POS, region.getStop(), region.getStop());
		SortedSet<Gene> overlappingGenes = genes.subSet(geneStart, true, geneStop, true);
		ArrayList<Gene> overlappingGenesList = new ArrayList<>();
		overlappingGenesList.addAll(overlappingGenes);

//
//		Gene gene1 = null;
//		Collection<Gene> allGenes = annot.getGenes();
//		for (Gene g : allGenes) {
//			if (g.getGeneId().equals("HLA-DQB1")) {
//				gene1 = g;
//			}
//		}
//
//
//		int start = gene1.getStart() - 5000;
//		int stop = gene1.getStop() + 5000;
//		region.setStart(start);
//		region.setStop(stop);
//
//		System.out.println(start);
//		System.out.println(stop);


		Grid grid = new Grid(1200, 300, 2, 1, 300, 300);

		// add
		TextFile tf = new TextFile(bedGraphListFile, TextFile.R);
		String[] bfiles = tf.readAsArray();
		tf.close();

		// plotBinaryTrait gene annotations
		GenePanel genePanel = new GenePanel(1, 1);
		genePanel.setData(region, overlappingGenesList);
		grid.addPanel(genePanel, 0, 0);


		ArrayList<ArrayList<BedGraphFeature>> data = readFeatures(region, bedGraphListFile);
		GraphAnnotationPanel gp = new GraphAnnotationPanel(1, 1);
		gp.setData(region, data, bfiles);
		grid.addPanel(gp, 1, 0);

		try {
			grid.draw(outdir);
		} catch (DocumentException e) {
			e.printStackTrace();
		}

	}

	public void makePlotSpecificRegionAnnotation(String gtf, String listOfAnnotationsFile, String bedGraphListFile, String outdir, Feature region) throws IOException {
//		String gtf = "/Data/Ref/Annotation/UCSC/genes.gtf";
		GTFAnnotation annot = new GTFAnnotation(gtf);

//		String listOfAnnotationsFile = "/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/Peaks/allatacpeaks.txt";
//
//		String bedGraphListFile = "/Data/Projects/2015-cd4timelinepilot/BedGraphsSeqProj/samples2.txt";
//		String outdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/PlotsWAnnotations/REL-peakcenter.pdf";
//		Gpio.createDir(outdir);

		boolean usepeakcenter = true;
		int bptoextendannotation = 150;
		boolean mergeoverlapping = true;


//		Gene gene1 = null;
//		Gene gene2 = null;
//		Collection<Gene> allGenes = annot.getGenes();
//
//		for (Gene g : allGenes) {
//			if (g.getGeneId().equals("FLJ16341")) {
//				gene1 = g;
//			}
//			if (g.getGeneId().equals("REL")) {
//				gene2 = g;
//			}
//		}
//
//
//		int start = gene1.getStop() - 1000;
//		int stop = gene2.getStart() + 1000;

		Grid grid = new Grid(400, 300, 3, 1, 300, 300);

//		Feature region = new Feature(gene1.getChromosome(), start, stop);

		System.out.println(region.toString());

		TreeSet<Gene> genes = annot.getGeneTree();
		Gene geneStart = new Gene("", region.getChromosome(), Strand.POS, region.getStart(), region.getStart());
		Gene geneStop = new Gene("", region.getChromosome(), Strand.POS, region.getStop(), region.getStop());
		SortedSet<Gene> overlappingGenes = genes.subSet(geneStart, true, geneStop, true);

		ArrayList<Gene> overlappingGenesList = new ArrayList<>();
		overlappingGenesList.addAll(overlappingGenes);

		// get the subregions that overlap with this region
		// make a subplot for the zoomregion

		// plotBinaryTrait gene annotations
		GenePanel genePanel = new GenePanel(1, 1);
		genePanel.setData(region, overlappingGenesList);
		grid.addPanel(genePanel, 0, 0);

		// add annotations
		ArrayList<Track> annotations = loadAnnotations(listOfAnnotationsFile, usepeakcenter, bptoextendannotation, mergeoverlapping);
		AnnotationPanel annotPanel = new AnnotationPanel(1, 1);
		annotPanel.setData(region, annotations);
		grid.addPanel(annotPanel, 2, 0);


		// assocpvals
		PosteriorPValFile assocfile = new PosteriorPValFile();
		String pvalfile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/Posteriors/RA/seq/Chr2_60462160-62303311.txt";

		Feature origRegion = new Feature(Chromosome.TWO, 60462160, 62303311);

		ApproximateBayesPosterior abp = new ApproximateBayesPosterior();
		System.out.println("in: " + pvalfile);
		ArrayList<AssociationResult> assoc = assocfile.readVariantPValues(pvalfile, origRegion);
		ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(assoc, 0.95);


		ArrayList<Feature> credibleSetSNPList = new ArrayList<>();

		int minPosCredSet = Integer.MAX_VALUE;
		int maxPosCredSet = 0;
		for (AssociationResult r : credibleSet) {
			Feature snp = r.getSnp();
			if (snp.getStart() < minPosCredSet) {
				minPosCredSet = snp.getStart();

			}
			if (snp.getStart() > maxPosCredSet) {
				maxPosCredSet = snp.getStart();
			}
			credibleSetSNPList.add(snp);
		}

		if (maxPosCredSet > region.getStop()) {
			maxPosCredSet = region.getStop();
		}
		if (minPosCredSet < region.getStart()) {
			minPosCredSet = region.getStart();
		}

		HashSet<AssociationResult> hashCredible = new HashSet<>();
		hashCredible.addAll(credibleSet);
		boolean[] mark = new boolean[assoc.size()];
		for (int i = 0; i < assoc.size(); i++) {
			AssociationResult r = assoc.get(i);
			if (hashCredible.contains(r)) {
				mark[i] = true;
			}
		}


		// plotBinaryTrait pvalues and posteriors
		Pair<boolean[], ArrayList<AssociationResult>> filtered = filterAssoc(region, assoc, mark);
		mark = filtered.getLeft();
		assoc = filtered.getRight();

		ArrayList<Pair<Integer, Double>> pvals = convertToPairs(assoc, true);
		ArrayList<SNPFeature> snpFeatures = convertToFeatures(assoc);
		double maxp = getMax(pvals);
		ArrayList<Pair<Integer, Double>> posteriors = convertToPairs(assoc, false);
		double maxposterior = getMax(posteriors);
		AssociationPanel pvalpanel = new AssociationPanel(1, 1);
		pvalpanel.setTitle("Pvalues");
		pvalpanel.setDataSingleDs(region, null, pvals, "Pvalues");
		pvalpanel.setMaxPVal(maxp);
		pvalpanel.setMarkDifferentColor(mark);
		grid.addPanel(pvalpanel, 1, 0);

		try {
			grid.draw(outdir);
		} catch (DocumentException e) {
			e.printStackTrace();
		}


	}

	public ArrayList<ArrayList<BedGraphFeature>> readFeatures(Feature region, String bedGraphListFile) throws IOException {

		TextFile tf = new TextFile(bedGraphListFile, TextFile.R);
		String[] bfiles = tf.readAsArray();
		tf.close();

		ArrayList<ArrayList<BedGraphFeature>> output = new ArrayList<>();
		for (String file : bfiles) {
			TextFile in = new TextFile(file, TextFile.R);
			TextFile out = new TextFile(file + "-filtered.bedGraph.gz", TextFile.W);
			ArrayList<BedGraphFeature> tmp = new ArrayList<>();
			String[] elems = in.readLineElems(TextFile.tab);
			while (elems != null) {
				Chromosome chr = Chromosome.parseChr(elems[0]);
				if (chr.equals(region.getChromosome())) {
					Integer pos = Integer.parseInt(elems[1]);
					if (pos >= region.getStart() && pos <= region.getStop()) {
						double val = Double.parseDouble(elems[3]);
						BedGraphFeature f = new BedGraphFeature(chr, pos, pos + 1);
						f.setValue(val);
						tmp.add(f);
					}
				}
				elems = in.readLineElems(TextFile.tab);
			}
			output.add(tmp);
			in.close();
			out.close();
		}
		return output;
	}

	public void makePlotDec3() throws IOException {


		String[] ds = new String[]{"T1D", "RA"};
		String gtf = "/Data/Ref/Annotation/UCSC/genes.gtf";
		GTFAnnotation annot = new GTFAnnotation(gtf);
		for (int d = 0; d < ds.length; d++) {
			String dsname = ds[d];
			String zoomregionfile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/2015-03-30-allRegionsWith50PercentOfSamplesAbove20x.bed";
			String assocdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/Posteriors/" + dsname + "/seq/";
			String regionfile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";

			String listOfAnnotationsFile = "/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/Peaks/allatacpeaks.txt";

			String outdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/PlotsWAnnotations/" + dsname + "/";
			Gpio.createDir(outdir);

			boolean usepeakcenter = true;
			int bptoextendannotation = 150;
			boolean mergeoverlapping = true;


			// load posteriors and association p-vals
			// plotBinaryTrait  [ p-values | posteriors ]

			// load annotations
			// plotBinaryTrait annotations underneath

			BedFileReader bf = new BedFileReader();
			ArrayList<Feature> regions = bf.readAsList(regionfile);
			ArrayList<Feature> zoomregions = bf.readAsList(zoomregionfile);


			HashSet<Feature> zoomregionHash = new HashSet<Feature>();
			zoomregionHash.addAll(zoomregions);

			ArrayList<Track> annotations = loadAnnotations(listOfAnnotationsFile, usepeakcenter, bptoextendannotation, mergeoverlapping);


			ApproximateBayesPosterior abp = new ApproximateBayesPosterior();

			// produce a plotBinaryTrait for each region
			for (Feature region : regions) {

				// make a plotBinaryTrait for the whole region
				PosteriorPValFile assocfile = new PosteriorPValFile();
				String pvalfile = assocdir + region.toString() + ".txt";
				if (!Gpio.exists(pvalfile)) {
					System.out.println("Could not find file: " + pvalfile);
				} else {
					System.out.println("in: " + pvalfile);
					ArrayList<AssociationResult> assoc = assocfile.readVariantPValues(pvalfile, region);
					ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(assoc, 0.95);
					ArrayList<Feature> credibleSetSNPList = new ArrayList<>();

					if (dsname.equals("RA") && region.getChromosome().equals(Chromosome.TWO) && region.getStart() == 60462160) {
						System.out.println("gotit");
					}

					int minPosCredSet = Integer.MAX_VALUE;
					int maxPosCredSet = 0;
					for (AssociationResult r : credibleSet) {
						Feature snp = r.getSnp();
						if (snp.getStart() < minPosCredSet) {
							minPosCredSet = snp.getStart();

						}
						if (snp.getStart() > maxPosCredSet) {
							maxPosCredSet = snp.getStart();
						}
						credibleSetSNPList.add(snp);
					}

					if (maxPosCredSet > region.getStop()) {
						maxPosCredSet = region.getStop();
					}
					if (minPosCredSet < region.getStart()) {
						minPosCredSet = region.getStart();
					}

					HashSet<AssociationResult> hashCredible = new HashSet<>();
					hashCredible.addAll(credibleSet);
					boolean[] mark = new boolean[assoc.size()];
					for (int i = 0; i < assoc.size(); i++) {
						AssociationResult r = assoc.get(i);
						if (hashCredible.contains(r)) {
							mark[i] = true;
						}
					}

					String outfile = outdir + region.toString() + ".pdf";
					try {
						System.out.println(outfile);
						annotationPlot(region, annotations, assoc, mark, annot, zoomregionHash, outfile);
					} catch (DocumentException e) {
						e.printStackTrace();
					}

					// get overlapping sequenced regions
//				ArrayList<Feature> overlappingSequencedRegions = getOverlappingRegions(region, zoomregions);
//				int max = 0;
//				int min = Integer.MAX_VALUE;
//				for (Feature queryVariantFeature : overlappingSequencedRegions) {
//					if (queryVariantFeature.getStart() < min) {
//						min = queryVariantFeature.getStart();
//					}
//					if (queryVariantFeature.getStop() > max) {
//						max = queryVariantFeature.getStop();
//					}
//				}


					int sta = minPosCredSet - 5000;
					int sto = maxPosCredSet + 5000;
					if (sta < 0) {
						sta = 0;
					}
					if (sto > Integer.MAX_VALUE) {
						sto = Integer.MAX_VALUE;
					}

					Feature newRegion = new Feature(region.getChromosome(), sta, sto);
					Pair<boolean[], ArrayList<AssociationResult>> filtered = filterAssoc(newRegion, assoc, mark);

					String zoomout = outdir + "zoom/";
					Gpio.createDir(zoomout);
					outfile = zoomout + region.toString() + ".pdf";
					try {
						System.out.println(outfile);
						annotationPlot(newRegion, annotations, filtered.getRight(), filtered.getLeft(), annot, zoomregionHash, outfile);
					} catch (DocumentException e) {
						e.printStackTrace();
					}


				}


			}
		}


	}


	private Pair<boolean[], ArrayList<AssociationResult>> filterAssoc(Feature newRegion, ArrayList<AssociationResult> assoc, boolean[] mark) {

		ArrayList<Boolean> newmark = new ArrayList<>();
		ArrayList<AssociationResult> newAssoc = new ArrayList<>();
		for (int i = 0; i < assoc.size(); i++) {
			AssociationResult r = assoc.get(i);
			if (r.getSnp().overlaps(newRegion)) {
				newmark.add(mark[i]);
				newAssoc.add(r);
			}
		}
		boolean[] outmark = new boolean[newAssoc.size()];
		for (int i = 0; i < newAssoc.size(); i++) {
			outmark[i] = newmark.get(i);
		}
		return new Pair<boolean[], ArrayList<AssociationResult>>(outmark, newAssoc);

	}

	private ArrayList<Feature> getOverlappingRegions(Feature region, ArrayList<Feature> zoomregions) {
		ArrayList<Feature> output = new ArrayList<>();
		for (Feature f : zoomregions) {
			if (f.overlaps(region)) {
				output.add(f);
			}
		}
		return output;
	}

	private void annotationPlot(Feature region, ArrayList<Track> annotations,
	                            ArrayList<AssociationResult> assoc,
	                            boolean[] mark,
	                            GTFAnnotation annot,
	                            HashSet<Feature> zoomregionhash,

	                            String outfile) throws IOException, DocumentException {

		Grid grid = new Grid(400, 300, 3, 2, 300, 300);

		// convert to arraylist of pairs
		ArrayList<Pair<Integer, Double>> pvals = convertToPairs(assoc, true);
		ArrayList<SNPFeature> snpFeatures = convertToFeatures(assoc);
		double maxp = getMax(pvals);
		ArrayList<Pair<Integer, Double>> posteriors = convertToPairs(assoc, false);
		double maxposterior = getMax(posteriors);


		TreeSet<Gene> genes = annot.getGeneTree();
		Gene geneStart = new Gene("", region.getChromosome(), Strand.POS, region.getStart(), region.getStart());
		Gene geneStop = new Gene("", region.getChromosome(), Strand.POS, region.getStop(), region.getStop());
		SortedSet<Gene> overlappingGenes = genes.subSet(geneStart, true, geneStop, true);

		ArrayList<Gene> overlappingGenesList = new ArrayList<>();
		overlappingGenesList.addAll(overlappingGenes);

		// get the subregions that overlap with this region
		// make a subplot for the zoomregion

		// plotBinaryTrait gene annotations
		GenePanel genePanel = new GenePanel(1, 1);
		genePanel.setData(region, overlappingGenesList);
		grid.addPanel(genePanel, 0, 0);
		grid.addPanel(genePanel, 0, 1);

		// plotBinaryTrait pvalues and posteriors
		AssociationPanel pvalpanel = new AssociationPanel(1, 1);
		pvalpanel.setTitle("Pvalues");
		pvalpanel.setDataSingleDs(region, zoomregionhash, pvals, "Pvalues");
		pvalpanel.setMaxPVal(maxp);
		pvalpanel.setMarkDifferentColor(mark);
		grid.addPanel(pvalpanel, 1, 0);

		AssociationPanel abfpanel = new AssociationPanel(1, 1);
		abfpanel.setTitle("Posteriors");
		abfpanel.setDataSingleDs(region, zoomregionhash, posteriors, "Posteriors");
		abfpanel.setMaxPVal(maxposterior);
		abfpanel.setMarkDifferentColor(mark);
		grid.addPanel(abfpanel, 1, 1);

		// add annotations
		AnnotationPanel annotPanel = new AnnotationPanel(1, 1);
		annotPanel.setData(region, annotations);
		annotPanel.setOverlappingFeatures(snpFeatures);
		grid.addPanel(annotPanel, 2, 0);
		grid.addPanel(annotPanel, 2, 1);

		grid.draw(outfile);
	}

	private ArrayList<SNPFeature> convertToFeatures(ArrayList<AssociationResult> assoc) {
		ArrayList<SNPFeature> output = new ArrayList<SNPFeature>();
		for (int i = 0; i < assoc.size(); i++) {
			SNPFeature f = new SNPFeature();
			Feature f2 = assoc.get(i).getSnp();
			f.setChromosome(f2.getChromosome());
			f.setName(f2.getName());
			f.setStart(f2.getStart());
			f.setStop(f2.getStop());
			f.setP(assoc.get(i).getPosterior());
			output.add(f);
		}
		return output;
	}

	private double getMax(ArrayList<Pair<Integer, Double>> pairs) {
		double max = -Double.MAX_VALUE;
		if (pairs == null) {
			System.out.println("meh.");

			System.exit(-1);
		}
		for (Pair<Integer, Double> p : pairs) {
			if (p.getRight() > max) {
				max = p.getRight();
			}
		}
		return max;
	}


	private ArrayList<Track> loadAnnotations(String listOfAnnotationsFile,
	                                         boolean usepeakcenter,
	                                         int bptoextendannotation,
	                                         boolean mergeverlapping) throws IOException {

		TextFile tf = new TextFile(listOfAnnotationsFile, TextFile.R);
		String[] list = tf.readAsArray();
		tf.close();

		ArrayList<Track> output = new ArrayList<>();
		AnnotationLoader loader = new AnnotationLoader();
		for (int i = 0; i < list.length; i++) {
			Track t = loader.loadAnnotations(list[i], usepeakcenter, bptoextendannotation, mergeverlapping, null);
			output.add(t);
		}
		return output;
	}


	public void makeplotsOct27() throws IOException {

		String sequencedRegionsFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/2015-03-30-allRegionsWith50PercentOfSamplesAbove20x.bed";
		String assocDir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/";
		String regionsFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-05-18-AllRegions/allLoci.bed";

		String rsquareDir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/ImpQScores/";
		String[] refs = new String[]{"1kg", "seq", "1kg-seq-merged"};
		String gtf = "/Data/Ref/Annotation/UCSC/genes.gtf";
		String outdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/Plots/";

		String condiontalassocDir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/";
		String conditionaloutdir = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-10-25-Assoc/Conditional/Plots/";

		String immunoChipT1D = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/T1D-Onengut/hg19_gwas_ic_t1d_onengut_cc_4_18_1.tab";
		String immunoChipRA = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-06-04-Assoc/RA-Eyre/hg19_gwas_ic_ra_eyre_4_18_0.tab";


		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(regionsFile);
		ArrayList<Feature> sequencedRegions = reader.readAsList(sequencedRegionsFile);
		HashSet<Feature> sequencedRegionsHash = new HashSet<Feature>();
		sequencedRegionsHash.addAll(sequencedRegions);

		double bayesthreshold = 0.95;

		GTFAnnotation annot = new GTFAnnotation(gtf);


		for (int r = 0; r < regions.size(); r++) {
			Feature region = regions.get(r);
			TreeSet<Gene> genes = annot.getGeneTree();
			Gene geneStart = new Gene("", region.getChromosome(), Strand.POS, region.getStart(), region.getStart());
			Gene geneStop = new Gene("", region.getChromosome(), Strand.POS, region.getStop(), region.getStop());
			SortedSet<Gene> overlappingGenes = genes.subSet(geneStart, true, geneStop, true);

			ArrayList<Gene> overlappingGenesList = new ArrayList<>();
			overlappingGenesList.addAll(overlappingGenes);

			// make conditional plots
			int nrRuns = 5;
			for (int d = 0; d < 2; d++) {
				String ds = "T1D";
				String ic = immunoChipT1D;
				if (d == 1) {
					ds = "RA";
					ic = immunoChipRA;
				}


				GenePanel genePanel = new GenePanel(1, 1);
				genePanel.setData(region, overlappingGenesList);

				// genes / pvals / rsquared / maf
				Grid grid = new Grid(400, 300, 4, refs.length + 1, 100, 100);
				grid.addPanel(genePanel, 0, refs.length);

				// load IC p-values

				AssociationFile f = new AssociationFile();
				ArrayList<AssociationResult> icData = f.readVariantPValues(ic, region);

				ArrayList<Pair<Integer, Double>> icPvals = convertToPairs(icData, true);

				AssociationPanel icPanel = new AssociationPanel(1, 1);
				icPanel.setDataSingleDs(region, sequencedRegionsHash, icPvals, "Association stats ImmunoChip study");


				Grid pvalgrid = new Grid(400, 300, nrRuns + 1, refs.length, 100, 100);
				Grid bfgrid = new Grid(400, 300, nrRuns + 1, refs.length, 100, 100);
				Grid abfgrid = new Grid(400, 300, nrRuns + 1, refs.length, 100, 100);


				ArrayList<AssociationPanel> assocPanels = new ArrayList<>();
				ArrayList<AssociationPanel> assocPanelsBayes = new ArrayList<>();
				ArrayList<AssociationPanel> assocPanelsABF = new ArrayList<>();

				double maxICP = -Double.MAX_VALUE;
				for (Pair<Integer, Double> p : icPvals) {
					if (p.getRight() > maxICP) {
						maxICP = p.getRight();
					}
				}

				double[] maxBayesPerRun = new double[nrRuns];
				double[] maxPosteriorPerRun = new double[nrRuns];
				double[] maxPvalPerRun = new double[nrRuns];
				maxPvalPerRun[0] = maxICP;

				ApproximateBayesPosterior abp = new ApproximateBayesPosterior();

				for (int refId = 0; refId < refs.length; refId++) {
					// read gwas file
					String ref = refs[refId];

					Chromosome chr = region.getChromosome();

					String regionStr = region.toString();


					grid.addPanel(genePanel, 0, refId);
					pvalgrid.addPanel(genePanel, 0, refId);
					abfgrid.addPanel(genePanel, 0, refId);
					bfgrid.addPanel(genePanel, 0, refId);


					for (int run = 0; run < nrRuns; run++) {
						double maxP = -Double.MAX_VALUE;
						double maxBayes = -Double.MAX_VALUE;
						double maxAverageBayes = -Double.MAX_VALUE;
						String assocFile = condiontalassocDir + "/Conditional/" + ds + "/" + ref + "/" + chr.toString() + "-" + regionStr + "-gwas-" + run + ".txt";

						if (Gpio.exists(assocFile)) {
//							String ldFile = condiontalassocDir + "/Conditional/" + ds + "/" + ref + "/" + chr.toString() + "-" + regionStr + "-ld-" + run + ".txt";
//							ArrayList<Triple<Integer, Double, Double>> ld = loadLDInfo(ldFile);
//							ArrayList<Pair<Integer, Double>> ldsqr = new ArrayList<Pair<Integer, Double>>();
//							for(Triple<Integer, Double, Double> l: ld){
//								ldsqr.add(new Pair<Integer, Double>(l.getLeft(),l.getMiddle()));
//							}

							AssociationFile associationFile = new AssociationFile();
							ArrayList<AssociationResult> associationResults = associationFile.read(assocFile, region);

							String model = associationFile.getModel();

							ArrayList<Pair<Integer, Double>> pvals = new ArrayList<Pair<Integer, Double>>();
							for (int q = 0; q < associationResults.size(); q++) {
								AssociationResult p = associationResults.get(q);
								Pair<Integer, Double> pair = new Pair<Integer, Double>(p.getSnp().getStart(), p.getPval());
								double pval = p.getPval();
								if (pval > maxP) {
									maxP = pval;
								}
								p.setBf(0);
								p.setPosterior(0);
								pvals.add(pair);
							}

							// calculate abf
							abp.calculatePosterior(associationResults);

							// make credible sets
							ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(associationResults, bayesthreshold);
							HashSet<AssociationResult> credibleSetSet = new HashSet<AssociationResult>();


							credibleSetSet.addAll(credibleSet);
							ArrayList<Pair<Integer, Double>> averageBayesFactors = new ArrayList<Pair<Integer, Double>>();
							boolean[] mark = new boolean[associationResults.size()];

							ArrayList<Pair<Integer, Double>> bayesFactors = new ArrayList<Pair<Integer, Double>>();
							for (int q = 0; q < associationResults.size(); q++) {
								AssociationResult p = associationResults.get(q);
								Pair<Integer, Double> data = new Pair<Integer, Double>(p.getSnp().getStart(), p.getBf());
								double bf = p.getBf();
								if (bf > maxBayes) {
									maxBayes = bf;
								}
								bayesFactors.add(data);

								double abf = p.getPosterior();
								if (run > 0 && abf > 0.9) {
									System.out.println("ERROR");
								}
								if (!Double.isNaN(abf)) {
									Pair<Integer, Double> abfdata = new Pair<Integer, Double>(p.getSnp().getStart(), abf);
									if (abf > maxAverageBayes) {
										maxAverageBayes = abf;
									}
									if (credibleSetSet.contains(p)) {
										mark[q] = true;
									}
									averageBayesFactors.add(abfdata);
								}

							}


							String modelStr = "";
							if (model != null) {
								modelStr = model;
								modelStr = "Conditional " + (run);//modelStr.replaceAll("#model: status ~ snp + covar +", "");
							}

							AssociationPanel pvalpanel = new AssociationPanel(1, 1);
							pvalpanel.setTitle("Association stats " + ref + " - " + modelStr);
//							pvalpanel.setLdInfo(ldsqr);
							pvalpanel.setDataSingleDs(region, sequencedRegionsHash, pvals, "Association stats " + ref + " - " + modelStr);
							pvalgrid.addPanel(pvalpanel, run + 1, refId);


							// plotBinaryTrait rsquared and maf
							if (run == 0) {
								grid.addPanel(pvalpanel, 1, refId);

								String rsquareFile = rsquareDir + ds + "/" + ref + "/merged-variantinfo-Chr" + chr.getNumber() + ".txt";
								VCFRSquares q = new VCFRSquares();
								ArrayList<Pair<Feature, Double>> rsquareds = q.loadRSquareds(rsquareFile, region);
								RSquaredPanel p3 = new RSquaredPanel(1, 1);

								ArrayList<Feature> inputBeforeImputation = null;
								p3.setData(region, sequencedRegionsHash, "Imputation qual scores " + ref, inputBeforeImputation, rsquareds);
								grid.addPanel(p3, 2, refId);

								ArrayList<Pair<Integer, Double>> mafs = new ArrayList<Pair<Integer, Double>>();
								for (AssociationResult p : associationResults) {
									Pair<Integer, Double> maf = new Pair<Integer, Double>(p.getSnp().getStart(), p.getMaf());
									mafs.add(maf);
								}

								AssociationPanel mafPanel = new AssociationPanel(1, 1);
								mafPanel.setDataSingleDs(region, sequencedRegionsHash, mafs, "Minor Allele Frequency " + ref);
								mafPanel.setPlotGWASSignificance(false);
								grid.addPanel(mafPanel, 3, refId);

							}

							AssociationPanel bayespanel = new AssociationPanel(1, 1);
							bayespanel.setDataSingleDs(region, sequencedRegionsHash, bayesFactors, "Bayes factors " + ref + " - " + modelStr);
							bayespanel.setPlotGWASSignificance(false);
//							bayespanel.setLdInfo(ldsqr);
							bfgrid.addPanel(bayespanel, run + 1, refId);

							AssociationPanel averagebayespanel = new AssociationPanel(1, 1);
							averagebayespanel.setPlotGWASSignificance(false);
//							averagebayespanel.setLdInfo(ldsqr);
							averagebayespanel.setMarkDifferentColor(mark);
							averagebayespanel.setDataSingleDs(region, sequencedRegionsHash, averageBayesFactors, "Posteriors " + ref + " - " + modelStr);
							abfgrid.addPanel(averagebayespanel, run + 1, refId);

							pvalpanel.setMarkDifferentColor(mark);
							assocPanels.add(pvalpanel);
							assocPanelsBayes.add(bayespanel);
							assocPanelsABF.add(averagebayespanel);

							if (maxAverageBayes > maxPosteriorPerRun[run]) {
								maxPosteriorPerRun[run] = maxAverageBayes;
							}
							if (maxBayes > maxBayesPerRun[run]) {
								maxBayesPerRun[run] = maxBayes;
							}
							if (maxP > maxPvalPerRun[run]) {
								maxPvalPerRun[run] = maxP;
							}


						} else {
							System.out.println("Could not find file: " + assocFile);
							// add spacer panels...

							pvalgrid.addPanel(new SpacerPanel(1, 1), run + 1, refId);
							bfgrid.addPanel(new SpacerPanel(1, 1), run + 1, refId);
							abfgrid.addPanel(new SpacerPanel(1, 1), run + 1, refId);
						}
					}
				}

				for (int refId = 0; refId < refs.length; refId++) {
					for (int run = 0; run < nrRuns; run++) {

						Panel pvalPanel = pvalgrid.getPanel(run + 1, refId);
						if (pvalPanel instanceof AssociationPanel) {
							AssociationPanel p = (AssociationPanel) pvalPanel;
							p.setMaxPVal(maxPvalPerRun[run]);
						}
						Panel bfPanel = bfgrid.getPanel(run + 1, refId);
						if (bfPanel instanceof AssociationPanel) {
							AssociationPanel p = (AssociationPanel) bfPanel;
							p.setMaxPVal(maxBayesPerRun[run]);
						}

						Panel abfPanel = abfgrid.getPanel(run + 1, refId);
						if (abfPanel instanceof AssociationPanel) {
							AssociationPanel p = (AssociationPanel) abfPanel;
							p.setMaxPVal(maxPosteriorPerRun[run]);
						}
					}
				}


				icPanel.setTitle("IcPanel");
				icPanel.setMaxPVal(maxPvalPerRun[0]);


				grid.addPanel(icPanel, 1, refs.length);

				try {
					Gpio.createDir(conditionaloutdir);
					String outfile = conditionaloutdir + ds + "-" + region.toString() + "-pvals-rsquares-mafs.pdf";
					grid.draw(outfile);
					outfile = conditionaloutdir + ds + "-" + region.toString() + "-pvals.pdf";
					pvalgrid.draw(outfile);
					outfile = conditionaloutdir + ds + "-" + region.toString() + "-bf.pdf";
					bfgrid.draw(outfile);
					outfile = conditionaloutdir + ds + "-" + region.toString() + "-posterior.pdf";
					abfgrid.draw(outfile);
				} catch (DocumentException e) {
					e.printStackTrace();
				}
			}

		}

	}

	private ArrayList<Pair<Integer, Double>> convertToPairs(ArrayList<AssociationResult> icData, boolean pval) {
		ArrayList<Pair<Integer, Double>> output = new ArrayList<>();
		if (pval) {
			for (AssociationResult s : icData) {
				output.add(new Pair<Integer, Double>(s.getSnp().getStart(), s.getPval()));
			}

		} else {
			for (AssociationResult s : icData) {
				output.add(new Pair<Integer, Double>(s.getSnp().getStart(), s.getPosterior()));
			}
		}
		return output;
	}

	public ArrayList<Triple<Integer, Double, Double>> loadLDInfo(String file) throws IOException {
		TextFile tf = new TextFile(file, TextFile.R);

		ArrayList<Triple<Integer, Double, Double>> output = new ArrayList<>();
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length > 3) {
				String snp = elems[1];
				String rsq = elems[2];
				String dpr = elems[3];

				Double rsqd = Double.parseDouble(rsq);
				Double dprd = Double.parseDouble(dpr);

				String[] snpelems = snp.split("-");
				String[] posElems = snpelems[0].split(":");
				Integer pos = Integer.parseInt(posElems[1]);


				Triple<Integer, Double, Double> t = new Triple<Integer, Double, Double>(pos, rsqd, dprd);
				output.add(t);
			}
			elems = tf.readLineElems(TextFile.tab);
		}

		tf.close();
		return output;

	}
}
