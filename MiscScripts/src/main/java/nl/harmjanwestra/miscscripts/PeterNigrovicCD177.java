package nl.harmjanwestra.miscscripts;

import com.itextpdf.text.DocumentException;
import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Strand;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.AssociationPanel;
import nl.harmjanwestra.utilities.graphics.panels.GenePanel;
import nl.harmjanwestra.utilities.graphics.panels.LDPanel;
import nl.harmjanwestra.utilities.annotation.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.*;

/**
 * Created by hwestra on 9/29/16.
 */
public class PeterNigrovicCD177 {


	public static void main(String[] args) {
		String tabixrefprefix = "/Data/Ref/beagle_1kg/1kg.phase3.v5a.chr";
		String ib = "/Data/ImmunoChip/ImmunoBase/ImmunoChip/hg19_gwas_ra_okada_4_19_1.tab";
//		String ib = "/Data/ImmunoChip/ImmunoBase/ImmunoChip/hg19_gwas_ic_ra_eyre_4_19_1.tab";

		String bed = "/Data/Projects/PeterNigrovic/CD177/cd177.bed";
		String annotationfile = "/Data/Ref/Annotation/UCSC/genes.gtf";
		String tabixsamplelimit = "/Data/Ref/1kg-europeanpopulations.txt.gz";
		String out = "/Data/Projects/PeterNigrovic/CD177/cd177out";
		String eqtlfile = "/Data/Projects/PeterNigrovic/CD177/cd177bioseqtl.txt";

		PeterNigrovicCD177 c = new PeterNigrovicCD177();
		try {
			c.runImmunoBase(ib, bed, annotationfile, tabixrefprefix, tabixsamplelimit, out, "rs201821720", 0.1, eqtlfile);

		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}


	}

	public void runImmunoBase(String ib,
							  String bed,
							  String annotationfile,
							  String tabixrefprefix,
							  String tabixsamplelimit,
							  String out,
							  String queryVariant,
							  double ldthreshold,
							  String eqtlfile) throws IOException, DocumentException {

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bed);

		GTFAnnotation annotation = new GTFAnnotation(annotationfile);
		for (Feature f : regions) {
			AssociationFile associationFile = new AssociationFile();
			ArrayList<AssociationResult> associationresults = associationFile.readRegion(ib, f);

			ArrayList<Pair<Integer, Double>> allPvalue = new ArrayList<Pair<Integer, Double>>();
			TextFile outtf = new TextFile(out + f.toString() + "-assoc.txt", TextFile.W);
			outtf.writeln(associationFile.getHeader());

//			HashSet<String> allSNPs = new HashSet<String>();
			for (AssociationResult r : associationresults) {
				outtf.writeln(r.toString());
//				allSNPs.add(r.getSnp().getName());
				allPvalue.add(new Pair<>(r.getSnp().getStart(), r.getLog10Pval()));
			}
			outtf.close();


			// plot the region
			TreeSet<Gene> genes = annotation.getGeneTree();
			Gene geneStart = new Gene("", f.getChromosome(), Strand.POS, f.getStart(), f.getStart());
			Gene geneStop = new Gene("", f.getChromosome(), Strand.POS, f.getStop(), f.getStop());
			SortedSet<Gene> overlappingGenes = genes.subSet(geneStart, true, geneStop, true);

			ArrayList<Gene> genesList = new ArrayList<>();
			genesList.addAll(overlappingGenes);


			Grid grid = new Grid(400, 400, 2, 1, 100, 100);
			GenePanel genepanel = new GenePanel(1, 1);
			genepanel.setData(f, genesList);
			grid.addPanel(genepanel);


			AssociationPanel assocPanel = new AssociationPanel(1, 1);
			assocPanel.setDataSingleDs(f, null, allPvalue, "");

			grid.addPanel(assocPanel);
			grid.draw(out + f.toString() + ".pdf");


			// plot the LD
			grid = new Grid(400, 300, 2, 1, 100, 100);


			String tabixfile = tabixrefprefix + f.getChromosome().getNumber() + ".vcf.gz";
			boolean[] samplesToInclude = null;

			if (tabixsamplelimit != null) {
				VCFGenotypeData d = new VCFGenotypeData(tabixfile);
				ArrayList<String> tabixSamples = d.getSamples();
				samplesToInclude = new boolean[tabixSamples.size()];

				TextFile tf = new TextFile(tabixsamplelimit, TextFile.R);
				String[] elems = tf.readLineElems(Strings.whitespace);
				HashSet<String> allSamplesInListSet = new HashSet<String>();
				while (elems != null) {
					allSamplesInListSet.add(elems[0]);
					elems = tf.readLineElems(Strings.whitespace);
				}
				tf.close();

				for (int i = 0; i < tabixSamples.size(); i++) {
					if (allSamplesInListSet.contains(tabixSamples.get(i))) {
						samplesToInclude[i] = true;
					}
				}
			}

			TabixReader treader = new TabixReader(tabixfile);
			TabixReader.Iterator window = treader.query(f.getChromosome().getNumber() + ":" + (f.getStart() - 10) + "-" + (f.getStop() + 10));

			ArrayList<VCFVariant> vcfVariants = new ArrayList<VCFVariant>();
			String next = window.next();
			TextFile allvarsout = new TextFile(out + f.toString() + "-all1kgvars.txt", TextFile.W);
			allvarsout.writeln("chr\tpos\tid\talleles\tallelefrequencies\tminorAllele\tMAF\thwep");
			while (next != null) {
				VCFVariant variant = new VCFVariant(next, VCFVariant.PARSE.ALL, samplesToInclude);
				vcfVariants.add(variant);
				variant.recalculateMAFAndCallRate();
				variant.calculateHWEP();
				allvarsout.writeln(variant.getChr().toString()
						+ "\t" + variant.getPos()
						+ "\t" + variant.getId()
						+ "\t" + Strings.concat(variant.getAlleles(), Strings.semicolon)
						+ "\t" + Strings.concat(variant.getAlleleFrequencies(), Strings.semicolon)
						+ "\t" + variant.getMinorAllele()
						+ "\t" + variant.getMAF()
						+ "\t" + variant.getHwep()
				);
				next = window.next();
			}
			allvarsout.close();

			ArrayList<Pair<Integer, Integer>> positions = new ArrayList<>();
			ArrayList<Double> ld = new ArrayList<>();

			TextFile allvarsoutld = new TextFile(out + f.toString() + "-all1kgvars-ld.txt", TextFile.W);
			allvarsoutld.writeln("Pos1\tId1\tPos2\tId2\tR-squared");
			DetermineLD ldcalc = new DetermineLD();

			VCFVariant querySNP = null;
			HashMap<String, Double> proxyvariants = new HashMap<String, Double>();
			for (int i = 0; i < vcfVariants.size(); i++) {
				VCFVariant var1 = vcfVariants.get(i);
				for (int j = i + 1; j < vcfVariants.size(); j++) {
					VCFVariant var2 = vcfVariants.get(j);
					Pair<Double, Double> d = ldcalc.getLD(var1, var2);
					if (queryVariant != null && (queryVariant.equals(var1.getId()) || queryVariant.equals(var2.getId()))) {
						double r = d.getRight();

						if (queryVariant.equals(var1.getId())) {
							querySNP = var1;
							proxyvariants.put(var2.getId(), r);
						} else {
							querySNP = var2;
							proxyvariants.put(var1.getId(), r);
						}
					}


					if (d != null && !Double.isNaN(d.getRight())) {
						positions.add(new Pair<>(var1.getPos(), var2.getPos()));
						ld.add(d.getRight());
						allvarsoutld.writeln(var1.getPos() + "\t" + var1.getId() + "\t" + var2.getPos() + "\t" + var2.getId() + "\t" + d.getRight());
					} else {
						allvarsoutld.writeln(var1.getPos() + "\t" + var1.getId() + "\t" + var2.getPos() + "\t" + var2.getId() + "\tNaN");
//						System.out.println("edgecase:\t" + var1.getPos()
//								+ "\t" + Strings.concat(var1.getAlleles(), Strings.semicolon)
//								+ "\t" + var1.getId()
//								+ "\t" + var2.getPos()
//								+ "\t" + Strings.concat(var2.getAlleles(), Strings.semicolon)
//								+ "\t" + var2.getId()
//								+ "\tNaN");
					}
				}
			}
			allvarsoutld.close();

			LDPanel ldPanel = new LDPanel(1, 1);
			ldPanel.setData(f, positions, ld);
			grid.addPanel(ldPanel);
			grid.draw(out + f.toString() + "-ld.pdf");

			int alleleAA = 0;
			int alleleAB = 0;
			int alleleBB = 0;
			byte[] gtvector = querySNP.getGenotypesAsByteVector();
			for (int i = 0; i < querySNP.getNrSamples(); i++) {
				if (gtvector[i] == 0) {
					alleleAA++;
				} else if (gtvector[i] == 1) {
					alleleAB++;
				} else {
					alleleBB++;
				}
			}


			outtf = new TextFile(out + f.toString() + "-assoc-proxyfinder.txt", TextFile.W);
			outtf.writeln(associationFile.getHeader() + "\tLD");


//			HashSet<String> allSNPs = new HashSet<String>();
			for (AssociationResult r : associationresults) {
				if (proxyvariants.containsKey(r.getSnp().getName())) {
					outtf.writeln(r.toString() + "\t" + proxyvariants.get(r.getSnp().getName()));
				}
			}

			outtf.close();

			System.out.println();
			System.out.println("Allele counts");
			System.out.println(querySNP.getAlleles()[0] + querySNP.getAlleles()[0] + "\t" + alleleAA);
			System.out.println(querySNP.getAlleles()[0] + querySNP.getAlleles()[1] + "\t" + alleleAB);
			System.out.println(querySNP.getAlleles()[1] + querySNP.getAlleles()[1] + "\t" + alleleBB);

			TextFile tfein = new TextFile(eqtlfile, TextFile.R);
			TextFile tfeout = new TextFile(out + f.toString() + "-assoc-eqtls.txt", TextFile.W);
			tfeout.writeln("Pval\tSNP\tLD");
			String[] eelems = tfein.readLineElems(TextFile.tab);
			while (eelems != null) {
				String pval = eelems[0];
				String snp = eelems[1];
				if (proxyvariants.containsKey(snp)) {
					tfeout.writeln(pval + "\t" + snp + "\t" + proxyvariants.get(snp));
				}
				eelems = tfein.readLineElems(TextFile.tab);
			}
			tfeout.close();

			tfein.close();
		}


	}

}
