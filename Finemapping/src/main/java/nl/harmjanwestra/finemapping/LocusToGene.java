package nl.harmjanwestra.finemapping;

import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.features.Strand;
import nl.harmjanwestra.utilities.gtf.GTFAnnotation;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Created by Harm-Jan on 06/19/16.
 */
public class LocusToGene {

	public static void main(String[] args) {

		LocusToGene t = new LocusToGene();
		String annot = "d:/Data/Annotation/UCSC/genes.gtf";
		String bed = "D:\\Cloud\\Dropbox\\2016-03-RAT1D-Finemappng\\Data\\AllICLoci.bed";
		String output = "D:\\Cloud\\Dropbox\\2016-03-RAT1D-Finemappng\\Data\\AllLoci-GenesPerLocus.txt";
		try {
			t.run(annot, bed, output);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String annot, String bed, String output) throws IOException {

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> loci = reader.readAsList(bed);

		GTFAnnotation annotation = new GTFAnnotation(annot);
		TreeSet<Gene> genes = annotation.getGeneTree();

		TextFile out = new TextFile(output, TextFile.W);
		for (Feature region : loci) {
			Gene geneStart = new Gene("", region.getChromosome(), Strand.POS, region.getStart(), region.getStart());
			Gene geneStop = new Gene("", region.getChromosome(), Strand.POS, region.getStop(), region.getStop());
			SortedSet<Gene> overlappingGenes = genes.subSet(geneStart, true, geneStop, true);

			ArrayList<String> genesInlocus = new ArrayList<>();

			for (Gene g : overlappingGenes) {
				genesInlocus.add(g.getGeneId());
			}

			out.writeln(region.toString() + "\t" + Strings.concat(genesInlocus, Strings.semicolon));

		}
		out.close();

	}
}
