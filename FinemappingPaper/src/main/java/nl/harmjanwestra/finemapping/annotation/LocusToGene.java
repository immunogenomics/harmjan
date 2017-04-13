package nl.harmjanwestra.finemapping.annotation;

import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Strand;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.annotation.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
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
//			t.determineRegionSignificanceThresholds(annot, bed, output);

			output = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/AllLoci-GenesPerLocus.txt";
			String file = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/META-assoc0.3-COSMO-merged-posterior.txt.gz-credibleSets-sets.txt";
			String fileout = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/META-assoc0.3-COSMO-merged-posterior.txt.gz-credibleSets-sets-wGenes.txt";
			t.runset(output, file, fileout);

			file = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/RA-assoc0.3-COSMO-merged-posterior.txt.gz-credibleSets-sets.txt";
			fileout = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/RA-assoc0.3-COSMO-merged-posterior.txt.gz-credibleSets-sets-wGenes.txt";
			t.runset(output, file, fileout);

			file = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/T1D-assoc0.3-COSMO-merged-posterior.txt.gz-credibleSets-sets.txt";
			fileout = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-07-25-SummaryStats/Normal/T1D-assoc0.3-COSMO-merged-posterior.txt.gz-credibleSets-sets-wGenes.txt";
			t.runset(output, file, fileout);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void runset(String output, String file, String fileout) throws IOException {

		HashMap<String, String> locusToGene = new HashMap<String, String>();
		TextFile tf = new TextFile(output, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length > 1) {
				locusToGene.put(elems[0], elems[1]);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		TextFile out = new TextFile(fileout, TextFile.W);
		TextFile in = new TextFile(file, TextFile.R);
		out.writeln(in.readLine() + "\tGenes");
		elems = in.readLineElems(Strings.whitespace);
		while (elems != null) {
			if (elems.length > 1) {
				String loc = elems[0];
				String genes = locusToGene.get(loc);
				out.writeln(loc + "\t" + elems[1] + "\t" + genes);
			}
			elems = in.readLineElems(Strings.whitespace);
		}
		out.close();
		in.close();
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
				genesInlocus.add(g.getGeneSymbol());
			}

			out.writeln(region.toString() + "\t" + Strings.concat(genesInlocus, Strings.semicolon));

		}
		out.close();

	}
}
