package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.FeatureMerger;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.annotation.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.*;

/**
 * Created by hwestra on 5/13/15.
 */
public class BedFileFromGTF {

	public static void main(String[] args) {

		String[] genes = new String[]{"CD8", "CD4", "CD3", "CD28", "IFNG", "IL2", "TNFA", "STAT1", "VWCE"};
		String gtfFile = "/Data/Annotation/UCSC/genes.gtf";
		String outfile = "/Data/tmp/genes.bed";
		int margin = 1000;
		BedFileFromGTF f = new BedFileFromGTF();
		try {
			f.createBed(genes, gtfFile, outfile, margin);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void createBed(String[] genes, String gtfFile, String outfile, int margin) throws IOException {

		GTFAnnotation annot = new GTFAnnotation(gtfFile);
		Collection<Gene> genesList = annot.getGenes();

		HashSet<String> queryGenes = new HashSet<String>();
		queryGenes.addAll(Arrays.asList(genes));

		ArrayList<Feature> features = new ArrayList<Feature>();

		for (Gene g : genesList) {
			if (queryGenes.contains(g.getGeneSymbol())) {
				String chrstr = g.getChromosome().toString();
				int start = g.getStart();
				int stop = g.getStop();

				Feature f = new Feature();
				f.setChromosome(g.getChromosome());
				f.setStart(start - margin);
				f.setStop(stop + margin);
				features.add(f);


			}
		}

		features = FeatureMerger.merge(features, true);
		Collections.sort(features, new FeatureComparator(false));

		TextFile tf = new TextFile(outfile, TextFile.W);
		for (Feature f : features) {
			System.out.println(f.toString());
			tf.writeln(f.getChromosome().toString() + "\t" + f.getStart() + "\t" + f.getStop());
		}
		tf.close();

	}
}
