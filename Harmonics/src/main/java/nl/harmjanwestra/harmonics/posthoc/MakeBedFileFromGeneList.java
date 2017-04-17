package nl.harmjanwestra.harmonics.posthoc;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import umcg.genetica.ensembl.Features;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 10/14/15.
 */
public class MakeBedFileFromGeneList {

	public static void main(String[] args) {


		MakeBedFileFromGeneList m = new MakeBedFileFromGeneList();

		String annot = "/Data/Ref/Annotation/Ensembl/structures.txt";
		String genelist = "/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/ExpressionMaria/gene_RPKMs_TMMnorm_effectiveLength_exprCMP5in1_cd4tlpAllSamples_woExtraCols.QuantileNormalized.Log2Transformed-topGenes.txt";
		String bedout = "/Data/Projects/2015-cd4timelinepilot/2015-08-31-Peaks/ExpressionMaria/gene_RPKMs_TMMnorm_effectiveLength_exprCMP5in1_cd4tlpAllSamples_woExtraCols.QuantileNormalized.Log2Transformed-topGenes.saf";
		try {
			m.run(annot, genelist, bedout);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void run(String ensemblAnnotation, String genelist, String bedout) throws IOException {
		Features feat = new Features();
		feat.loadAnnotation(ensemblAnnotation);
		HashMap<String, umcg.genetica.containers.Gene> genHash = feat.getGeneHash();

		HashSet<String> genes = new HashSet<String>();
		TextFile tf = new TextFile(genelist, TextFile.R);

		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {

			String gene = elems[0];
			genes.add(gene);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		ArrayList<String> allGenes = new ArrayList<String>();
		allGenes.addAll(genes);

		TextFile bedFileout = new TextFile(bedout, TextFile.W);

		for (int i = 0; i < allGenes.size(); i++) {
			String gene = allGenes.get(i);

			umcg.genetica.containers.Gene g = genHash.get(gene);
			int start = g.getStart();
			int stop = g.getEnd();
			int strand = g.getStrand();

			if (strand < 0) {
				System.out.println("neg strand: " + start + "\t" + stop);
				bedFileout.writeln(g.getName() + "-" + g.getAnnotation() + "\tchr" + g.getParentChromosome().getName() + "\t" + (stop) + "\t" + (stop + 1000) + "\t+");
			} else {
				bedFileout.writeln(g.getName() + "-" + g.getAnnotation() + "\tchr" + g.getParentChromosome().getName() + "\t" + (start - 1000) + "\t" + start + "\t+");
			}


		}
		bedFileout.close();


	}
}
