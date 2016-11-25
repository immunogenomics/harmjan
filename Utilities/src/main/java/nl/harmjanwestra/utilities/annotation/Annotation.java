package nl.harmjanwestra.utilities.annotation;

import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.features.Transcript;

import java.util.Collection;
import java.util.HashMap;
import java.util.TreeSet;

/**
 * Created by hwestra on 11/14/16.
 */
public class Annotation {

	protected Collection<Gene> genes;
	protected HashMap<String, Gene> strToGene;
	protected String annotationLocation;


	public HashMap<String, Gene> getStrToGene() {
		return strToGene;
	}

	public TreeSet<Gene> getGeneTree() {
		FeatureComparator comp = new FeatureComparator(false);
		TreeSet<Gene> geneTree = new TreeSet<Gene>(comp);
		for (Gene g : genes) {
			g.getBounds();
			geneTree.add(g);
		}
		comp.setAllowOverlap(true);
		return geneTree;
	}

	public Collection<Gene> getGenes() {
		return genes;
	}

	public TreeSet<Transcript> getTranscriptTree() {
		TreeSet<Transcript> transcriptTree = new TreeSet<Transcript>(new FeatureComparator(true));
		for (Gene g : genes) {
			transcriptTree.addAll(g.getTranscripts());
		}
		return transcriptTree;
	}

}
