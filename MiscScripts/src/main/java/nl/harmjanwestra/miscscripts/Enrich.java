package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.annotation.ensembl.EnsemblStructures;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.features.Transcript;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.FisherExactTest;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

public class Enrich {
	
	
	public static void main(String[] args) {
		// meh
		
		Enrich e = new Enrich();
		String ensembl = "C:\\Data\\Ref\\Ensembl\\GrCH37-b86-Structures.txt.gz";
		String efile = "C:\\Data\\LudeFranke\\SumOfPrunedSNPEffectOnGene.txt";
		String pfile = "C:\\Data\\LudeFranke\\fordist_cleaned_nonpsych_z_pli_rec_null_data.txt";
		String outfile = "C:\\Data\\LudeFranke\\Comp.txt";
		try {
			e.run(ensembl, efile, pfile, outfile);
		} catch (IOException e1) {
			e1.printStackTrace();
		}
	}
	
	public void run(String ensembl, String efile, String pfile, String outf) throws IOException {
		
		
		EnsemblStructures s = new EnsemblStructures(ensembl);
		
		
		TextFile tf = new TextFile(efile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		
		HashMap<String, Double> geneToEQTLZ = new HashMap<String, Double>();
		HashMap<String, Double> geneToEQTLE = new HashMap<String, Double>();
		
		while (elems != null) {
			String gene = elems[0];
			String z = elems[1];
			String e = elems[2];
			
			Double zd = Double.parseDouble(z);
			Double ed = Double.parseDouble(e);
			
			geneToEQTLE.put(gene, ed);
			geneToEQTLZ.put(gene, zd);
			
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		
		
		HashMap<String, String> tTog = new HashMap<String, String>();
		Collection<Gene> genes = s.getGenes();
		for (Gene g : genes) {
			ArrayList<Transcript> transcripts = g.getTranscripts();
			for (Transcript t : transcripts) {
				String tname = t.getName();
				String gene = t.getGene().getName();
				System.out.println(tname + "\t" + gene);
				tTog.put(tname, gene);
			}
		}
		
		System.out.println("");
		
		TextFile tf2 = new TextFile(pfile, TextFile.R);
		elems = tf2.readLineElems(TextFile.tab);
		int bpcol = -1;
		int plicol = -1;
		int plizcol = -1;
		int synzcol = -1;
		int miszcol = -1;
		
		// syn_z	mis_z
		for (int i = 0; i < elems.length; i++) {
			if (elems[i].toLowerCase().equals("bp")) {
				bpcol = i;
			} else if (elems[i].toLowerCase().equals("pli")) {
				plicol = i;
			} else if (elems[i].toLowerCase().equals("lof_z")){
				plizcol = i;
			}else if (elems[i].toLowerCase().equals("syn_z")){
				synzcol = i;
			}else if (elems[i].toLowerCase().equals("mis_z")){
				miszcol = i;
			}
		}
		
		System.out.println(plicol);
		System.out.println(bpcol);
		
		elems = tf2.readLineElems(TextFile.tab);
		
		HashMap<String, Double> geneToPLI = new HashMap<String, Double>();
		HashMap<String, Double> geneToBP = new HashMap<String, Double>();
		HashMap<String, Double> geneToPLIZ = new HashMap<String, Double>();
		HashMap<String, Double> geneToMisZ = new HashMap<String, Double>();
		HashMap<String, Double> geneToSynZ = new HashMap<String, Double>();
		
		while (elems != null) {
			String transcript = elems[0];
			String[] telems = Strings.dot.split(transcript);
			transcript = telems[0];
			String gene = tTog.get(transcript);
//			System.out.println(transcript + "\t" + gene);
			
			if (gene != null) {
				if (geneToPLI.containsKey(gene)) {
					System.out.println("Duplicate transcript: " + transcript + " for gene " + gene);
				} else {
					String bp = elems[bpcol];
					String pli = elems[plicol];
					String pliz = elems[plizcol];
					String misz = elems[miszcol];
					String synz = elems[synzcol];
					double plid = Double.parseDouble(pli);
					double bpd = Double.parseDouble(bp);
					double plizd = Double.parseDouble(pliz);
					double synzd = Double.parseDouble(synz);
					double miszd = Double.parseDouble(misz);
					
					geneToPLI.put(gene, plid);
					geneToBP.put(gene, bpd);
					geneToPLIZ.put(gene, plizd);
					geneToSynZ.put(gene, synzd);
					geneToMisZ.put(gene, miszd);
				}
			}
			
			
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		
		
		// merge
		TextFile out = new TextFile(outf, TextFile.W);
		out.writeln("Gene\tTrans-Squared-Z\tGeneExp\tGeneLen\tpLI\tLIZ\tSYNZ\tMISZ");
		for (String gene : geneToEQTLZ.keySet()) {
			if (geneToBP.containsKey(gene)) {
				out.writeln(gene
						+ "\t" + geneToEQTLZ.get(gene)
						+ "\t" + geneToEQTLE.get(gene)
						+ "\t" + geneToBP.get(gene)
						+ "\t" + geneToPLI.get(gene)
						+ "\t" + geneToPLIZ.get(gene)
						+ "\t" + geneToSynZ.get(gene)
						+ "\t" + geneToMisZ.get(gene)
						
				);
			}
		}
		out.close();
	}
}
