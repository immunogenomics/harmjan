package nl.harmjanwestra.finemapping.rebuttal;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

public class SNPAnnotator {
	
	
	public static void main(String[] args) {
		String cyto = "D:\\Data\\cytoBand.txt.gz";
		String snps = "D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\2018-01-16-ListOfAssocIds.txt";
		SNPAnnotator q = new SNPAnnotator();
		try {
			q.run(cyto, snps);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String cyto, String snplist) throws IOException {
		
		ArrayList<Feature> bands = new ArrayList<Feature>();
		TextFile tf = new TextFile(cyto, TextFile.R);
		String ln = tf.readLine();
		while (ln != null) {
			String[] elems = ln.split("\t");
			Chromosome chr = Chromosome.parseChr(elems[0]);
			Integer start = Integer.parseInt(elems[1]);
			Integer stop = Integer.parseInt(elems[2]);
			Feature f = new Feature(chr, start, stop);
			f.setName(elems[3]);
			bands.add(f);
			ln = tf.readLine();
		}
		tf.close();
		
		ArrayList<SNPFeature> f = new ArrayList<>();
		tf = new TextFile(snplist, TextFile.R);
		ln = tf.readLine();
		while (ln != null) {
			SNPFeature feature = SNPFeature.parseSNPFeature(ln);
			
			f.add(feature);
			
			ln = tf.readLine();
		}
		tf.close();
		
		Collections.sort(f, new FeatureComparator(true));
		
		for (int i = 0; i < f.size(); i++) {
			Feature band = null;
			SNPFeature feat = f.get(i);
			for (int q = 0; q < bands.size(); q++) {
				if (bands.get(q).overlaps(feat)) {
					band = bands.get(q);
				}
			}
			
			System.out.println(feat.getChromosome().getNumber() + "" + band.getName() + "\t" + feat.getChromosome().toString() + "\t" + feat.getStart() + "\t" + feat.toString());
			
		}
		
	}
}
