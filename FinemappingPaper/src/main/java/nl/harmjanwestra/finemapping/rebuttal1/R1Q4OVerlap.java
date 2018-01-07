package nl.harmjanwestra.finemapping.rebuttal1;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import org.broad.igv.bbfile.*;

import java.io.IOException;
import java.util.ArrayList;

public class R1Q4OVerlap {
	
	public static void main(String[] args) {
	
		
	
	}
	
	
	public void overlap(String variantfile, String annotationFile, String outfile) throws IOException {
		
		ArrayList<Feature> snps = new ArrayList<Feature>();
		ArrayList<Pair<String, String>> annotations = new ArrayList<Pair<String, String>>();
		
		TextFile t1 = new TextFile(annotationFile, TextFile.R);
		String[] elems = t1.readLineElems(TextFile.tab);
		while (elems != null) {
			String desc = elems[0];
			String file = elems[1];
			annotations.add(new Pair<String, String>(desc, file));
			elems = t1.readLineElems(TextFile.tab);
		}
		t1.close();
		
		TextFile t2 = new TextFile(variantfile, TextFile.R);
		String ln = t2.readLine();
		while (ln != null) {
			SNPFeature f = SNPFeature.parseSNPFeature(ln);
			snps.add(f);
			ln = t2.readLine();
		}
		t2.close();
		
		boolean[][] overlap = new boolean[annotations.size()][snps.size()];
		
		for (int i = 0; i < annotations.size(); i++) {
			
			String annot = annotations.get(i).getRight();
			
			System.out.println("Parsing: " + annotations.get(i));
			
			if (annotations.get(i).getRight().endsWith(".bb")) {
				
				// bigwig file
				BBFileReader reader = new BBFileReader(annotations.get(i).getRight());
				
				
				//get the big header
				BBFileHeader bbFileHdr = reader.getBBFileHeader();
				if (!bbFileHdr.isHeaderOK()) {
					throw new IOException("bad header for " + annot);
				}
				//is it wig or bed ?
				if (!(bbFileHdr.isBigBed() || bbFileHdr.isBigWig())) {
					throw new IOException("undefined big type for " + annot);
				}
				
				
				if (bbFileHdr.isBigBed()) {
					BigBedIterator iter = reader.getBigBedIterator();
					while (iter.hasNext()) {
						BedFeature b = iter.next();
						Feature f = new Feature(Chromosome.parseChr(b.getChromosome()), b.getStartBase(), b.getEndBase());
						for (int s = 0; s < snps.size(); s++) {
							Feature snp = snps.get(s);
							if (snp.overlaps(f)) {
								overlap[i][s] = true;
							}
						}
					}
					
				} else if (bbFileHdr.isBigWig()) {
					BigWigIterator iter = reader.getBigWigIterator();
					while (iter.hasNext()) {
						WigItem b = iter.next();
						Feature f = new Feature(Chromosome.parseChr(b.getChromosome()), b.getStartBase(), b.getEndBase());
						for (int s = 0; s < snps.size(); s++) {
							Feature snp = snps.get(s);
							if (snp.overlaps(f)) {
								overlap[i][s] = true;
							}
						}
					}
				}
			} else {
				TextFile tf = new TextFile(annotations.get(i).getRight(), TextFile.R);
				elems = tf.readLineElems(Strings.whitespace);
				while (elems != null) {
					Feature f = new Feature();
					f.setChromosome(Chromosome.parseChr(elems[0]));
					f.setStart(Integer.parseInt(elems[1]));
					f.setStop(Integer.parseInt(elems[2]));
					for (int s = 0; s < snps.size(); s++) {
						Feature snp = snps.get(s);
						if (snp.overlaps(f)) {
							overlap[i][s] = true;
						}
					}
					elems = tf.readLineElems(Strings.whitespace);
				}
				tf.close();
			}
			
			
		}
		
		TextFile out = new TextFile(outfile, TextFile.W);
		
		String header = "Annotation";
		for (
				int s = 0; s < snps.size(); s++)
		
		{
			header += "\t" + snps.get(s).toString();
		}
		out.writeln(header);
		for (
				int a = 0; a < annotations.size(); a++)
		
		{
			String lnout = annotations.get(i).getLeft();
			for (int s = 0; s < snps.size(); s++) {
				if (overlap[a][s]) {
					lnout += "\t1";
				} else {
					lnout += "\t0";
				}
			}
			out.writeln(lnout);
		}
		out.close();
		
	}
	
}
