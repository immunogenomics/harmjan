package nl.harmjanwestra.finemapping.rebuttal1;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import org.broad.igv.bbfile.*;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.stream.Stream;

public class R1Q4OVerlap {
	
	public static void main(String[] args) {

//		String[] dirs = new String[]{
//				"D:\\Data\\ENCODE\\histone_macs",
//				"D:\\Data\\ENCODE\\peakSeq",
//				"D:\\Data\\ENCODE\\spp\\hub",
//				"D:\\Data\\ENCODE\\wgEncodeAwgDnaseUniform",
//				"D:\\Data\\ENCODE\\wgEncodeAwgTfbsUniform"
//		};
//
//		String[] outfiles = new String[]{
//				"D:\\Data\\ENCODE\\histone_macs.txt",
//				"D:\\Data\\ENCODE\\peakSeq.txt",
//				"D:\\Data\\ENCODE\\spp.txt",
//				"D:\\Data\\ENCODE\\wgEncodeAwgDnaseUniform.txt",
//				"D:\\Data\\ENCODE\\wgEncodeAwgTfbsUniform.txt"
//		};

//		for (int d = 0; d < dirs.length; d++) {
//			try {
//				Stream<Path> list = Files.list(Paths.get(dirs[d]));
//
//				TextFile out = new TextFile(outfiles[d], TextFile.W);
//				Iterator<Path> paths = list.iterator();
//				while (paths.hasNext()) {
//					Path p = paths.next();
//
//					if (p.toString().endsWith(".bb") || p.toString().toLowerCase().endsWith("narrowpeak.gz")) {
//						out.writeln(p.toString());
//					}
//				}
//				out.close();
//
//
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//
//		}
		
		
		String[] maps = new String[]{
//				"D:\\Data\\ENCODE\\histone_macs.txt",
//				"D:\\Data\\ENCODE\\peakSeq.txt",
//				"D:\\Data\\ENCODE\\spp.txt",
//				"D:\\Data\\ENCODE\\wgEncodeAwgDnaseUniform.txt",
//				"D:\\Data\\ENCODE\\wgEncodeAwgTfbsUniform.txt",
				"D:\\Data\\GTRD\\g.txt"
		};
		
		String[] outs = new String[]{
//				"D:\\Data\\ENCODE\\histone_macs_overlap.txt",
//				"D:\\Data\\ENCODE\\peakSeq_overlap.txt",
//				"D:\\Data\\ENCODE\\spp_overlap.txt",
//				"D:\\Data\\ENCODE\\wgEncodeAwgDnaseUniform_overlap.txt",
//				"D:\\Data\\ENCODE\\wgEncodeAwgTfbsUniform_overlap.txt",
				"D:\\Data\\ENCODE\\GTRD_overlap.txt"
		};
		
		
		R1Q4OVerlap o = new R1Q4OVerlap();
		String variants = "D:\\Data\\ENCODE\\variants.txt";
		
		for (int m = 0; m < maps.length; m++) {
			
			try {
				o.overlap(variants, maps[m], outs[m]);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
	}
	
	
	public void overlap(String variantfile, String annotationFile, String outfile) throws IOException {
		
		ArrayList<Feature> snps = new ArrayList<Feature>();
		ArrayList<Pair<String, String>> annotations = new ArrayList<Pair<String, String>>();
		
		TextFile t1 = new TextFile(annotationFile, TextFile.R);
		String[] elems = t1.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length >= 2) {
				String desc = elems[0];
				String file = elems[1];
				annotations.add(new Pair<String, String>(desc, file));
			}
			
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
		ArrayList<String> celtypes = null;
		for (int i = 0; i < annotations.size(); i++) {
			
			String annot = annotations.get(i).getRight();
			
			System.out.println("Parsing: " + annotations.get(i).getLeft());
			if (annot.endsWith("interval.gz")) {
				
				HashMap<String, Integer> celtypemap = new HashMap<>();
				celtypes = new ArrayList<>();
				TextFile tf = new TextFile(annot, TextFile.R);
				tf.readLine();
				elems = tf.readLineElems(Strings.whitespace);
				int ctr = 0;
				while (elems != null) {
					if (!celtypemap.containsKey(elems[6])) {
						celtypemap.put(elems[6], ctr);
						celtypes.add(elems[6]);
						System.out.println(elems[6]);
						ctr++;
					}
					elems = tf.readLineElems(Strings.whitespace);
				}
				
				System.out.println(celtypes.size() + " cell types in total!");
				tf.close();
				overlap = new boolean[celtypes.size()][snps.size()];
				tf.open();
				
				tf.readLine();
				elems = tf.readLineElems(Strings.whitespace);
				while (elems != null) {
					Feature f = new Feature();
					f.setChromosome(Chromosome.parseChr(elems[0]));
					f.setStart(Integer.parseInt(elems[1]));
					f.setStop(Integer.parseInt(elems[2]));
					Integer celid = celtypemap.get(elems[6]);
					for (int s = 0; s < snps.size(); s++) {
						Feature snp = snps.get(s);
						if (snp.overlaps(f)) {
							overlap[celid][s] = true;
						}
					}
					elems = tf.readLineElems(Strings.whitespace);
				}
				tf.close();
				
				
			} else if (annotations.get(i).getRight().endsWith(".bb")) {
				
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
		for (int s = 0; s < snps.size(); s++) {
			header += "\t" + snps.get(s).toString();
		}
		out.writeln(header);
		if (celtypes != null) {
			for (int a = 0; a < celtypes.size(); a++) {
				String lnout = celtypes.get(a);
				boolean overlaponce = false;
				for (int s = 0; s < snps.size(); s++) {
					if (overlap[a][s]) {
						overlaponce = true;
						lnout += "\t1";
					} else {
						lnout += "\t0";
					}
				}
				if (overlaponce) {
					out.writeln(lnout);
				}
			}
		} else {
			for (int a = 0; a < annotations.size(); a++) {
				String lnout = annotations.get(a).getLeft();
				boolean overlaponce = false;
				for (int s = 0; s < snps.size(); s++) {
					if (overlap[a][s]) {
						overlaponce = true;
						lnout += "\t1";
					} else {
						lnout += "\t0";
					}
				}
				if (overlaponce) {
					out.writeln(lnout);
				}
			}
		}
		
		out.close();
		
	}
	
}
