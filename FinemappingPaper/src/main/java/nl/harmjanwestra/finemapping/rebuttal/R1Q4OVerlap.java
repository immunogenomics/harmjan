package nl.harmjanwestra.finemapping.rebuttal;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import org.broad.igv.bbfile.*;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.stream.Stream;

public class R1Q4OVerlap {
	
	public static void main(String[] args) {

//		String[] dirs = new String[]{
////				"D:\\Data\\ENCODE\\histone_macs",
////				"D:\\Data\\ENCODE\\peakSe
//////				"D:\\Data\\ENCODE\\spp\\hub",
////				"D:\\Data\\ENCODE\\wgEncodeAwgDnaseUniform",
////				"D:\\Data\\ENCODE\\wgEncodeAwgTfbsUniform"
//
//				"C:\\Data\\Enhancers\\CD4TimelinePilot",
//				"C:\\Data\\Enhancers\\ChromHMM\\12ImputedMarks",
//
//		};
//
//		String[] outfiles = new String[]{
////				"D:\\Data\\ENCODE\\histone_macs.txt",
////				"D:\\Data\\ENCODE\\peakSeq.txt",
////				"D:\\Data\\ENCODE\\spp.txt",
////				"D:\\Data\\ENCODE\\wgEncodeAwgDnaseUniform.txt",
////				"D:\\Data\\ENCODE\\wgEncodeAwgTfbsUniform.txt"
//		};

//
		
		R1Q4OVerlap o = new R1Q4OVerlap();
//
//		String[] dirs = new String[]{
//				"C:\\Data\\Enhancers\\FAMTOM5\\faced_expressed_enhancers",
//				"C:\\Data\\Enhancers\\FAMTOM5\\facet_differentially_expressed_0.05"
//		};
//		String[] outfiles = new String[]{
//				"C:\\Data\\Enhancers\\FAMTOM5\\enhancers_exp.txt",
//				"C:\\Data\\Enhancers\\FAMTOM5\\enhancers_diffex.txt"
//		};
//
//		o.lsof(dirs, outfiles);
//
//		System.exit(-1);
		String[] maps = new String[]{
				"C:\\Data\\ENCODE\\histone_macs.txt",
				"C:\\Data\\ENCODE\\peakSeq.txt",
				"C:\\Data\\ENCODE\\spp.txt",
				"C:\\Data\\ENCODE\\wgEncodeAwgDnaseUniform.txt",
				"C:\\Data\\ENCODE\\wgEncodeAwgTfbsUniform.txt",
				"C:\\Data\\Enhancers\\TrynkaEpigenetics\\hg19\\alldnase-windows.txt",
				"C:\\Data\\Enhancers\\ChromHMM\\allChromHMM-desc-windows.txt",
				"C:\\Data\\Enhancers\\FAMTOM5\\enhancers_diffex.txt",
				"C:\\Data\\Enhancers\\FAMTOM5\\enhancers_exp.txt",
				"C:\\Data\\Enhancers\\CD4TimelinePilot\\listwindows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\dnase-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H2A.Z-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H2AK5ac-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H2BK120ac-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H2BK12ac-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H2BK15ac-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H2BK20ac-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H2BK5ac-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H3K14ac-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H3K18ac-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H3K23ac-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H3K27ac-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\h3k27me3-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\h3k36me3-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H3K4ac-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\h3k4me1-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H3K4me2-groups-windows.txt",
				"C:\\Data\\Enhancers\\TrynkaEpigenetics\\hg19\\h3k4me3files-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H3K79me1-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H3K79me2-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\h3k9ac-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\h3k9me3-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H4K12ac-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H4K20me1-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H4K5ac-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H4K8ac-groups-windows.txt",
				"C:\\Data\\Enhancers\\Roadmap\\H4K91ac-groups-windows.txt"
//				"D:\\Data\\GTRD\\g.txt"
		};
		
		String[] outs = new String[]{
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\histone_macs.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\peakSeq.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\spp.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\wgEncodeAwgDnaseUniform.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\wgEncodeAwgTfbsUniform.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\trynkadnase.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\allChromHMM-desc-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\enhancers_diffex.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\enhancers_exp.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\listwindows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\dnase-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H2A.Z-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H2AK5ac-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H2BK120ac-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H2BK12ac-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H2BK15ac-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H2BK20ac-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H2BK5ac-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H3K14ac-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H3K18ac-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H3K23ac-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H3K27ac-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\h3k27me3-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\h3k36me3-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H3K4ac-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\h3k4me1-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H3K4me2-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\h3k4me3-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H3K79me1-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H3K79me2-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\h3k9ac-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\h3k9me3-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H4K12ac-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H4K20me1-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H4K5ac-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H4K8ac-groups-windows.txt",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\H4K91ac-groups-windows.txt"
//				"D:\\Data\\ENCODE\\GTRD_overlap.txt"
		};
		
		maps = new String[]{
				"d:\\Data\\CD4TimelinePilot\\listwindows.txt"
		};
		outs = new String[]{
				"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\enhanceroverlap\\atac.txt"
		};
		
		String variants = "d:\\Sync\\Dropbox\\FineMap\\2018-01-Rebuttal\\tables\\listofsnpswithposterior0.2.txt";
		
		for (int m = 0; m < maps.length; m++) {
			try {
				o.overlap(variants, maps[m], 0, outs[m]);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
	}
	
	
	public void lsof(String[] dirs, String[] outfiles) {
		for (int d = 0; d < dirs.length; d++) {
			try {
				Stream<Path> list = Files.list(Paths.get(dirs[d]));
				
				TextFile out = new TextFile(outfiles[d], TextFile.W);
				Iterator<Path> paths = list.iterator();
				while (paths.hasNext()) {
					Path p = paths.next();
					
					if (p.toString().endsWith(".bed.gz") || p.toString().endsWith(".bb") || p.toString().toLowerCase().endsWith("narrowpeak.gz")) {
						out.writeln(p.toString());
					}
				}
				out.close();
				
				
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		}
		
	}
	
	public void overlap(String variantfile, String annotationFile, int window, String outfile) throws IOException {
		
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
		
		System.out.println(annotations.size() + "  annotations.");
		TextFile t2 = new TextFile(variantfile, TextFile.R);
		String ln = t2.readLine();
		while (ln != null) {
			SNPFeature f = SNPFeature.parseSNPFeature(ln);
			f.setStart(f.getStart() - window);
			f.setStop(f.getStop() + window);
			snps.add(f);
			ln = t2.readLine();
		}
		t2.close();
		
		Collections.sort(snps, new FeatureComparator(true));
		System.out.println(snps.size() + " snps");
		
		boolean[][] overlap = new boolean[annotations.size()][snps.size()];
		boolean[] overlapsAny = new boolean[snps.size()];
		ArrayList<String> celtypes = null;
		for (int i = 0; i < annotations.size(); i++) {
			
			String annot = annotations.get(i).getRight();
			if (Gpio.exists(annot)) {
				System.out.println("Parsing: " + annotations.get(i).getLeft() + "\t" + annotations.get(i).getRight());
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
								overlapsAny[s] = true;
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
									overlapsAny[s] = true;
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
									overlapsAny[s] = true;
								}
							}
						}
					}
				} else {
					TextFile tf = new TextFile(annotations.get(i).getRight(), TextFile.R);
					
					
					elems = tf.readLineElems(Strings.whitespace);
					while (elems != null) {
						boolean parseln = true;
						if (tf.getFileName().endsWith(".xls") || tf.getFileName().endsWith(".xls.gz")) {
							if (elems[0].equals("chr")) {
								parseln = false;
							}
						}
						if (elems[0].startsWith("#") || elems.length < 3) {
							parseln = false;
						}
						if (parseln) {
							Feature f = new Feature();
							f.setChromosome(Chromosome.parseChr(elems[0]));
							f.setStart(Integer.parseInt(elems[1]));
							f.setStop(Integer.parseInt(elems[2]));
							for (int s = 0; s < snps.size(); s++) {
								Feature snp = snps.get(s);
								if (snp.overlaps(f)) {
									overlap[i][s] = true;
									overlapsAny[s] = true;
								}
							}
						}
						
						elems = tf.readLineElems(Strings.whitespace);
					}
					tf.close();
				}
			}
		}
		
		TextFile out = new TextFile(outfile, TextFile.W);
		
		String header = "Annotation";
		for (int s = 0; s < snps.size(); s++) {
			if (overlapsAny[s]) {
				header += "\t" + snps.get(s).toString();
			}
		}
		out.writeln(header);
		if (celtypes != null) {
			for (int a = 0; a < celtypes.size(); a++) {
				String lnout = celtypes.get(a);
				boolean overlaponce = false;
				for (int s = 0; s < snps.size(); s++) {
					if (overlapsAny[s]) {
						if (overlap[a][s]) {
							overlaponce = true;
							lnout += "\t1";
						} else {
							lnout += "\t0";
						}
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
					if (overlapsAny[s]) {
						if (overlap[a][s]) {
							overlaponce = true;
							lnout += "\t1";
						} else {
							lnout += "\t0";
						}
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
