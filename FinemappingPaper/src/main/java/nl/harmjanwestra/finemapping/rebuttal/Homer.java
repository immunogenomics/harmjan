package nl.harmjanwestra.finemapping.rebuttal;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.enums.Strand;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class Homer {
	
	public static void main(String[] args) {
		Homer s = new Homer();
		String variantfile = "C:\\Sync\\Dropbox\\FineMap\\2018-01-Rebuttal\\tables\\listofsnpswithposterior0.2.txt";
		
		int window = 20;
		String output = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\homer\\bed\\";
		
		try {
			s.getBed(variantfile, window, output);
//			output = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\homer\\seq\\";
//
//			String fasta = "C:\\Data\\HumanGenome\\ucsc.hg19.fasta.gz";
//			s.getSeq(variantfile, fasta, 20, output);
			
			String filein = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\homer\\bed\\output-text.txt";
			String fileout = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\homer\\bed\\output-text-overlapping.txt";
			s.parseOutput(filein, fileout, 40, Strand.POS);
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void getBed(String variantfile, int window, String output) throws IOException {
		TextFile tf = new TextFile(variantfile, TextFile.R);
		TextFile allout = new TextFile(output + "all.txt", TextFile.W);
		TextFile allbed = new TextFile(output + "all.bed", TextFile.W);
		String ln = tf.readLine();
		
		ArrayList<SNPFeature> feats = new ArrayList<>();
		while (ln != null) {
			SNPFeature f = SNPFeature.parseSNPFeature(ln);
			
			if (feats == null) {
				feats = new ArrayList<>();
			}
			feats.add(f);
			ln = tf.readLine();
		}
		
		
		Collections.sort(feats, new FeatureComparator(true));
		
		for (SNPFeature f : feats) {
//			TextFile out = new TextFile(output + f.toString() + ".bed", TextFile.W);
//			out.writeln(f.getChromosome().toString().toLowerCase() + "\t" + (f.getStart() - window) + "\t" + (f.getStop() + window) + "\t" + f.toString() + "\t-\t0");
//			out.writeln(f.getChromosome().toString().toLowerCase() + "\t" + (f.getStart() - window) + "\t" + (f.getStop() + window) + "\t" + f.toString() + "\t-\t1");
			allbed.writeln(f.getChromosome().toString().toLowerCase() + "\t" + (f.getStart() - window) + "\t" + (f.getStop() + window) + "\t" + f.toString() + "\t+\t+");
//			allbed.writeln(f.getChromosome().toString().toLowerCase() + "\t" + (f.getStart() - window) + "\t" + (f.getStop() + window) + "\t" + f.toString() + "\t-\t0");

//			out.close();
			allout.writeln("mkdir -p ./output/" + f.toString() + "/");
			allout.writeln("findMotifsGenome.pl ./" + f.toString() + ".bed hg19 ./output/" + f.toString() + "/ -p 4 -nomotif");
			
		}
		
		allbed.close();
		allout.close();
		tf.close();
		
		System.out.println("Done writing bed files");
		
	}
	
	public void getSeq(String variantfile, String genomeFasta, int window, String output) throws IOException {
		
		
		TextFile tf = new TextFile(variantfile, TextFile.R);
		String ln = tf.readLine();
		HashMap<Chromosome, ArrayList<SNPFeature>> snpsPerchr = new HashMap<Chromosome, ArrayList<SNPFeature>>();
		while (ln != null) {
			SNPFeature f = SNPFeature.parseSNPFeature(ln);
			ArrayList<SNPFeature> feats = snpsPerchr.get(f.getChromosome());
			if (feats == null) {
				feats = new ArrayList<>();
			}
			feats.add(f);
			snpsPerchr.put(f.getChromosome(), feats);
			ln = tf.readLine();
		}
		tf.close();
		
		TextFile allout = new TextFile(output + "all.txt", TextFile.W);
		TextFile allseq = new TextFile(output + "all.fa", TextFile.W);
		
		htsjdk.samtools.reference.FastaSequenceFile fastaFile = new htsjdk.samtools.reference.FastaSequenceFile(new File(genomeFasta), false);
		htsjdk.samtools.reference.ReferenceSequence seq = fastaFile.nextSequence();
		
		while (seq != null) {
			Chromosome chr = Chromosome.parseChr(seq.getName());
			System.out.println(seq.getName() + "\t" + chr.toString());
			byte[] bases = seq.getBases();
			
			ArrayList<SNPFeature> feats = snpsPerchr.get(chr);
			if (feats != null) {
				for (SNPFeature feat : feats) {
					int sta = feat.getStart() - window;
					int sto = feat.getStop() + window;
					if (sta < 0) {
						sta = 0;
					}
					if (sto > bases.length - 1) {
						sto = bases.length - 1;
					}
					byte[] select = new byte[sto - sta];
					System.arraycopy(bases, sta, select, 0, sto - sta);
					String sequence = new String(select).toUpperCase();
					TextFile out = new TextFile(output + feat.toString() + ".fa", TextFile.W);
					out.writeln(">" + feat.toString());
					out.writeln(sequence);
					out.close();
					allout.writeln(output + feat.toString() + ".fa");
					
					allseq.writeln(">" + feat.toString());
					allseq.writeln(sequence);
					System.out.println(output + feat.toString() + ".fa");
					
				}
			}
			
			seq = fastaFile.nextSequence();
		}
		
		allout.close();
		allseq.close();
		fastaFile.close();
	}
	
	
	public void parseOutput(String input, String output, int len, Strand snpStrand) throws IOException {
		TextFile tf = new TextFile(input, TextFile.R);
		TextFile out = new TextFile(output, TextFile.W);
		out.writeln(tf.readLine());
		
		String[] elems = tf.readLineElems(TextFile.tab);
		HashSet<String> uniquesnps = new HashSet<String>();
		
		
		ArrayList<HomerPair> outp = new ArrayList<>();
		
		while (elems != null) {
			// PositionID	Offset	Sequence	Motif	Name	Strand	MotifScore
			if (elems.length > 5) {
				SNPFeature snp = SNPFeature.parseSNPFeature(elems[0]);
				
				Integer offset = Integer.parseInt(elems[1]);
				String seq = elems[2];
				String motif = elems[3];
				Strand motifStrand = Strand.parseStr(elems[4]);
				Double score = Double.parseDouble(elems[5]);
				
				// determine overlap
				boolean overlap = false;
				if (elems[0].contains("rs117701653")) {
					System.out.println("found");
				}
				if (snpStrand.equals(Strand.POS)) {
					int motiflen = seq.length();
					
					if (motifStrand.equals(Strand.POS)) {
						// do stuff
						int snppos = 20;
						int motifstart = offset;
						int motifend = offset + motiflen;
						if (snppos >= motifstart && snppos <= motifend) {
							overlap = true;
						}
					} else {
						// do stuff differently
						int snppos = 20;
						int motifend = offset;
						int motifstart = offset - motiflen;
						if (snppos >= motifstart && snppos <= motifend) {
							overlap = true;
						}
					}
				}
				if (overlap) {
					uniquesnps.add(elems[0]);
					HomerPair p = new HomerPair();
					p.s = snp;
					p.ln = Strings.concat(elems, Strings.tab);
					outp.add(p);
				}
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		
		
		tf.close();
		
		Collections.sort(outp);
		
		for (HomerPair p : outp) {
			out.writeln(p.ln);
		}
		out.writeln(Strings.concat(elems, Strings.tab));
		out.close();
		System.out.println(uniquesnps.size() + " overlapping variants");
	}
	
	class HomerPair implements Comparable<HomerPair> {
		SNPFeature s;
		String ln;
		
		
		@Override
		public boolean equals(Object o) {
			if (this == o) return true;
			if (o == null || getClass() != o.getClass()) return false;
			
			HomerPair homerPair = (HomerPair) o;
			
			return s != null ? s.equals(homerPair.s) : homerPair.s == null;
		}
		
		@Override
		public int hashCode() {
			return s != null ? s.hashCode() : 0;
		}
		
		@Override
		public int compareTo(HomerPair o) {
			FeatureComparator c = new FeatureComparator(true);
			return c.compare(this.s, o.s);
		}
	}
}
