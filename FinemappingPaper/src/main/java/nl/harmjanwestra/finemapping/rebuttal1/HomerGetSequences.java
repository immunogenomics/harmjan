package nl.harmjanwestra.finemapping.rebuttal1;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class HomerGetSequences {
	
	public static void main(String[] args) {
		HomerGetSequences s = new HomerGetSequences();
		String variantfile = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\2018-01-16-ListOfAssocIds.txt";
		
		int window = 20;
		String output = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\homer\\bed\\";
		
		try {
			s.getBed(variantfile, window, output);
			output = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\homer\\seq\\";
			
			String fasta = "C:\\Data\\HumanGenome\\ucsc.hg19.fasta.gz";
//			s.getSeq(variantfile, fasta, 20, output);
			
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
			TextFile out = new TextFile(output + f.toString() + ".bed", TextFile.W);
			out.writeln(f.getChromosome().toString().toLowerCase() + "\t" + (f.getStart() - window) + "\t" + (f.getStop() + window) + "\t" + f.toString() + "\t-\t0");
			out.writeln(f.getChromosome().toString().toLowerCase() + "\t" + (f.getStart() - window) + "\t" + (f.getStop() + window) + "\t" + f.toString() + "\t-\t1");
			allbed.writeln(f.getChromosome().toString().toLowerCase() + "\t" + (f.getStart() - window) + "\t" + (f.getStop() + window) + "\t" + f.toString() + "\t-\t1");
//			allbed.writeln(f.getChromosome().toString().toLowerCase() + "\t" + (f.getStart() - window) + "\t" + (f.getStop() + window) + "\t" + f.toString() + "\t-\t0");
			
			out.close();
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
	
}
