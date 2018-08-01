package nl.harmjanwestra.finemapping.rebuttal;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.io.trityper.util.BaseAnnot;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class FlipGenotypedVariants {
	
	public static void main(String[] args) {
		
		String[] assoc = new String[]{
				"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
				"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
		};
		String[] logloc = new String[]{
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\matchlogs\\RA-chrCHR-rsId1kg.vcf.gz-log.txt",
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\imputeinput\\matchlogs\\T1D-chrCHR-rsId1kg.vcf.gz-log.txt",
		};
		String[] out = new String[]{
				"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-flipgenotyped.txt.gz",
				"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-flipgenotyped.txt.gz",
		};
		
		FlipGenotypedVariants g = new FlipGenotypedVariants();
		try {
			g.run(assoc, logloc, out);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String[] in, String[] log, String[] out) throws IOException {
		for (int q = 0; q < in.length; q++) {
			HashMap<String, Boolean> flip = new HashMap<String, Boolean>();
			HashMap<String, Boolean> kill = new HashMap<String, Boolean>();
			HashMap<String, Boolean> complement = new HashMap<String, Boolean>();
			
			for (int chr = 1; chr < 23; chr++) {
				String logfile = log[q].replaceAll("CHR", "" + chr);
				if (Gpio.exists(logfile)) {
					System.out.println("Parsing: " + logfile);
					TextFile tf = new TextFile(logfile, TextFile.R);
					tf.readLine();
					String[] elems = tf.readLineElems(TextFile.tab);
					while (elems != null) {
						String snpstr = Chromosome.parseChr(elems[0]).toString() + "_" + elems[1] + "_" + elems[3];
						String conclusion = elems[elems.length - 1];
						if (conclusion.toLowerCase().contains("notok")) {
							kill.put(snpstr, true);
						}
						if (conclusion.toLowerCase().contains("flippedalleles")) {
							flip.put(snpstr, true);
						} else {
							flip.put(snpstr, false);
						}
						if (conclusion.toLowerCase().contains("complement")) {
							complement.put(snpstr, true);
						} else {
							complement.put(snpstr, false);
						}
						
						elems = tf.readLineElems(TextFile.tab);
					}
					tf.close();
				}
			}
			
			
			AssociationFile f = new AssociationFile();
			ArrayList<AssociationResult> assoc = f.read(in[q]);
			
			TextFile outf = new TextFile(out[q], TextFile.W);
			outf.writeln(f.getHeader());
			for (int a = 0; a < assoc.size(); a++) {
				AssociationResult r = assoc.get(a);
				
				String snpstr = r.getSnp().getChromosome().toString() + "_" + r.getSnp().getStart() + "_" + r.getSnp().getName();
				if (snpstr.contains("rs6698586")) {
					System.out.println("Found it");
				}
				Boolean flipsnp = flip.get(snpstr);
				Boolean complementsnp = complement.get(snpstr);
				if (flipsnp != null && flipsnp) {
					r.flip();
					
					
				}
				if (complementsnp != null && complementsnp) {
					String[] allelestmp = r.getSnp().getAlleles();
					String[] alleles = new String[allelestmp.length];
					for (int z = 0; z < alleles.length; z++) {
						alleles[z] = BaseAnnot.getComplement(allelestmp[z]);
					}
					r.getSnp().setAlleles(alleles);
					r.getSnp().setMinorAllele(BaseAnnot.getComplement(r.getSnp().getMinorAllele()));
				}
				if (!kill.containsKey(snpstr)) {
					outf.writeln(r.toString());
				}
			}
			outf.close();
			
			
		}
		
		
	}
	
	
}
