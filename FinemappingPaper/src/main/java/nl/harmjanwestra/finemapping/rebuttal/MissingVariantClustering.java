package nl.harmjanwestra.finemapping.rebuttal;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.vcf.VCFTabix;
import nl.harmjanwestra.utilities.vcf.VCFVariant;

import java.io.IOException;

public class MissingVariantClustering {
	
	public static void main(String[] args) {
		
		String ref = "C:\\Data\\Ref\\1kg\\ALL.chrCHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";
		String samplefile = "C:\\Data\\Ref\\1kg-europeanpopulations.txt.gz";
		String out = "C:\\Data\\Ref\\1kg-maf\\allvars.txt";
		MissingVariantClustering v = new MissingVariantClustering();
		try {
			
			v.variantStats(ref, samplefile, out);
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run(){
		/*
		1. determine which veriants are included
		2. determine which variants are missing after imputation
		3. measure distance between variants in regions that are missing
		4. compare to random genomic snps with similar maf and LD?
		 */
	}
	
	public void variantStats(String ref, String samplefile, String out) throws IOException {
		TextFile tfo = new TextFile(out, TextFile.W);
		VCFTabix t = new VCFTabix();
		for (int c = 1; c < 23; c++) {
			String vcf = ref.replaceAll("CHR", "" + c);
			
			boolean[] samplestoinclude = null;
			if (samplefile != null) {
				samplestoinclude = t.getSampleFilter(samplefile, vcf);
			}
			
			TextFile tf = new TextFile(vcf, TextFile.R);
			System.out.println(vcf);
			String ln = tf.readLine();
			int ctr = 0;
			while (ln != null) {
				if (!ln.startsWith("#")) {
					VCFVariant v = new VCFVariant(ln, VCFVariant.PARSE.ALL, samplestoinclude);
					if (v.getMAF() > 0) {
						String outln = v.asSNPFeature().toString() + "\t" + v.getMAF();
						tfo.writeln(outln);
						ctr++;
						if (ctr % 10000 == 0) {
							System.out.print("\r" + ctr + " parsed");
						}
					}
				}
				ln = tf.readLine();
				
			}
			tf.close();
			
		}
		tfo.close();
		
	}
	
}
