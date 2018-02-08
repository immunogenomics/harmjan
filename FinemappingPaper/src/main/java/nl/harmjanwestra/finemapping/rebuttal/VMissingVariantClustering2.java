package nl.harmjanwestra.finemapping.rebuttal;

import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;

public class VMissingVariantClustering2 {
	
	public static void main(String[] args) {
	
	
	}
	
	
	public void determineImputationOutputInfoScores(String[] refpanels,
													String[] ds,
													String[] dsinputfromIC,
													String[] dsnames,
													double mafthresholdref,
													double mafthresholdds,
													double infothreshold,
													String outputfileloc
	) {
		
		// gtcaller ds referencevariants referencevariantsnotonIC variantsinds info>0.8 maf>1% info>0.8+maf>1%
		for (int d = 0; d < ds.length; d++) {
			// get a list of variants present in this dataset
			
			KgVariant kgVariant = new KgVariant();
			
			
			// get a list of variants for the input
			
			for (int r = 0; r < refpanels.length; r++) {
				// get a list of variants present in this reference panel
				
				// determine how many variants total
				// determine how many variants maf>threshold
				// determine how many variants info>threshold
				// determine how many variants maf>threshold+info>threshold
				// repeat for variants not on IC
				
				
			}
		}
		
	}
	
	
	protected int idcol = 0;
	protected int minorallele1col = 1;
	protected int aleleles1col = 2;
	protected int maf1col = 3;
	protected int cr1col = 4;
	protected int minorallele2col = 5;
	protected int aleleles2col = 6;
	protected int maf2col = 7;
	protected int cr2col = 8;
	protected int dfcol = 9;
	protected int samplecol = 10;
	protected int rcol = 11;
	protected int rsqlcol = 12;
	protected int betacol = 13;
	protected int secol = 14;
	protected int impqual1 = 15;
	protected int impqual2 = 16;
	
	public void readCorrelationFile(String f) throws IOException {
		
		TextFile tf = new TextFile(f, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			
			if (elems.length > 1) {
				SNPFeature feat = SNPFeature.parseSNPFeature(elems[0]);
				String alleles = elems[2];
				String minor = elems[1];
				
				
				
			}
			
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
	}
}
