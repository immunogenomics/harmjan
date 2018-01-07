package nl.harmjanwestra.finemappingtools;


import com.itextpdf.text.DocumentException;

import nl.harmjanwestra.finemappingtools.gwas.*;
import nl.harmjanwestra.finemappingtools.gwas.CLI.*;

import java.io.IOException;

/**
 * Created by hwestra on 11/23/15.
 */
public class Main {
	
	
	public static void main(String[] args) {
		
		try {
			MainOptions options = new MainOptions(args);
			if (options.mode.equals(MainOptions.MODE.NA)) {
				System.out.println("Please specify a mode");
				
			} else if (options.mode.equals(MainOptions.MODE.BEDFILTER)) {
				new BedAssocFilter(new BedAssocFilterOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.POSTERIORPVAL)) {
				new PosteriorPvalues(new PosteriorPvalueOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.PLOTPOSTERIORS)) {
				new AssociationPosteriorPlotter(new AssociationPlotterOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.CAVIAR)) {
				new Caviar(new CaviarOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.MERGE)) {
				new AssociationResultMerger(new AssociationResultMergerOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.UPDATERS)) {
				new ReannotateRSIds(new ReannotateRSIdsOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.FILTERASSOC)) {
				new AssociationResultFilter(new AssociationFilterOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.QTL)) {
				new QTLTest(new QTLTestOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.GUESS)) {
				new GUESS(new LRTestOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.FINEMAP)) {
				new FINEMAP(new LRTestOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.ASSOC)) {
				LRTestOptions optionsLR = new LRTestOptions(args);
				if (optionsLR.getAnalysisType().equals(LRTestOptions.ANALYSIS.HAPLOTYPE) ||
						optionsLR.getAnalysisType().equals(LRTestOptions.ANALYSIS.CONDITIONALHAPLOTYPE)) {
					new LRTestHaplotype(optionsLR);
				} else if (optionsLR.getAnalysisType().equals(LRTestOptions.ANALYSIS.STATS)) {
					new VariantStats(optionsLR);
				} else if (optionsLR.getAnalysisType().equals(LRTestOptions.ANALYSIS.FILTERLIST)) {
					new VariantFilterList(optionsLR);
				} else {
					new LRTest(optionsLR);
				}
			} else if (options.mode.equals(MainOptions.MODE.COUNTVARIANTS)) {
				new CountVariants(new CountVariantsOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.PLINK)) {
				new PLINKPCAConvert(args);
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
		
	}
	
	
}
