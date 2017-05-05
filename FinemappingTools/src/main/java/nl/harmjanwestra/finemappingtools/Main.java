package nl.harmjanwestra.finemappingtools;


import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.finemappingtools.broshifter.AnnotationOverlapPlot;
import nl.harmjanwestra.finemappingtools.broshifter.BroShifter;
import nl.harmjanwestra.finemappingtools.broshifter.CLI.BroShifterOptions;
import nl.harmjanwestra.finemappingtools.broshifter.CLI.GoShifterOptions;
import nl.harmjanwestra.finemappingtools.broshifter.CLI.MainOptions;
import nl.harmjanwestra.finemappingtools.broshifter.GoShifter;
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
			} else if (options.mode.equals(MainOptions.MODE.BROSHIFTER)) {
				new BroShifter(new BroShifterOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.GOSHIFTER)) {
				new GoShifter(new GoShifterOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.BEDFILTER)) {
				new BedAssocFilter(new BedAssocFilterOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.ANNOTATIONOVERLAPPLOT)) {
				new AnnotationOverlapPlot(new BroShifterOptions(args));
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
			} else if (options.mode.equals(MainOptions.MODE.QTL)) {
				new QTLTest(new QTLTestOptions(args));
			} else if (options.mode.equals(MainOptions.MODE.ASSOC)) {
				LRTestOptions optionsLR = new LRTestOptions(args);
				if (optionsLR.getAnalysisType().equals(LRTestOptions.ANALYSIS.HAPLOTYPE) ||
						optionsLR.getAnalysisType().equals(LRTestOptions.ANALYSIS.CONDITIONALHAPLOTYPE)) {
					new LRTestHaplotype(optionsLR);
				} else if (optionsLR.getAnalysisType().equals(LRTestOptions.ANALYSIS.STATS)) {
					new VariantStats(optionsLR);
				} else {
					new LRTest(optionsLR);
				}
			} else if (options.mode.equals(MainOptions.MODE.COUNTVARIANTS)) {
				new CountVariants(new CountVariantsOptions(args));
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}

	}


}
