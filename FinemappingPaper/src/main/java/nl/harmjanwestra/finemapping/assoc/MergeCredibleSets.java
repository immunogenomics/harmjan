package nl.harmjanwestra.finemapping.assoc;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.finemapping.annotation.EQTL;
import nl.harmjanwestra.finemappingtools.gwas.PosteriorPvalues;
import nl.harmjanwestra.utilities.annotation.Annotation;
import nl.harmjanwestra.utilities.annotation.ensembl.EnsemblStructures;
import nl.harmjanwestra.utilities.annotation.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.AssociationResultPValuePair;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.enums.SNPClass;
import nl.harmjanwestra.utilities.enums.Strand;
import nl.harmjanwestra.utilities.features.*;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.panels.AnnotationTrackPanel;
import nl.harmjanwestra.utilities.graphics.panels.AssociationPanel;
import nl.harmjanwestra.utilities.graphics.panels.CircularHeatmapPanel;
import nl.harmjanwestra.utilities.graphics.panels.GenePanel;
import nl.harmjanwestra.utilities.legacy.genetica.io.Gpio;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFTabix;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Triple;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.lang.reflect.Array;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by hwestra on 9/7/16.
 */
public class MergeCredibleSets {
	
	static boolean windows = false;
	static int promotordistance = 1000;
	
	public static void main(String[] args) {
		
		
		try {
			MergeCredibleSets c = new MergeCredibleSets();


//			c.determineRegionSignificanceThresholds(bedregions, assocfiles, datasetnames, genenames, outfile);
			
			String bedregions = "d:/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
//			bedregions = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/overlapplots/regions.bed";
//		String bedregions = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2016-09-06-SummaryStats/NormalHWEP1e4/RA-significantloci-75e7.bed";
			String genenames = "d:/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/AllLoci-GenesPerLocus.txt";
			String geneAnnotation = "d:/Data/Ref/Ensembl/GrCH37-b86-Structures.txt.gz"; //"/Data/Ref/Annotation/UCSC/genes.gtf.gz";
			
			String outdir = "d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\normal\\";
			String outfile = outdir + "mergedCredibleSets.txt";
			String outeqtlfile = outdir + "mergedCredibleSets-eqtls.txt";
			String outeqtlfileblueprint = outdir + "mergedCredibleSets-eqtls-blueprint";
			String outpqtlfile = outdir + "mergedCredibleSets-pqtls-sun";
			String outoverlapfile = outdir + "mergedCredibleSets-overlap.txt";
			String outoverlapplot = outdir + "annotationplots/mergedCredibleSets-overlapplot-";
			String outoverlapfileatac = outdir + "mergedCredibleSets-overlap-atac.txt";
			String outoverlapplotatac = outdir + "annotationplots/mergedCredibleSets-overlapplot-atac-";
			String outplot = outdir + "mergedCredibleSets-plot-promoter" + promotordistance + ".pdf";
			String outgtex = outdir + "mergedCredibleSets-eqtls-gtex.txt";

//			String[] assocfiles = new String[]{
//					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-15-ReImpute3\\normal\\1kg-RA\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
//					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-15-ReImpute3\\normal\\1kg-T1D\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
//
//			};
			
			
			double significanceThreshold = 7.5E-7;
			int nrVariantsInCredibleSet = 10;
			double maxPosteriorCredibleSet = 0.95;
			boolean includeAllLoci = true;
			
			String[] assocfiles = new String[]{
					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-03-25-SummaryStats\\normal\\mergedcrediblesets\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\normal\\1kg-RA\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\normal\\1kg-PBWT-RA-wGenotypes\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\normal\\HRC-RA\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-03-25-SummaryStats\\normal\\mergedcrediblesets\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\normal\\1kg-T1D\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\normal\\1kg-PBWT-T1D-wGenotypes\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\normal\\HRC-T1D\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-03-25-SummaryStats\\normal\\mergedcrediblesets\\META-assoc0.3-COSMO-merged-posterior.txt.gz",
					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\normal\\1kg-META\\META-assoc0.3-COSMO-merged-posterior.txt.gz",
					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\normal\\1kg-PBWT-META-wGenotypes\\META-assoc0.3-COSMO-merged-posterior.txt.gz",
			};
			String[] datasetnames = new String[]{
					"1kg-RA-orig",
					"1kg-RA",
					"1kg-PBWT-RA",
					"HRC-RA",
					"1kg-T1D-orig",
					"1kg-T1D",
					"1kg-PBWT-T1D",
					"HRC-T1D",
					"1kg-META-orig",
					"1kg-META",
					"1kg-PBWT-META"
			};
			
			// merge in TYK2
			assocfiles = new String[]{
					"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\normal\\1kg-RA\\RA-assoc0.3-COSMO-tyk2-chr19-gwas-0.txt",
					"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\normal\\1kg-T1D\\T1D-assoc0.3-COSMO-tyk2-chr19-gwas-0.txt",
					"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\normal\\1kg-META\\META-assoc0.3-COSMO-tyk2-chr19-gwas-0.txt"
			};
			String[] assocfiles2 = new String[]{
					"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\normal\\1kg-RA\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
					"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\normal\\1kg-T1D\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
					"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\normal\\1kg-META\\META-assoc0.3-COSMO-merged-posterior.txt.gz"
			};

//			c.replaceLocus(assocfiles, assocfiles2, new Feature(Chromosome.NINETEEN, 10396336, 10628468), bedregions, maxPosteriorCredibleSet);
			
			
			assocfiles = new String[]{
					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\normal\\1kg-RA\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\normal\\1kg-T1D\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\normal\\1kg-META\\META-assoc0.3-COSMO-merged-posterior.txt.gz"
			};
			
			
			assocfiles = new String[]{
					"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
					"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
					"d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\META-assoc0.3-COSMO-merged-posterior.txt.gz"
			};
			datasetnames = new String[]{
					"1kg-RA",
					"1kg-T1D",
					"1kg-META"
				
			};
			
			outfile += "missp.txt";
			
			c.mergeCredibleSets(bedregions,
					assocfiles,
					datasetnames,
					genenames,
					outfile,
					maxPosteriorCredibleSet,
					significanceThreshold,
					nrVariantsInCredibleSet,
					geneAnnotation,
					includeAllLoci);
			System.out.println(outfile + "missp.txt");
//
			boolean onlyincludevariantsbelowsignificancethreshold = true;
			c.makeCircularPlot(bedregions,
					assocfiles,
					datasetnames,
					genenames,
					outplot,
					significanceThreshold,
					onlyincludevariantsbelowsignificancethreshold,
					maxPosteriorCredibleSet,
					nrVariantsInCredibleSet,
					geneAnnotation);
			System.exit(-1);
//
////			System.exit(-1);
//
////
////			/*
////			eQTL overlap
////		*/
			
			String[] eqtlfiles = new String[]{
					"/Data/eQTLs/ImmVar/Raj/tableS12_meta_cd4T_cis_fdr05-upd.tab",
//					"/Data/eQTLs/Milani/CD4-cis-eQTLs-ProbeLevelFDR0.5.txt.gz",
//					"/Data/eQTLs/Milani/CD4-trans-eQTLs-ProbeLevelFDR0.5.txt.gz",
//					"/Data/eQTLs/Milani/CD8-cis-eQTLs-ProbeLevelFDR0.5.txt.gz",
//					"/Data/eQTLs/Milani/CD8-trans-eQTLs-ProbeLevelFDR0.5.txt.gz",
					"/Data/eQTLs/BiosEQTLs/eQTLsFDR0.05-ProbeLevel.txt.gz"
			};
			String[] eqtlfilenames = new String[]{
					"Raj",
//					"Milani-CD4",
//					"Milani-CD4-trans",
//					"Milani-CD8",
//					"Milani-CD8-trans",
					"Bios"
			};
////
////			String dbsnpvcf = "/Data/Ref/dbSNP/human_9606_b147_GRCh37p13/00-All.vcf.gz";
////			c.rewriteEQTLFile("/Data/eQTLs/ImmVar/Raj/tableS12_meta_cd4T_cis_fdr05.tsv", dbsnpvcf,
////					"/Data/eQTLs/ImmVar/Raj/tableS12_meta_cd4T_cis_fdr05-upd.tab");
////
			String tabixprefix = "/Data/Ref/beagle_1kg/1kg.phase3.v5a.chrCHR.vcf.gz";
			String tabixfilter = "/Data/Ref/1kg-europeanpopulations.txt.gz";
//			c.eQTLOverlap(bedregions, eqtlfiles, eqtlfilenames, tabixprefix, tabixfilter, assocfiles, datasetnames, genenames, outeqtlfile, maxPosteriorCredibleSet, nrVariantsInCredibleSet, geneAnnotation);
////			System.exit(-1);
			eqtlfilenames = new String[]{
					"Monocyte-eQTL",
					"Monocyte-hQTL-K27AC",
					"Monocyte-hQTL-K4ME1",
					"Monocyte-mQTL",
					"Monocyte-sQTL",
			};
			
			eqtlfiles = new String[]{
					"/Data/eQTLs/BluePrint/mono_gene_nor_combat_peer_10_all_summary-fdr005.tab.gz",
					"/Data/eQTLs/BluePrint/mono_K27AC_log2rpm_peer_10_all_summary-fdr005.tab.gz",
					"/Data/eQTLs/BluePrint/mono_K4ME1_log2rpm_peer_10_all_summary-fdr005.tab.gz",
					"/Data/eQTLs/BluePrint/mono_meth_M_peer_10_all_summary-fdr005.tab.gz",
					"/Data/eQTLs/BluePrint/mono_psi_peer_10_all_summary-fdr005.tab.gz",
			};
//			c.eQTLOverlap(bedregions, eqtlfiles, eqtlfilenames, tabixprefix, tabixfilter, assocfiles, datasetnames, genenames, outeqtlfileblueprint + "-monocyte.txt", maxPosteriorCredibleSet, nrVariantsInCredibleSet, geneAnnotation);
			eqtlfilenames = new String[]{
					"Neutrophil-eQTL",
					"Neutrophil-hQTL-K27AC",
					"Neutrophil-hQTL-K4ME1",
					"Neutrophil-mQTL",
					"Neutrophil-sQTL",
			};
			
			eqtlfiles = new String[]{
					"/Data/eQTLs/BluePrint/neut_gene_nor_combat_peer_10_all_summary-fdr005.tab.gz",
					"/Data/eQTLs/BluePrint/neut_K27AC_log2rpm_peer_10_all_summary-fdr005.tab.gz",
					"/Data/eQTLs/BluePrint/neut_K4ME1_log2rpm_peer_10_all_summary-fdr005.tab.gz",
					"/Data/eQTLs/BluePrint/neut_meth_M_peer_10_all_summary-fdr005.tab.gz",
					"/Data/eQTLs/BluePrint/neut_psi_peer_10_all_summary-fdr005.tab.gz",
			};
//			c.eQTLOverlap(bedregions, eqtlfiles, eqtlfilenames, tabixprefix, tabixfilter, assocfiles, datasetnames, genenames, outeqtlfileblueprint + "-neutrophils.txt", maxPosteriorCredibleSet, nrVariantsInCredibleSet, geneAnnotation);
			
			eqtlfilenames = new String[]{
					"TCell-eQTL",
					"TCell-hQTL-K27AC",
					"TCell-hQTL-K4ME1",
					"TCell-mQTL",
					"TCell-sQTL"
				
			};
			
			eqtlfiles = new String[]{
					"/Data/eQTLs/BluePrint/tcel_gene_nor_combat_peer_10_all_summary-fdr005.tab.gz",
					"/Data/eQTLs/BluePrint/tcel_K27AC_log2rpm_peer_10_all_summary-fdr005.tab.gz",
					"/Data/eQTLs/BluePrint/tcel_K4ME1_log2rpm_peer_10_all_summary-fdr005.tab.gz",
					"/Data/eQTLs/BluePrint/tcel_meth_M_peer_10_all_summary-fdr005.tab.gz",
					"/Data/eQTLs/BluePrint/tcel_psi_peer_10_all_summary-fdr005.tab.gz"
			};
//			c.eQTLOverlap(bedregions, eqtlfiles, eqtlfilenames, tabixprefix, tabixfilter, assocfiles, datasetnames, genenames, outeqtlfileblueprint + "-tcell.txt", maxPosteriorCredibleSet, nrVariantsInCredibleSet, geneAnnotation);
			
			// pQTL
//			c.rewritepQTLFiles();
			eqtlfilenames = new String[]{
					"Sun-pQTL",
					"Sun-pQTL-Conditional",
			};
			eqtlfiles = new String[]{
					
					"/Data/eQTLs/Sun-pQTL/table1.tab.gz",
					"/Data/eQTLs/Sun-pQTL/table2.tab.gz"
			};
			c.eQTLOverlap(bedregions, eqtlfiles, eqtlfilenames, tabixprefix, tabixfilter, assocfiles, datasetnames, genenames, outpqtlfile + "-tcell.txt", maxPosteriorCredibleSet, nrVariantsInCredibleSet, geneAnnotation);
			
			
			// GTEX

//			c.rewriteGTEX();
			eqtlfilenames = new String[]{
					"GTEX",
			};
			eqtlfiles = new String[]{
					"/Data/eQTLs/GTEx/v6peQTLsTab.txt"
			};
			boolean eQTLFileIsList = true;
			c.eQTLOverlapGTEX(bedregions, eqtlfiles, eqtlfilenames, tabixprefix, tabixfilter, assocfiles, datasetnames, genenames, outgtex, maxPosteriorCredibleSet, nrVariantsInCredibleSet, geneAnnotation);
//
//			/*
//			Bed file overlap
//		*/
//
//
			String[] bedfiles = new String[]{
					"/Data/Enhancers/Roadmap/dnase-groups.txt",
					"/Data/Enhancers/Roadmap/h3k4me3-groups.txt",
					"/Data/Enhancers/ChromHMM/ChromHMMEnhancers-groups.txt",
					"/Data/Enhancers/ChromHMM/ChromHMMPromotors-groups.txt"
			};
			String[] bedfilenames = new String[]{
					"DNASE", "H3K4me3", "ChromHMM-Enhancers", "ChromHMM-Promotors"
			};
//
//
//			c.bedOverlap(bedregions,
//					bedfiles,
//					bedfilenames,
//					assocfiles, datasetnames, genenames,
//					outoverlapfile, outoverlapplot,
//					maxPosteriorCredibleSet, nrVariantsInCredibleSet, geneAnnotation, false);
//
			bedfiles = new String[]{
					"/Data/Enhancers/CD4TimelinePilot/list.txt"
			};
			bedfilenames = new String[]{
					"Atac"
			};

//			c.bedOverlap(bedregions,
//					bedfiles,
//					bedfilenames,
//					assocfiles, datasetnames, genenames,
//					outoverlapfileatac, outoverlapplotatac,
//					maxPosteriorCredibleSet, nrVariantsInCredibleSet, geneAnnotation, false);
//
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void replaceLocus(String[] assocfiles, String[] assocfiles2, Feature feature, String regionfile, double bayesthreshold) throws IOException {
		
		for (int q = 0; q < assocfiles2.length; q++) {
			if (Gpio.exists(assocfiles[q]) && Gpio.exists(assocfiles2[q])) {
				AssociationFile f = new AssociationFile();
				ArrayList<AssociationResult> results1 = f.read(assocfiles2[q]);
				ArrayList<AssociationResult> filteredResults = new ArrayList<>();
				System.out.println(results1.size() + " results in " + assocfiles2[q]);
				for (AssociationResult r : results1) {
					if (!r.getSnp().overlaps(feature)) {
						filteredResults.add(r);
					}
				}
				System.out.println(filteredResults.size() + " results after filtering");
				ArrayList<AssociationResult> results2 = f.read(assocfiles[q]);
				// calculate posteriors
				ApproximateBayesPosterior abf = new ApproximateBayesPosterior();
				abf.calculatePosterior(results2);
				filteredResults.addAll(results2);
				System.out.println(filteredResults.size() + " results after merging.");
				// overwrite original
				TextFile tf = new TextFile(assocfiles2[q] + "_tmp", TextFile.W);
				tf.writeln(f.getHeader());
				for (AssociationResult r : filteredResults) {
					tf.writeln(r.toString());
				}
				tf.close();
				PosteriorPvalues post = new PosteriorPvalues();
				post.determinePosteriors(assocfiles2[q] + "_tmp", regionfile, assocfiles2[q], bayesthreshold);
				Gpio.delete(assocfiles2[q] + "_tmp");
			}
		}
		
		
	}
	
	private void rewritepQTLFiles() throws IOException {
		String pQTLFile1 = "/Data/eQTLs/Sun-pQTL/table1.txt";
		String pQTLFile1Out = "/Data/eQTLs/Sun-pQTL/table1.tab.gz";
		
		// in:
		// Variant	chr	pos	regionsta	regionsto	al1	al2	EAF	INFO	SomaID	Target	FullName	UniProt	cistrans	gene	coh1beta	co1se	co1p	coh2beta	co2se	co2p	metabeta	metase	metap	nrConditional	prevreported	locusid
		// out:
		// Chr1    168921757       rs12095636      ENSG00000000460.11      1.047e-02
		
		TextFile tf = new TextFile(pQTLFile1, TextFile.R);
		TextFile tfout = new TextFile(pQTLFile1Out, TextFile.W);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String pos = elems[2].replaceAll(",", "");
			pos = pos.replaceAll("\"", "");
			String out = elems[1]
					+ "\t" + pos
					+ "\t" + elems[0]
					+ "\t" + elems[12] + "/" + elems[10]
					+ "\t" + elems[23];
			tfout.writeln(out);
			elems = tf.readLineElems(TextFile.tab);
		}
		
		tf.close();
		tfout.close();
		
		// input:
		// SOMAmer ID	Target	UniProt	Sentinel variant	Conditional variant	LD with sentinel variant (r2)	Conditional variant Chr	Conditional variant Pos	Conditional variant EA	Conditional variant OA	Conditional variant EAF	Conditional variant INFO	Joint model beta	Joint model SE	Joint model p	Univariate beta	Univariate SE	Univariate p
		String pQTLFile2 = "/Data/eQTLs/Sun-pQTL/table2.txt";
		String pQTLFile2Out = "/Data/eQTLs/Sun-pQTL/table2.tab.gz";
		tf = new TextFile(pQTLFile2, TextFile.R);
		tfout = new TextFile(pQTLFile2Out, TextFile.W);
		tf.readLine();
		elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String pos = elems[7].replaceAll(",", "");
			pos = pos.replaceAll("\"", "");
			String out = elems[6]
					+ "\t" + pos
					+ "\t" + elems[4]
					+ "\t" + elems[2] + "/" + elems[1]
					+ "\t" + elems[14];
			tfout.writeln(out);
			elems = tf.readLineElems(TextFile.tab);
		}
		
		tf.close();
		tfout.close();
	}
	
	public void rewriteGTEX() throws IOException {
		String listfile = "/Data/eQTLs/GTEx/v6peQTLs.txt";
		String listfileout = "/Data/eQTLs/GTEx/v6peQTLsTab.txt";
		TextFile tf = new TextFile(listfile, TextFile.R);
		TextFile tfout1 = new TextFile(listfileout, TextFile.W);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			
			String file = elems[0];
			
			String outfile = elems[0] + ".tab.gz";
			TextFile tfin = new TextFile(file, TextFile.R);
			TextFile tfout = new TextFile(outfile, TextFile.W);
			tfin.readLine();
			String[] elems2 = tfin.readLineElems(TextFile.tab);
			int ctr = 0;
			while (elems2 != null) {
// out:
				// Chr1    168921757       rs12095636      ENSG00000000460.11      1.047e-02
				
				//in
				// variant_id      gene_id tss_distance    pval_nominal    slope   slope_se        slope_fpkm      slope_fpkm_se   pval_nominal_threshold  min_pval_nominal        pval_beta
				// 1_662622_G_A_b37
				
				if (elems2.length < 4) {
					System.err.println("Unexpected end of line? " + elems2.length + " elems found...");
					System.out.println(Strings.concat(elems2, Strings.tab));
					
					System.exit(-1);
				}
				String var = elems2[0];
				String[] varelems = var.split("_");
				String outln = varelems[0]
						+ "\t" + varelems[1]
						+ "\t" + var
						+ "\t" + elems2[1]
						+ "\t" + elems2[3];
				tfout.writeln(outln);
				ctr++;
				elems2 = tfin.readLineElems(TextFile.tab);
			}
			tfin.close();
			tfout.close();
			tfout1.writeln(outfile + "\t" + elems[1]);
			System.out.println(ctr + "\tin\t" + outfile);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		tfout1.close();
	}
	
	
	private void bedOverlap(String bedregions,
							String[] bedfiles,
							String[] bedfilenames,
							String[] assocFiles,
							String[] datasetnames,
							String genenamefile,
							String outoverlapfile,
							String outplot,
							double maxPosteriorCredibleSet,
							int maxNrVariantsInCredibleSet,
							String geneAnnotation,
							boolean produceplot) throws IOException, DocumentException {
		
		// load regions
		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bedregions);
		
		HashMap<String, String> locusToGene = loadLocusToGene(genenamefile);
		
		// load annotations
		AnnotationData[] annotationdata = new AnnotationData[bedfiles.length];
		int totalGroups = 0;
		for (int i = 0; i < annotationdata.length; i++) {
			System.out.println("Loading: " + bedfiles[i]);
			annotationdata[i] = new AnnotationData(bedfiles[i], regions);
			totalGroups += annotationdata[i].getUniqueGroups().size();
		}
		
		//
		// region variant1 pval1 posterior1 variant2 pval2 posterior2 variant3 pval3 posterior3
		AssociationFile f = new AssociationFile();
		
		
		// [dataset][region][variants]
		AssociationResult[][][] crediblesets = new AssociationResult[assocFiles.length][regions.size()][];
		
		AssociationResult[][][] data = new AssociationResult[assocFiles.length][regions.size()][];
		
		ApproximateBayesPosterior abp = new ApproximateBayesPosterior();
		ArrayList<Feature> regionsWithCredibleSets = new ArrayList<>();
		
		
		ArrayList<ArrayList<AssociationResult>> associationResults = new ArrayList<>();
		for (int i = 0; i < assocFiles.length; i++) {
			associationResults.add(f.read(assocFiles[i]));
		}
		
		for (int d = 0; d < regions.size(); d++) {
			boolean hasSet = false;
			Feature region = regions.get(d);
			for (int i = 0; i < assocFiles.length; i++) {
				ArrayList<AssociationResult> allDatasetData = filterAssocResults(associationResults.get(i), region);
				
				data[i][d] = allDatasetData.toArray(new AssociationResult[0]);
				ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(allDatasetData, maxPosteriorCredibleSet);
				crediblesets[i][d] = credibleSet.toArray(new AssociationResult[0]);
				if (credibleSet.size() <= maxNrVariantsInCredibleSet) {
					hasSet = true;
				}
			}
			if (hasSet) {
				regionsWithCredibleSets.add(regions.get(d));
			}
		}
		
		HashSet<Feature> regionsWithCredibleSetsHash = new HashSet<Feature>();
		regionsWithCredibleSetsHash.addAll(regionsWithCredibleSets);
		
		int len = maxNrVariantsInCredibleSet;
		
		
		// initiate output..
		TextFile out = new TextFile(outoverlapfile, TextFile.W);
		String header2 = "\t\t";
		
		String header1 = "region\tgene";
		for (int i = 0; i < data.length; i++) {
			header2 += datasetnames[i]
					+ "\t\t\t";
			header1 += "\tNrVariantsInCredibleSet" +
					"\tVariants" +
					"\tPosterior";
			for (int e = 0; e < annotationdata.length; e++) {
				
				// get groups
				ArrayList<String> annotationGroups = annotationdata[e].getUniqueGroups();
				for (int s = 0; s < annotationGroups.size(); s++) {
					header1 += "\t" + bedfilenames[e] + "-" + annotationGroups.get(s) + "(" + annotationdata[e].getNrAnnotationsInGroup(s) + ")";
					header2 += "\t";
				}
			}
		}
		out.writeln(header2);
		out.writeln(header1);
		
		
		// determine overlap
		int nrDatasets = data.length;
		int ctr = 0;
		DecimalFormat format = new DecimalFormat("#.###");
		for (int regionId = 0; regionId < regions.size(); regionId++) {
			Feature region = regions.get(regionId);
			if (regionsWithCredibleSetsHash.contains(region)) {
				ctr++;
				ArrayList<ArrayList<AssociationResult>> resultsPerDs = new ArrayList<>();
				for (int datasetId = 0; datasetId < nrDatasets; datasetId++) {
					ArrayList<AssociationResult> topResults = getTopVariants(data, datasetId, regionId, len);
					resultsPerDs.add(topResults);
				}
				
				double[] regionsums = new double[data.length];
				for (int snpId = 0; snpId < len; snpId++) {
					String ln = "";
					boolean allSNPsPrinted = true;
					for (int datasetId = 0; datasetId < data.length; datasetId++) {
						if (regionsums[datasetId] < maxPosteriorCredibleSet) {
							allSNPsPrinted = false;
						}
					}
					
					if (!allSNPsPrinted) {
						if (snpId == 0) {
							ln = region.toString() + "\t" + locusToGene.get(region.toString());
							
							
							for (int datasetId = 0; datasetId < data.length; datasetId++) {
								AssociationResult r = resultsPerDs.get(datasetId).get(snpId); //data[datasetId][regionId][snpId];
								ln += "\t" + crediblesets[datasetId][regionId].length
										+ "\t" + r.getSnp().toString()
										+ "\t" + r.getPosterior();
								
								String lnblock = "";
								for (int an = 0; an < bedfilenames.length; an++) {
									int nrgroups = annotationdata[an].getUniqueGroups().size();
									for (int group = 0; group < nrgroups; group++) {
										int overlapdata = annotationdata[an].countOverlappingAnnotations(r.getSnp(), group);
										int nrids = annotationdata[an].getNrAnnotationsInGroup(group);
										double perc = (double) overlapdata / nrids;
										lnblock += "\t" + format.format(perc);
									}
								}
								
								ln += lnblock;
								regionsums[datasetId] += r.getPosterior();
							}
						} else {
							ln = "\t";
							for (int datasetId = 0; datasetId < data.length; datasetId++) {
								AssociationResult r = resultsPerDs.get(datasetId).get(snpId); //data[datasetId][regionId][snpId];
								
								String lnblock = "";
								for (int an = 0; an < bedfilenames.length; an++) {
									int nrgroups = annotationdata[an].getUniqueGroups().size();
									for (int group = 0; group < nrgroups; group++) {
										int overlapdata = annotationdata[an].countOverlappingAnnotations(r.getSnp(), group);
										int nrids = annotationdata[an].getNrAnnotationsInGroup(group);
										double perc = (double) overlapdata / nrids;
										lnblock += "\t" + format.format(perc);
									}
								}
								
								if (regionsums[datasetId] < maxPosteriorCredibleSet) {
									ln += "\t"
											+ "\t" + r.getSnp().toString()
											+ "\t" + r.getPosterior();
									ln += lnblock;
								} else {
									ln += "\t"
											+ "\t"
											+ "\t";
									for (int e = 0; e < totalGroups; e++) {
										ln += "\t";
									}
								}
								regionsums[datasetId] += r.getPosterior();
							}
							
						}
						out.writeln(ln);
					}
					
				}
			}
			
			
		}
		
		out.close();
		
		
		//
		System.out.println(ctr + " regions processed.. ");
		int qctr = 0;
		if (produceplot) {
			
			Annotation annotation = null;
			if (geneAnnotation.endsWith(".gtf.gz") || geneAnnotation.endsWith(".gtf")) {
				annotation = new GTFAnnotation(geneAnnotation);
			} else {
				annotation = new EnsemblStructures(geneAnnotation);
			}
			
			for (int regionId = 0; regionId < regions.size(); regionId++) {
				Feature region = regions.get(regionId);
				if (regionsWithCredibleSetsHash.contains(region)) {
					int bp = region.getStop() - region.getStart();
					int bpperbin = 10;
					int nrBins = bp / bpperbin;
					
					System.out.println(nrBins);
					
					
					double[][][] overlap = new double[bedfilenames.length][][]; // annotation, group, bp
					for (int an = 0; an < bedfilenames.length; an++) {
						int nrgroups = annotationdata[an].getUniqueGroups().size();
						double[][] tmpoverlap = new double[nrgroups][nrBins];
						overlap[an] = tmpoverlap;
					}
					
					String[][] groupnames = new String[bedfilenames.length][];
					for (int an = 0; an < bedfilenames.length; an++) {
						int nrgroups = annotationdata[an].getUniqueGroups().size();
						groupnames[an] = new String[nrgroups];
						for (int group = 0; group < nrgroups; group++) {
							groupnames[an][group] = annotationdata[an].getUniqueGroups().get(group) + "-" + bedfilenames[an];
						}
					}
					
					
					int binno = 0;
					for (int s = region.getStart(); s < region.getStop(); s += bpperbin) {
						Feature feat = new Feature();
						feat.setChromosome(region.getChromosome());
						feat.setStart(s);
						feat.setStop(s + bpperbin);
						
						for (int an = 0; an < bedfilenames.length; an++) {
							int nrgroups = annotationdata[an].getUniqueGroups().size();
							for (int group = 0; group < nrgroups; group++) {
								if (binno < nrBins) {
									overlap[an][group][binno] = annotationdata[an].countOverlappingAnnotations(feat, group);
//								System.out.println(overlap[an][group][binno]);
								}
							}
						}
						binno++;
					}
//				System.exit(-1);

//				// plot output
//
//
					for (int ds = 0; ds < datasetnames.length; ds++) {
						Grid grid = new Grid(1000, 100, 3, 1, 100, 100);
						AnnotationTrackPanel p = new AnnotationTrackPanel(1, 1);
						AssociationPanel assocP = new AssociationPanel(1, 1);
						
						
						GenePanel gp = new GenePanel(1, 1);
						ArrayList<Gene> genes = new ArrayList<Gene>();
						for (Gene g : annotation.getGenes()) {
							if (g.overlaps(region)) {
								genes.add(g);
							}
						}
						gp.setData(region, genes);
						grid.addPanel(gp);
						
						ArrayList<AssociationResult> results = associationResults.get(ds);
						ArrayList<Pair<Integer, Double>> pvalues = new ArrayList<>();
						
						for (int i = 0; i < results.size(); i++) {
							AssociationResult r = results.get(i);
							if (r.getSnp().overlaps(region)) {
								pvalues.add(new Pair<Integer, Double>(r.getSnp().getStart(), r.getLog10Pval()));
							}
						}
						
						boolean[] mark = new boolean[pvalues.size()];
						for (int i = 0; i < mark.length; i++) {
							int pos = pvalues.get(i).getLeft();
							for (int j = 0; j < crediblesets[ds][regionId].length; j++) {
								if (pos == crediblesets[ds][regionId][j].getSnp().getStart()) {
									mark[j] = true;
								}
							}
						}
						
						assocP.setDataSingleDs(region, null, pvalues, datasetnames[ds]);
						assocP.setMarkDifferentShape(mark);
						grid.addPanel(assocP);
						
						int[] highlight = new int[crediblesets[ds][regionId].length];
						for (int i = 0; i < highlight.length; i++) {
							highlight[i] = crediblesets[ds][regionId][i].getSnp().getStart();
						}
						
						p.setData(overlap, highlight, region, groupnames);
//						p.setMarginX(10);
//						p.setMarginY(10);
						grid.addPanel(p);
						
						grid.draw(outplot + datasetnames[ds] + "-" + region.toString() + ".pdf");
					}

//
//				qctr++;
//				System.out.println(ctr + " regions processed.. ");
				
				}
			}
		}
		
		
	}
	
	
	private void rewriteEQTLFile(String filename, String dbsnpVCF, String out) throws IOException {
		
		TextFile tf = new TextFile(filename, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, VCFVariant> rsidToVCFVariant = new HashMap<String, VCFVariant>();
		while (elems != null) {
			rsidToVCFVariant.put(elems[0], null);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println(rsidToVCFVariant.size() + " variants to check for rs id positions");
		
		TextFile tf2 = new TextFile(dbsnpVCF, TextFile.R);
		String vcfln = tf2.readLine();
		int ctr = 0;
		int found = 0;
		while (vcfln != null) {
			if (!vcfln.startsWith("#")) {
				VCFVariant var = new VCFVariant(vcfln, VCFVariant.PARSE.HEADER);
				if (rsidToVCFVariant.containsKey(var.getId())) {
					rsidToVCFVariant.put(var.getId(), var);
					found++;
				}
			}
			ctr++;
			if (ctr % 100000 == 0) {
				System.out.println(ctr + " lines read from: dbsnp." + found + " found out of " + rsidToVCFVariant.size());
			}
			vcfln = tf2.readLine();
		}
		tf2.close();
		
		// SNP	PROBESET_ID	GENE	CHR	BETA(META)	SE(META)	P-VALUE
		tf.open();
		tf.readLine();
		ArrayList<EQTL> allEQTLs = new ArrayList<>();
		elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[0];
			VCFVariant var = rsidToVCFVariant.get(snp);
			if (var != null) {
				EQTL eqtl = new EQTL();
				SNPFeature feature = var.asSNPFeature();
				eqtl.setSnp(feature);
				eqtl.setGenename(elems[2]);
				eqtl.setName(elems[2]);
				eqtl.setPval(Double.parseDouble(elems[6]));
				allEQTLs.add(eqtl);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		TextFile outf = new TextFile(out, TextFile.W);
		
		for (EQTL e : allEQTLs) {
			if (e.getSnp().getChromosome() != null) {
				String ln = e.getSnp().getChromosome().toString()
						+ "\t" + e.getSnp().getStart()
						+ "\t" + e.getSnp().getName()
						+ "\t" + e.getGenename()
						+ "\t" + e.getPval();
				outf.writeln(ln);
			}
		}
		
		outf.close();
	}
	
	private HashMap<String, String> loadLocusToGene(String genenamefile) throws IOException {
		HashMap<String, String> locusToGene = new HashMap<String, String>();
		TextFile genefiletf = new TextFile(genenamefile, TextFile.R);
		String[] genelems = genefiletf.readLineElems(TextFile.tab);
		while (genelems != null) {
			if (genelems.length > 1) {
				locusToGene.put(genelems[0], genelems[1]);
			} else {
				locusToGene.put(genelems[0], "");
			}
			genelems = genefiletf.readLineElems(TextFile.tab);
		}
		genefiletf.close();
		return locusToGene;
	}
	
	
	private boolean determineOverlap(Feature[] allannotation, SNPFeature snp) {
		for (Feature f : allannotation) {
			if (f.overlaps(snp)) {
				return true;
			}
		}
		return false;
	}
	
	private Pair<Feature[][][][], String[][]> loadAnnotations(String[] bedfilenames, ArrayList<Feature> regions) throws IOException {
		
		// [region][dataset][celltype][annotations]
		Feature[][][][] output = new Feature[regions.size()][bedfilenames.length][][];
		String[][] filenames = new String[bedfilenames.length][];
		for (int b = 0; b < bedfilenames.length; b++) {
			TextFile f = new TextFile(bedfilenames[b], TextFile.R);
			
			String[] elems = f.readLineElems(TextFile.tab);
			ArrayList<String> files = new ArrayList<>();
			ArrayList<String> bednames = new ArrayList<>();
			while (elems != null) {
				String name = elems[0];
				String file = elems[1];
				bednames.add(name);
				files.add(file);
				elems = f.readLineElems(TextFile.tab);
			}
			f.close();
			
			filenames[b] = bednames.toArray(new String[0]);
			System.out.println(filenames[b].length + " files in " + bedfilenames[b]);
			
			for (int p = 0; p < files.size(); p++) {
				String file = files.get(p);
				
				BedFileReader reader = new BedFileReader();
				ArrayList<Feature> feats = reader.readAsList(file);
				
				for (int r = 0; r < regions.size(); r++) {
					if (output[r][b] == null) {
						output[r][b] = new Feature[filenames[b].length][];
					}
					ArrayList<Feature> selected = feats.stream().filter(feat -> feat.overlaps(regions)).collect(Collectors.toCollection(ArrayList::new));
					output[r][b][p] = selected.toArray(new Feature[0]);
				}
				System.out.println("Done processing " + p + " out of " + filenames[b].length);
			}
		}
		
		return new Pair<>(output, filenames);
	}
	
	
	private void eQTLOverlap(String bedregions,
							 String[] eqtlfiles,
							 String[] eqtlfilenames,
	
							 String tabixPrefix,
							 String tabixFilter,
							 String[] assocFiles,
							 String[] datasetNames,
							 String genenamefile,
							 String outfile, double maxPosteriorCredibleSet, int maxNrVariantsInCredibleSet, String annot) throws
			IOException {
//		GTFAnnotation annotation = new GTFAnnotation(annot);
//		TreeSet<Gene> genes = annotation.getGeneTree();
		
		HashMap<String, String> locusToGene = loadLocusToGene(genenamefile);
		
		
		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bedregions);
		
		// region variant1 pval1 posterior1 variant2 pval2 posterior2 variant3 pval3 posterior3
		AssociationFile f = new AssociationFile();
		
		
		// [dataset][region][variants]
		AssociationResult[][][] crediblesets = new AssociationResult[assocFiles.length][regions.size()][];
		
		AssociationResult[][][] data = new AssociationResult[assocFiles.length][regions.size()][];
		
		ApproximateBayesPosterior abp = new ApproximateBayesPosterior();
		ArrayList<Feature> regionsWithCredibleSets = new ArrayList<>();
		
		
		ArrayList<ArrayList<AssociationResult>> associationResults = new ArrayList<>();
		for (int i = 0; i < assocFiles.length; i++) {
			associationResults.add(f.read(assocFiles[i]));
		}
		
		for (int d = 0; d < regions.size(); d++) {
			boolean hasSet = false;
			Feature region = regions.get(d);
			for (int i = 0; i < assocFiles.length; i++) {
				ArrayList<AssociationResult> allDatasetData = filterAssocResults(associationResults.get(i), region);
				
				data[i][d] = allDatasetData.toArray(new AssociationResult[0]);
				ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(allDatasetData, maxPosteriorCredibleSet);
				crediblesets[i][d] = credibleSet.toArray(new AssociationResult[0]);
				if (credibleSet.size() <= maxNrVariantsInCredibleSet) {
					hasSet = true;
				}
			}
			if (hasSet) {
				regionsWithCredibleSets.add(regions.get(d));
			}
		}
		
		HashSet<Feature> regionsWithCredibleSetsHash = new HashSet<Feature>();
		regionsWithCredibleSetsHash.addAll(regionsWithCredibleSets);
		
		int len = maxNrVariantsInCredibleSet;
		TextFile out = new TextFile(outfile, TextFile.W);
		TextFile out2 = new TextFile(outfile + "-pvalsonly.txt", TextFile.W);
		String header2 = "\t\t";
		
		String header1 = "region\tgene";
		String header3 = "region\tgene";
		for (int i = 0; i < data.length; i++) {
			header2 += datasetNames[i]
					+ "\t\t\t\t\t";
			header1 += "\tNrVariantsInCredibleSet" +
					"\tVariants" +
					"\tPosterior";
			header3 += "\tNrVariantsInCredibleSet" +
					"\tVariants" +
					"\tPosterior";
			for (int e = 0; e < eqtlfilenames.length; e++) {
				header1 += "\teSNP-" + eqtlfilenames[e]
						+ "\teGene-" + eqtlfilenames[e]
						+ "\tP(eQTL)-" + eqtlfilenames[e]
						+ "\tDistance-" + eqtlfilenames[e]
						+ "\tLD-" + eqtlfilenames[e];
				header3 += "\tP(eQTL)-" + eqtlfilenames[e];
				header2 += "\t\t\t\t";
			}
		}
		out.writeln(header2);
		out.writeln(header1);
		out2.writeln(header3);
		
		// load all top eQTLs per gene for all regions (+- 1mb)
		EQTL[][][] eqtls = loadEQTLs(eqtlfiles, regions); // [eqtldataset][region][eqtl]
		
		int nrDatasets = data.length;
		int ctr = 0;
		for (int regionId = 0; regionId < regions.size(); regionId++) {
			Feature region = regions.get(regionId);
			if (regionsWithCredibleSetsHash.contains(region)) {
				
				ctr++;
				
				// get all VCFVariants in region
				String tabixFile = tabixPrefix.replaceAll("CHR", "" + region.getChromosome().getNumber());
				
				VCFTabix tabix = new VCFTabix(tabixFile);
				boolean[] filter = tabix.getSampleFilter(tabixFilter);
				Feature eqtlregion = new Feature(region);
				eqtlregion.setStart(eqtlregion.getStart() - 10000);
				eqtlregion.setStop(eqtlregion.getStop() + 10000);
				ArrayList<VCFVariant> all1kgvariants = tabix.getAllVariants(region, filter);
				System.out.println(all1kgvariants.size() + " variants in LD reference for region: " + region.toString());
				
				ArrayList<ArrayList<AssociationResult>> resultsPerDs = new ArrayList<>();
				
				for (int datasetId = 0; datasetId < nrDatasets; datasetId++) {
					ArrayList<AssociationResult> topResults = getTopVariants(data, datasetId, regionId, len);
					resultsPerDs.add(topResults);
				}
				
				double[] regionsums = new double[data.length];
				for (int snpId = 0; snpId < len; snpId++) {
					String ln = "";
					String ln2 = "";
					boolean allSNPsPrinted = true;
					for (int datasetId = 0; datasetId < data.length; datasetId++) {
						if (regionsums[datasetId] < maxPosteriorCredibleSet) {
							allSNPsPrinted = false;
						}
					}
					
					if (!allSNPsPrinted) {
						if (snpId == 0) {
							ln = region.toString() + "\t" + locusToGene.get(region.toString());
							ln2 = region.toString() + "\t" + locusToGene.get(region.toString());
							for (int datasetId = 0; datasetId < data.length; datasetId++) {
								AssociationResult r = resultsPerDs.get(datasetId).get(snpId); //data[datasetId][regionId][snpId];
								ln += "\t" + crediblesets[datasetId][regionId].length
										+ "\t" + r.getSnp().toString()
										+ "\t" + r.getPosterior();
								ln2 += "\t" + crediblesets[datasetId][regionId].length
										+ "\t" + r.getSnp().toString()
										+ "\t" + r.getPosterior();
								
								
								String lnblock = eqtllineblock(all1kgvariants,
										r,
										eqtlfilenames,
										eqtls,
										regionId);
								
								ln += lnblock;
								String[] lnblockelems = lnblock.split("\t");
								String p = lnblockelems[4];
								ln2 += "\t" + p;
								regionsums[datasetId] += r.getPosterior();
							}
						} else {
							ln = "\t";
							ln2 = "\t";
							for (int datasetId = 0; datasetId < data.length; datasetId++) {
								AssociationResult r = resultsPerDs.get(datasetId).get(snpId); //data[datasetId][regionId][snpId];
								
								
								String lnblock = eqtllineblock(all1kgvariants,
										r,
										eqtlfilenames,
										eqtls,
										regionId);
								
								if (regionsums[datasetId] < maxPosteriorCredibleSet) {
									ln += "\t"
											+ "\t" + r.getSnp().toString()
											+ "\t" + r.getPosterior();
									ln2 += "\t"
											+ "\t" + r.getSnp().toString()
											+ "\t" + r.getPosterior();
									ln += lnblock;
									String[] lnblockelems = lnblock.split("\t");
									String p = lnblockelems[4];
									ln2 += "\t" + p;
								} else {
									ln += "\t"
											+ "\t"
											+ "\t";
									ln2 += "\t"
											+ "\t"
											+ "\t";
									for (int e = 0; e < eqtlfilenames.length; e++) {
										ln += "\t\t\t\t\t";
										ln2 += "\t";
									}
								}
								
								
								regionsums[datasetId] += r.getPosterior();
							}
							
						}
						out.writeln(ln);
						out2.writeln(ln2);
					}
					
				}
				System.out.println(ctr + "/" + regionsWithCredibleSets.size() + " regions processed....");
			}
			
			
		}
		
		out.close();
		out2.close();
	}
	
	private void eQTLOverlapGTEX(String bedregions,
								 String[] eqtlfiles,
								 String[] eqtlfilenames,
								 String tabixPrefix,
								 String tabixFilter,
								 String[] assocFiles,
								 String[] datasetNames,
								 String genenamefile,
								 String outfile, double maxPosteriorCredibleSet, int maxNrVariantsInCredibleSet, String annot) throws
			IOException {
//		GTFAnnotation annotation = new GTFAnnotation(annot);
//		TreeSet<Gene> genes = annotation.getGeneTree();
		
		HashMap<String, String> locusToGene = loadLocusToGene(genenamefile);
		
		
		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bedregions);
		
		// region variant1 pval1 posterior1 variant2 pval2 posterior2 variant3 pval3 posterior3
		AssociationFile f = new AssociationFile();
		
		
		// [dataset][region][variants]
		AssociationResult[][][] crediblesets = new AssociationResult[assocFiles.length][regions.size()][];
		
		AssociationResult[][][] data = new AssociationResult[assocFiles.length][regions.size()][];
		
		ApproximateBayesPosterior abp = new ApproximateBayesPosterior();
		ArrayList<Feature> regionsWithCredibleSets = new ArrayList<>();
		
		
		ArrayList<ArrayList<AssociationResult>> associationResults = new ArrayList<>();
		for (int i = 0; i < assocFiles.length; i++) {
			associationResults.add(f.read(assocFiles[i]));
		}
		
		for (int d = 0; d < regions.size(); d++) {
			boolean hasSet = false;
			Feature region = regions.get(d);
			for (int i = 0; i < assocFiles.length; i++) {
				ArrayList<AssociationResult> allDatasetData = filterAssocResults(associationResults.get(i), region);
				
				data[i][d] = allDatasetData.toArray(new AssociationResult[0]);
				ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(allDatasetData, maxPosteriorCredibleSet);
				crediblesets[i][d] = credibleSet.toArray(new AssociationResult[0]);
				if (credibleSet.size() <= maxNrVariantsInCredibleSet) {
					hasSet = true;
				}
			}
			if (hasSet) {
				regionsWithCredibleSets.add(regions.get(d));
			}
		}
		
		HashSet<Feature> regionsWithCredibleSetsHash = new HashSet<Feature>();
		regionsWithCredibleSetsHash.addAll(regionsWithCredibleSets);
		
		int len = maxNrVariantsInCredibleSet;
		TextFile out = new TextFile(outfile, TextFile.W);
		String header2 = "\t\t";
		
		String header1 = "region\tgene";
		
		for (int i = 0; i < data.length; i++) {
			header2 += datasetNames[i] + "" +
					"\t" +
					"\t" +
					"\t" +
					"\t" +
					"\t";
			header1 += "\tNrVariantsInCredibleSet" +
					"\tVariants" +
					"\tPosterior" +
					"\tNrTissues" +
					"\tTissues" +
					"\tGenes" +
					"\tSNPs" +
					"\tP(eQTL)" +
					"\tLD";
		}
		out.writeln(header2);
		out.writeln(header1);
		
		
		ArrayList<String> eQTLFileListArr = new ArrayList<>();
		ArrayList<String> eQTLFileListNameArr = new ArrayList<>();
		TextFile tfq = new TextFile(eqtlfiles[0], TextFile.R);
		String[] tfqel = tfq.readLineElems(TextFile.tab);
		while (tfqel != null) {
			eQTLFileListArr.add(tfqel[0]);
			eQTLFileListNameArr.add(tfqel[1]);
			tfqel = tfq.readLineElems(TextFile.tab);
		}
		tfq.close();
		
		eqtlfiles = eQTLFileListArr.toArray(new String[0]);
		eqtlfilenames = eQTLFileListNameArr.toArray(new String[0]);
		// load all top eQTLs per gene for all regions (+- 1mb)
		EQTL[][][] eqtls = loadEQTLs(eqtlfiles, regions); // [eqtldataset][region][eqtl]
		
		int nrDatasets = data.length;
		int ctr = 0;
		
		for (int regionId = 0; regionId < regions.size(); regionId++) {
//		for (int regionId = 0; regionId < 10; regionId++) {
			Feature region = regions.get(regionId);
			if (regionsWithCredibleSetsHash.contains(region)) {
				
				ctr++;
				
				// get all VCFVariants in region
				String tabixFile = tabixPrefix.replaceAll("CHR", "" + region.getChromosome().getNumber());
				
				VCFTabix tabix = new VCFTabix(tabixFile);
				boolean[] filter = tabix.getSampleFilter(tabixFilter);
				Feature eqtlregion = new Feature(region);
				eqtlregion.setStart(eqtlregion.getStart() - 10000);
				eqtlregion.setStop(eqtlregion.getStop() + 10000);
				ArrayList<VCFVariant> all1kgvariants = tabix.getAllVariants(region, filter);
				System.out.println(all1kgvariants.size() + " variants in LD reference for region: " + region.toString());
				
				ArrayList<ArrayList<AssociationResult>> resultsPerDs = new ArrayList<>();
				
				for (int datasetId = 0; datasetId < nrDatasets; datasetId++) {
					ArrayList<AssociationResult> topResults = getTopVariants(data, datasetId, regionId, len);
					resultsPerDs.add(topResults);
				}
				
				double[] regionsums = new double[data.length];
				for (int snpId = 0; snpId < len; snpId++) {
					String ln = "";
					
					boolean allSNPsPrinted = true;
					for (int datasetId = 0; datasetId < data.length; datasetId++) {
						if (regionsums[datasetId] < maxPosteriorCredibleSet) {
							allSNPsPrinted = false;
						}
					}
					
					if (!allSNPsPrinted) {
						if (snpId == 0) {
							ln = region.toString() + "\t" + locusToGene.get(region.toString());
							
							for (int datasetId = 0; datasetId < data.length; datasetId++) {
								AssociationResult r = resultsPerDs.get(datasetId).get(snpId); //data[datasetId][regionId][snpId];
								ln += "\t" + crediblesets[datasetId][regionId].length
										+ "\t" + r.getSnp().toString()
										+ "\t" + r.getPosterior();
								
								String lnblock = eqtllineblockGTEX(all1kgvariants,
										r,
										eqtlfilenames,
										eqtls,
										regionId);
								
								ln += lnblock;
								
								
								regionsums[datasetId] += r.getPosterior();
							}
						} else {
							ln = "\t";
							
							for (int datasetId = 0; datasetId < data.length; datasetId++) {
								AssociationResult r = resultsPerDs.get(datasetId).get(snpId); //data[datasetId][regionId][snpId];
								
								
								String lnblock = eqtllineblockGTEX(all1kgvariants,
										r,
										eqtlfilenames,
										eqtls,
										regionId);
								
								if (regionsums[datasetId] < maxPosteriorCredibleSet) {
									ln += "\t"
											+ "\t" + r.getSnp().toString()
											+ "\t" + r.getPosterior();
									ln += lnblock;
								} else {
									ln += "\t"
											+ "\t"
											+ "\t"
											+ "\t"
											+ "\t"
											+ "\t"
											+ "\t"
											+ "\t"
											+ "\t";
								}
								
								
								regionsums[datasetId] += r.getPosterior();
							}
							
						}
						out.writeln(ln);
						
					}
					
				}
				System.out.println(ctr + "/" + regionsWithCredibleSets.size() + " regions processed....");
			}
			
			
		}
		
		out.close();
		
	}
	
	private String eqtllineblock(ArrayList<VCFVariant> all1kgvariants,
								 AssociationResult r,
								 String[] eqtlfilenames,
								 EQTL[][][] eqtls,
								 int regionId) {
		String ln = "";
		VCFVariant finemapvariant = getVariant(all1kgvariants, r.getSnp());
		
		DetermineLD ldcal = new DetermineLD();
		for (int e = 0; e < eqtlfilenames.length; e++) {
			// get linked eQTLs for this dataset
			// get top linked eQTL for this eQTL dataset
			// print top effect
			
			EQTL[] eqtlsInRegion = eqtls[e][regionId];
			
			double maxLD = -1;
			double minpval = 2;
			EQTL maxEQTL = null;
			int maxdistance = Integer.MAX_VALUE;
//			System.out.println();
			for (EQTL eqtl : eqtlsInRegion) {
				if (eqtl.getSnp().getStart() == r.getSnp().getStart()) {
					// found it
					
					maxLD = 1;
					maxEQTL = eqtl;
					maxdistance = 0;
					double pval = eqtl.getPval();
//					System.out.println("select\t" + r.getSnp().toString() + "\t" + eqtl.getSnp().toString() + "\t" + 0 + "\t" + 1 + "\t" + pval);
					break;
				} else {
					if (finemapvariant != null) {
						VCFVariant eqtlvariant = getVariant(all1kgvariants, eqtl.getSnp());
						
						
						if (finemapvariant != null && eqtlvariant != null) {
							Pair<Double, Double> ld = ldcal.getLD(finemapvariant, eqtlvariant);
							double rsq = ld.getRight();
							int distance = Math.abs(finemapvariant.getPos() - eqtlvariant.getPos());
							
							double pval = eqtl.getPval();
							if (rsq > maxLD) {
								maxLD = rsq;
								maxEQTL = eqtl;
								maxdistance = distance;
								minpval = pval;
//								System.out.println("select\t" + r.getSnp().toString() + "\t" + eqtl.getSnp().toString() + "\t" + distance + "\t" + rsq + "\t" + pval);
							} else if (rsq == maxLD) {
								if (pval < minpval) {
									maxLD = rsq;
									maxEQTL = eqtl;
									maxdistance = distance;
									minpval = pval;
//									System.out.println("select\t" + r.getSnp().toString() + "\t" + eqtl.getSnp().toString() + "\t" + distance + "\t" + rsq + "\t" + pval);
								}
							}
						}
					}
				}
			}
			
			if (maxEQTL == null || maxLD < 0.8) {
				ln += "\t-"
						+ "\t-"
						+ "\t-"
						+ "\t-"
						+ "\t-";
//				System.out.println("No matching eQTL");
			} else {
				// esnp egene pval ld
				ln += "\t" + maxEQTL.getSnp().toString()
						+ "\t" + maxEQTL.getGenename()
						+ "\t" + maxEQTL.getPval()
						+ "\t" + maxdistance
						+ "\t" + maxLD;
			}
		}
		return ln;
	}
	
	private String eqtllineblockGTEX(ArrayList<VCFVariant> all1kgvariants,
									 AssociationResult r,
									 String[] eqtlfilenames,
									 EQTL[][][] eqtls,
									 int regionId) {
		
		VCFVariant finemapvariant = getVariant(all1kgvariants, r.getSnp());
		
		DetermineLD ldcal = new DetermineLD();
		
		// out:
		/*
		header1 += "\tNrVariantsInCredibleSet" +
					"\tVariants" +
					"\tPosterior" +
					"\tTissues" +
					"\tGenes" +
					"\tSNPs" +
					"\tP(eQTL)" +
					"\tLD";
		 */
		
		
		ArrayList<String> tissues = new ArrayList<>();
		ArrayList<String> genes = new ArrayList<>();
		ArrayList<String> SNPs = new ArrayList<>();
		ArrayList<String> peqtl = new ArrayList<>();
		ArrayList<String> LD = new ArrayList<>();
		
		
		for (int e = 0; e < eqtlfilenames.length; e++) {
			// get linked eQTLs for this dataset
			// get top linked eQTL for this eQTL dataset
			// print top effect
			
			EQTL[] eqtlsInRegion = eqtls[e][regionId];
			
			double maxLD = -1;
			double minpval = 2;
			EQTL maxEQTL = null;
			int maxdistance = Integer.MAX_VALUE;
//			System.out.println();
			for (EQTL eqtl : eqtlsInRegion) {
				if (eqtl.getSnp().getStart() == r.getSnp().getStart()) {
					// found it
					
					maxLD = 1;
					maxEQTL = eqtl;
					maxdistance = 0;
					double pval = eqtl.getPval();
//					System.out.println("select\t" + r.getSnp().toString() + "\t" + eqtl.getSnp().toString() + "\t" + 0 + "\t" + 1 + "\t" + pval);
					break;
				} else {
					if (finemapvariant != null) {
						VCFVariant eqtlvariant = getVariant(all1kgvariants, eqtl.getSnp());
						
						
						if (finemapvariant != null && eqtlvariant != null) {
							Pair<Double, Double> ld = ldcal.getLD(finemapvariant, eqtlvariant);
							double rsq = ld.getRight();
							int distance = Math.abs(finemapvariant.getPos() - eqtlvariant.getPos());
							
							double pval = eqtl.getPval();
							if (rsq > maxLD) {
								maxLD = rsq;
								maxEQTL = eqtl;
								maxdistance = distance;
								minpval = pval;
//								System.out.println("select\t" + r.getSnp().toString() + "\t" + eqtl.getSnp().toString() + "\t" + distance + "\t" + rsq + "\t" + pval);
							} else if (rsq == maxLD) {
								if (pval < minpval) {
									maxLD = rsq;
									maxEQTL = eqtl;
									maxdistance = distance;
									minpval = pval;
//									System.out.println("select\t" + r.getSnp().toString() + "\t" + eqtl.getSnp().toString() + "\t" + distance + "\t" + rsq + "\t" + pval);
								}
							}
						}
					}
				}
			}
			
			if (maxEQTL == null || maxLD < 0.8) {
			} else {
				// esnp egene pval ld
				tissues.add(eqtlfilenames[e]);
				genes.add(maxEQTL.getGenename());
				SNPs.add(maxEQTL.getSnp().toString());
				peqtl.add("" + maxEQTL.getPval());
				LD.add("" + maxLD);
			}
		}
		
		String ln = "\t" + tissues.size()
				+ "\t" + Strings.concat(tissues, Strings.semicolon)
				+ "\t" + Strings.concat(genes, Strings.semicolon)
				+ "\t" + Strings.concat(SNPs, Strings.semicolon)
				+ "\t" + Strings.concat(peqtl, Strings.semicolon)
				+ "\t" + Strings.concat(LD, Strings.semicolon);
		
		return ln;
	}
	
	private ArrayList<AssociationResult> filterAssocResults(ArrayList<AssociationResult> associationResults, Feature region) {
		ArrayList<AssociationResult> output = new ArrayList<>();
		for (AssociationResult a : associationResults) {
			if (a.getSnp().overlaps(region)) {
				output.add(a);
			}
		}
		return output;
	}
	
	private VCFVariant getVariant(ArrayList<VCFVariant> all1kgvariants, SNPFeature snp) {
		for (VCFVariant v : all1kgvariants) {
			if (v.asSNPFeature().overlaps(snp)) {
				return v;
			}
		}
		return null;
	}
	
	private EQTL[][][] loadEQTLs(String[] eqtlfilenames, ArrayList<Feature> regions) throws IOException {
		EQTL[][][] output = new EQTL[eqtlfilenames.length][regions.size()][];
		ArrayList<Feature> eqtlregions = new ArrayList<>();
		for (int r = 0; r < regions.size(); r++) {
			Feature eqtlregion = new Feature(regions.get(r));
			eqtlregion.setStart(eqtlregion.getStart() - 1000000);
			if (eqtlregion.getStart() < 0) {
				eqtlregion.setStart(0);
			}
			eqtlregion.setStop(eqtlregion.getStop() + 1000000);
			if (eqtlregion.getStop() > eqtlregion.getChromosome().getLength()) {
				eqtlregion.setStop(eqtlregion.getChromosome().getLength());
			}
			eqtlregions.add(eqtlregion);
		}
		
		
		for (int d = 0; d < eqtlfilenames.length; d++) {
			System.out.println("Parsing: " + eqtlfilenames[d]);
			
			
			ArrayList<EQTL> allEQTLs = new ArrayList<>();
			String filename = eqtlfilenames[d];
			if (filename.endsWith("tab") || filename.endsWith("tab.gz")) {
				
				TextFile tf = new TextFile(eqtlfilenames[d], TextFile.R);
				// Chr3	11604119	rs1000010	ATG7	3.91E-9
				String[] elems = tf.readLineElems(TextFile.tab);
				while (elems != null) {
					
					Chromosome chr = Chromosome.parseChr(elems[0]);
					Integer pos = Integer.parseInt(elems[1]);
					String snpname = elems[2];
					String gene = elems[3];
					double pval = Double.parseDouble(elems[4]);
					SNPFeature feature = new SNPFeature(chr, pos, pos);
					feature.setName(snpname);
					
					EQTL eqtl = new EQTL();
					eqtl.setSnp(feature);
					eqtl.setGenename(gene);
					eqtl.setName(gene);
					eqtl.setPval(pval);
					
					if (eqtloverlap(eqtl, eqtlregions)) {
						allEQTLs.add(eqtl);
					}
					
					elems = tf.readLineElems(TextFile.tab);
				}
				tf.close();
				
				
			} else {
				
				// BIOS eQTL file format
				TextFile tf = new TextFile(eqtlfilenames[d], TextFile.R);
				String[] header = tf.readLineElems(TextFile.tab);
				
				int snpcol = -1;
				int pvalcol = -1;
				int chrcol = -1;
				int poscol = -1;
				int genecol = -1;
				int fdrcol = -1;
				
				for (int q = 0; q < header.length; q++) {
					String h = header[q];
					if (h.equals("PValue")) {
						pvalcol = q;
					}
					if (h.equals("SNPName")) {
						snpcol = q;
					}
					if (h.equals("SNPChr")) {
						chrcol = q;
					}
					if (h.equals("SNPChrPos")) {
						poscol = q;
					}
					if (h.equals("HGNCName")) {
						genecol = q;
					}
					if (h.equals("FDR") && fdrcol == -1) {
						fdrcol = q;
					}
					
				}


				/*
				PValue  SNPName SNPChr  SNPChrPos       ProbeName       ProbeChr        ProbeCenterChrPos       CisTrans        SNPType AlleleAssessed  DatasetsNrSamples       OverallZScore   IncludedDatasetsCorrelationCoefficient  HGNCName        FDR

PValue  SNPName SNPChr  SNPChrPos       ProbeName       ProbeChr        ProbeCenterChrPos       CisTrans        SNPType AlleleAssessed  OverallZScore   DatasetsWhereSNPProbePairIsAvailableAndPassesQC DatasetsZScores DatasetsNrSamples
       IncludedDatasetsMeanProbeExpression     IncludedDatasetsProbeExpressionVariance HGNCName        IncludedDatasetsCorrelationCoefficient  Meta-Beta (SE)  Beta (SE)       FoldChange      FDR     FDR

				 */
				
				String[] elems = tf.readLineElems(TextFile.tab);
				int ctr = 0;
				while (elems != null) {
					
					String snp = elems[snpcol];
					Double pval = Double.parseDouble(elems[pvalcol]);
					Chromosome chr = Chromosome.parseChr(elems[chrcol]);
					Integer pos = Integer.parseInt(elems[poscol]);
					SNPFeature feature = new SNPFeature(chr, pos, pos);
					feature.setName(snp);
					
					EQTL eqtl = new EQTL();
					eqtl.setSnp(feature);
					eqtl.setGenename(elems[genecol]);
					eqtl.setName(elems[genecol]);
					eqtl.setPval(pval);
					
					boolean sig = true;
					if (fdrcol != -1) {
						Double fdr = Double.parseDouble(elems[fdrcol]);
						if (fdr > 0.05) {
							sig = false;
						}
					}
					if (sig && eqtloverlap(eqtl, eqtlregions)) {
						allEQTLs.add(eqtl);
					}
					
					ctr++;
					if (ctr % 100000 == 0) {
						System.out.println(ctr + " lines parsed. " + allEQTLs.size() + " stored");
					}
					elems = tf.readLineElems(TextFile.tab);
				}
				tf.close();
			}
			
			System.out.println(allEQTLs.size() + " eqtls finally loaded from: " + filename);
			for (int r = 0; r < eqtlregions.size(); r++) {
				ArrayList<EQTL> overlapping = new ArrayList<>();
				Feature eqtlregion = eqtlregions.get(r);
				for (int e = 0; e < allEQTLs.size(); e++) {
					EQTL eqtl = allEQTLs.get(e);
					
					if (eqtl.getSnp().overlaps(eqtlregion)) {
						overlapping.add(eqtl);
					}
				}
				output[d][r] = overlapping.toArray(new EQTL[0]);
			}
			System.out.println("done loading eQTLs");
		}
		
		return output;
		
	}
	
	private boolean eqtloverlap(EQTL eqtl, ArrayList<Feature> eqtlregions) {
		for (int r = 0; r < eqtlregions.size(); r++) {
			Feature eqtlregion = eqtlregions.get(r);
			if (eqtl.getSnp().overlaps(eqtlregion)) {
				return true;
			}
		}
		return false;
		
	}
	
	public void mergeCredibleSets(String bedregions,
								  String[] assocFiles,
								  String[] datasetNames,
								  String genenamefile,
								  String outfile,
								  double maxPosteriorCredibleSet,
								  double variantSignificanceThreshold,
								  int maxNrVariantsInCredibleSet,
								  String annot,
								  boolean includeAllLoci) throws IOException {
		
		
		Annotation annotation = null;
		if (annot.endsWith(".gtf.gz") || annot.endsWith(".gtf")) {
			annotation = new GTFAnnotation(annot);
		} else {
			annotation = new EnsemblStructures(annot);
		}
		
		TreeSet<Gene> genes = annotation.getGeneTree();
		
		HashMap<String, String> locusToGene = loadLocusToGene(genenamefile);
		
		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bedregions);
		
		// region variant1 pval1 posterior1 variant2 pval2 posterior2 variant3 pval3 posterior3
		AssociationFile f = new AssociationFile();
		
		
		// [dataset][region][variants]
		AssociationResult[][][] crediblesets = new AssociationResult[assocFiles.length][regions.size()][];
		AssociationResult[][][] data = new AssociationResult[assocFiles.length][regions.size()][];
		
		boolean[][] regionSignificant = new boolean[assocFiles.length][regions.size()];
		
		ApproximateBayesPosterior abp = new ApproximateBayesPosterior();
		ArrayList<Feature> regionsWithCredibleSets = new ArrayList<>();
		{
			ArrayList<ArrayList<AssociationResult>> allResults = new ArrayList<>();
			for (int i = 0; i < assocFiles.length; i++) {
				ArrayList<AssociationResult> allDatasetData = f.read(assocFiles[i]);
				allResults.add(allDatasetData);
			}
			
			for (int d = 0; d < regions.size(); d++) {
				boolean hasSet = false;
				for (int i = 0; i < assocFiles.length; i++) {
					ArrayList<AssociationResult> allDatasetData = filter(allResults.get(i), regions.get(d));
					data[i][d] = allDatasetData.toArray(new AssociationResult[0]);
					
					regionSignificant[i][d] = hasSignificantResult(allDatasetData, variantSignificanceThreshold);
					
					ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(allDatasetData, maxPosteriorCredibleSet);
					crediblesets[i][d] = credibleSet.toArray(new AssociationResult[0]);
					if (includeAllLoci || credibleSet.size() <= maxNrVariantsInCredibleSet) {
						hasSet = true;
					}
				}
				if (hasSet) {
					regionsWithCredibleSets.add(regions.get(d));
				}
			}
		}
		
		HashSet<Feature> regionsWithCredibleSetsHash = new HashSet<Feature>();
		regionsWithCredibleSetsHash.addAll(regionsWithCredibleSets);
		
		int len = maxNrVariantsInCredibleSet;
		TextFile out = new TextFile(outfile, TextFile.W);
		TextFile outg = new TextFile(outfile + "-genes.txt", TextFile.W);
		TextFile outa = new TextFile(outfile + "-coding.txt", TextFile.W);
		TextFile outi = new TextFile(outfile + "-indel.txt", TextFile.W);
		TextFile outr = new TextFile(outfile + "-regions.bed", TextFile.W);
		String header2 = "\t\t";
		
		String header1 = "region\tgene";
		for (int i = 0; i < data.length; i++) {
			header2 += datasetNames[i]
					+ "\t"
					+ "\t"
					+ "\t"
					+ "\t"
					+ "\t"
					+ "\t"
					+ "\t"
					+ "\t"
					+ "\t"
					+ "\t"
					+ "\t"
					+ "\t"
					+ "\t";
			header1 += "\tRegionSignificant" +
					"\tNrVariantsInCredibleSet" +
					"\tSumPosteriorTop" + maxNrVariantsInCredibleSet + "Variants" +
					"\tVariants" +
					"\tAlleles" +
					"\tINFO" +
					"\tAltAFControls" +
					"\tAltAFCases" +
					"\tOR (alt)" +
					"\tPval" +
					"\tPosterior" +
					"\tSignificant" +
					"\tAnnotation";
		}
		
		out.writeln(header2);
		out.writeln(header1);
		for (int regionId = 0; regionId < regions.size(); regionId++) {
			Feature region = regions.get(regionId);
			boolean outputregioninbed = false;
			if (regionsWithCredibleSetsHash.contains(region)) {
				
				// region nrCrediblesetVariants posteriorsumtop5 topvariants alleles or pval posterior
				ArrayList<ArrayList<AssociationResult>> resultsPerDs = new ArrayList<>();
				SNPClass[][][] variantAnnotations = new SNPClass[data.length][][];
				String[][][][] variantGeneAnnotations = new String[data.length][][][];
				for (int i = 0; i < data.length; i++) {
					ArrayList<AssociationResult> topResults = getTopVariants(data, i, regionId, len);
					resultsPerDs.add(topResults);
					
					ArrayList<SNPFeature> variants = new ArrayList<>();
					
					if (topResults == null) {
						System.out.println(assocFiles[i] + " produces null results for region " + regionId + "\t" + region);
						System.exit(-1);
					}
					
					for (int q = 0; q < topResults.size(); q++) {
						variants.add(topResults.get(q).getSnp());
					}
					Pair<SNPClass[][], String[][][]> annotationPair = annotateVariants(variants, region, genes);
					variantAnnotations[i] = annotationPair.getLeft();
					variantGeneAnnotations[i] = annotationPair.getRight(); // [variant][coding/promoter][genename]
					
					
					for (int q = 0; q < variantAnnotations[i].length; q++) {
						if (variantAnnotations[i][q][0] != null && variantAnnotations[i][q][0].equals(SNPClass.EXONIC)) {
							AssociationResult r = topResults.get(q);
							double p = r.getPosterior();
							if (p > 0.20) {
								String outln = datasetNames[i] + "\t" + region.toString() + "\t" + locusToGene.get(region.toString()) + "\t" + p + "\t" + r.getSnp().getName() + "\t" + r.getPval();
								outa.writeln(outln);
							}
						} else if (variantAnnotations[i][q][2] != null && variantAnnotations[i][q][2].equals(SNPClass.INDEL)) {
							AssociationResult r = topResults.get(q);
							double p = r.getPosterior();
							if (p > 0.20) {
								String outln = datasetNames[i] + "\t" + region.toString() + "\t" + locusToGene.get(region.toString()) + "\t" + p + "\t" + r.getSnp().getName() + "\t" + r.getPval();
								outi.writeln(outln);
							}
						}
					}
					
					
				}
				
				
				double[] sumsperregion = new double[datasetNames.length];
				for (int datasetId = 0; datasetId < datasetNames.length; datasetId++) {
					double sum = 0;
					for (int snpId = 0; snpId < len; snpId++) {
						AssociationResult r = resultsPerDs.get(datasetId).get(snpId); //data[datasetId][regionId][snpId];
						
						sum += r.getPosterior();
						
					}
					sumsperregion[datasetId] = sum;
				}
				
				
				outg.writeln(region.toString() + "\t" + locusToGene.get(region.toString()));
				double[] regionsums = new double[data.length];
				
				for (int snpId = 0; snpId < len; snpId++) {
					String ln = "";
					boolean allSNPsPrinted = true;
					for (int datasetId = 0; datasetId < data.length; datasetId++) {
						if (regionsums[datasetId] < maxPosteriorCredibleSet) {
							allSNPsPrinted = false;
						}
					}
					
					if (!allSNPsPrinted) {
						if (snpId == 0) {
							ln = region.toString() + "\t" + locusToGene.get(region.toString());
							for (int datasetId = 0; datasetId < data.length; datasetId++) {
								String annotStr = "";
								if (variantAnnotations[datasetId].length > snpId) {
									SNPClass[] snpAnnotation = variantAnnotations[datasetId][snpId];
									if (snpAnnotation != null) {
										for (int a = 0; a < snpAnnotation.length; a++) {
											String geneStr = null;
											if (a == 0) {
												geneStr = Strings.concat(variantGeneAnnotations[datasetId][snpId][0], Strings.semicolon);
											}
											if (a == 1) {
												geneStr = Strings.concat(variantGeneAnnotations[datasetId][snpId][1], Strings.semicolon);
											}
											if (snpAnnotation[a] != null) {
												if (annotStr.length() == 0) {
													annotStr += snpAnnotation[a].getName();
												} else {
													annotStr += ";" + snpAnnotation[a].getName();
												}
												if (geneStr != null && geneStr.length() > 0) {
													annotStr += " (" + geneStr + ")";
												}
											}
										}
									}
								}
								
								if (regionSignificant[datasetId][regionId] && crediblesets[datasetId][regionId].length <= maxNrVariantsInCredibleSet) {
									outputregioninbed = true;
								}
								
								AssociationResult r = resultsPerDs.get(datasetId).get(snpId); //data[datasetId][regionId][snpId];
								ln += "\t" + regionSignificant[datasetId][regionId]
										+ "\t" + crediblesets[datasetId][regionId].length
										+ "\t" + sumsperregion[datasetId]
										+ "\t" + r.getSnp().toString()
										+ "\t" + Strings.concat(r.getSnp().getAlleles(), Strings.comma)
										+ "\t" + r.getSnp().getImputationQualityScore()
										+ "\t" + (1 - r.getSnp().getAFControls())
										+ "\t" + (1 - r.getSnp().getAFCases())
										+ "\t" + Strings.concat(r.getORs()[0], Strings.semicolon)
										+ "\t" + r.getLog10Pval()
										+ "\t" + r.getPosterior()
										+ "\t" + (r.getPval() < variantSignificanceThreshold)
										+ "\t" + annotStr;
								regionsums[datasetId] += r.getPosterior();
							}
						} else {
							
							for (int datasetId = 0; datasetId < data.length; datasetId++) {
								String annotStr = "";
								if (variantAnnotations[datasetId].length > snpId) {
									SNPClass[] snpAnnotation = variantAnnotations[datasetId][snpId];
									if (snpAnnotation != null) {
										for (int a = 0; a < snpAnnotation.length; a++) {
											String geneStr = null;
											if (a == 0) {
												geneStr = Strings.concat(variantGeneAnnotations[datasetId][snpId][0], Strings.semicolon);
											}
											if (a == 1) {
												geneStr = Strings.concat(variantGeneAnnotations[datasetId][snpId][1], Strings.semicolon);
											}
											if (snpAnnotation[a] != null) {
												if (annotStr.length() == 0) {
													annotStr += snpAnnotation[a].getName();
												} else {
													annotStr += ";" + snpAnnotation[a].getName();
												}
												if (geneStr != null && geneStr.length() > 0) {
													annotStr += " (" + geneStr + ")";
												}
											}
											
										}
									}
								}
								
								if (datasetId == 0) {
									ln = "\t\t\t\t\t";
								} else {
									ln += "\t\t\t\t";
								}
								
								AssociationResult r = resultsPerDs.get(datasetId).get(snpId); //data[datasetId][regionId][snpId];
								if (regionsums[datasetId] < maxPosteriorCredibleSet) {
									ln += r.getSnp().toString()
											+ "\t" + Strings.concat(r.getSnp().getAlleles(), Strings.comma)
											+ "\t" + r.getSnp().getImputationQualityScore()
											+ "\t" + (1 - r.getSnp().getAFControls())
											+ "\t" + (1 - r.getSnp().getAFCases())
											+ "\t" + Strings.concat(r.getORs()[0], Strings.semicolon)
											+ "\t" + r.getLog10Pval()
											+ "\t" + r.getPosterior()
											+ "\t" + (r.getPval() < variantSignificanceThreshold)
											+ "\t" + annotStr;
								} else {
									ln += "\t"
											+ "\t"
											+ "\t"
											+ "\t"
											+ "\t"
											+ "\t"
											+ "\t"
											+ "\t"
											+ "\t";
								}
								regionsums[datasetId] += r.getPosterior();
							}
							
						}
						out.writeln(ln);
					}
					
				}
			}
			
			if (outputregioninbed) {
				outr.writeln(region.toBedString());
			}
			
			
		}
		outi.close();
		outa.close();
		outg.close();
		outr.close();
		out.close();
		
	}
	
	private boolean hasSignificantResult(ArrayList<AssociationResult> allDatasetData, double variantSignificanceThreshold) {
		for (AssociationResult r : allDatasetData) {
			if (r.getPval() < variantSignificanceThreshold) {
				return true;
			}
		}
		return false;
	}
	
	private ArrayList<AssociationResult> filter(ArrayList<AssociationResult> associationResults, Feature feature) {
		ArrayList<AssociationResult> output = new ArrayList<>();
		for (AssociationResult r : associationResults) {
			if (r.getSnp().overlaps(feature)) {
				output.add(r);
			}
		}
		return output;
	}
	
	
	private void makeCircularPlot(String bedregions,
								  String[] assocFiles,
								  String[] datasetNames,
								  String genenamefile,
								  String outfile,
								  double significanceThreshold,
								  boolean onlyIncludeVariantsBelowSignificanceThreshold,
								  double maxPosteriorCredibleSet,
								  int nrVariantsInCredibleSet,
								  String annot) throws IOException, DocumentException {
		
		
		AssociationFile assocFile = new AssociationFile();
		Annotation annotation = null;
		if (annot.endsWith(".gtf.gz") || annot.endsWith(".gtf")) {
			annotation = new GTFAnnotation(annot);
		} else {
			annotation = new EnsemblStructures(annot);
		}
		TreeSet<Gene> genes = annotation.getGeneTree();
		
		HashMap<String, String> locusToGene = loadLocusToGene(genenamefile);
		
		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bedregions);
		
		// region variant1 pval1 posterior1 variant2 pval2 posterior2 variant3 pval3 posterior3
		AssociationFile f = new AssociationFile();
		
		
		// [dataset][region][variants]
		AssociationResult[][][] crediblesets = new AssociationResult[assocFiles.length][regions.size()][];
		AssociationResult[][][] data = new AssociationResult[assocFiles.length][regions.size()][];
		
		
		ApproximateBayesPosterior abp = new ApproximateBayesPosterior();
		ArrayList<Feature> regionsWithCredibleSets = new ArrayList<>();
		
		
		{
			ArrayList<ArrayList<AssociationResult>> allResults = new ArrayList<>();
			for (int i = 0; i < assocFiles.length; i++) {
				ArrayList<AssociationResult> allDatasetData = f.read(assocFiles[i]);
				allResults.add(allDatasetData);
			}
			
			for (int regionNr = 0; regionNr < regions.size(); regionNr++) {
				boolean hasSet = false;
				for (int diseaseNr = 0; diseaseNr < assocFiles.length; diseaseNr++) {
					if (datasetNames[diseaseNr].equals("RA") && regions.get(regionNr).getChromosome().equals(Chromosome.THREE)
							&& regions.get(regionNr).getStart() == 58154177) {
						System.out.println("Got it");
					}
					ArrayList<AssociationResult> allDatasetData = filter(allResults.get(diseaseNr), regions.get(regionNr));
					data[diseaseNr][regionNr] = allDatasetData.toArray(new AssociationResult[0]);
					ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(allDatasetData, maxPosteriorCredibleSet);
					
					boolean abovethresh = false;
					for (AssociationResult r : credibleSet) {
						if (r.getPval() < significanceThreshold) {
							abovethresh = true;
						}
					}
					
					
					crediblesets[diseaseNr][regionNr] = credibleSet.toArray(new AssociationResult[0]);
					if (credibleSet.size() <= nrVariantsInCredibleSet && abovethresh) {
						hasSet = true;
					}
				}
				if (hasSet) {
					regionsWithCredibleSets.add(regions.get(regionNr));
				}
			}
		}
		
		HashSet<Feature> regionsWithCredibleSetsHash = new HashSet<Feature>();
		regionsWithCredibleSetsHash.addAll(regionsWithCredibleSets);
		
		int len = nrVariantsInCredibleSet;
		
		double[][][] dataForPlotting = new double[data.length][regionsWithCredibleSets.size()][];
		double[][][] dataForPlotting2 = new double[data.length][regionsWithCredibleSets.size()][];
		boolean[][] datasetHasCredibleSetInRegion = new boolean[data.length][regionsWithCredibleSets.size()];
		
		SNPClass[][][] snpClasses = new SNPClass[regionsWithCredibleSets.size()][][]; // [regions][snps][annotations]
		String[][] snpNames = new String[regionsWithCredibleSets.size()][];
		
		int regionCtr = 0;
		ArrayList<Triple<Integer, Integer, String>> groups = new ArrayList<Triple<Integer, Integer, String>>();
		
		ArrayList<String> groupnames = new ArrayList<>();
		int prevCtr = 0;
		for (int regionId = 0; regionId < regions.size(); regionId++) {
			Feature region = regions.get(regionId);
			if (regionsWithCredibleSetsHash.contains(region)) {
				// region nrCrediblesetVariants posteriorsumtop5 topvariants alleles or pval posterior
				ArrayList<ArrayList<AssociationResult>> resultsPerDs = new ArrayList<>();
				HashMap<String, Integer> variantToInt = new HashMap<String, Integer>();
				
				ArrayList<String> variantsNamesInRegion = new ArrayList<String>();
				ArrayList<SNPFeature> variantsInRegion = new ArrayList<SNPFeature>();
				for (int diseaseNr = 0; diseaseNr < data.length; diseaseNr++) {
//					ArrayList<AssociationResult> topResults = getTopVariants(data, diseaseNr, regionId, len);
//					resultsPerDs.add(topResults);
					AssociationResult[] credibleset = crediblesets[diseaseNr][regionId];
					ArrayList<AssociationResult> r = new ArrayList<>();
					for (int i = 0; i < credibleset.length; i++) {
						r.add(credibleset[i]);
					}
					resultsPerDs.add(r);
					
					double sum = 0;
					if (credibleset.length <= len) {
						for (int s = 0; (s < credibleset.length && s < len); s++) {
							double posterior = credibleset[s].getPosterior();
							if (sum < maxPosteriorCredibleSet) {
								String variant = credibleset[s].getSnp().getName();
								boolean isSignificant = (credibleset[s].getPval() < significanceThreshold);
								if (!variantToInt.containsKey(variant) && (isSignificant || !onlyIncludeVariantsBelowSignificanceThreshold)) {
									variantToInt.put(variant, variantToInt.size());
									variantsInRegion.add(credibleset[s].getSnp());
									variantsNamesInRegion.add(variant);
								}
							}
							sum += posterior;
						}
					}
//					for (int s = 0; s < topResults.size(); s++) {
//						double posterior = topResults.get(s).getPosterior();
//						if (sum < maxPosteriorCredibleSet) {
//							String variant = topResults.get(s).getSnp().getName();
//							boolean isSignificant = (topResults.get(s).getPval() < significanceThreshold);
//							if (!variantToInt.containsKey(variant) && isSignificant) {
//								variantToInt.put(variant, variantToInt.size());
//								variantsInRegion.add(topResults.get(s).getSnp());
//								variantsNamesInRegion.add(variant);
//							}
//						}
//						sum += posterior;
//					}
				}
				
				snpNames[regionCtr] = variantsNamesInRegion.toArray(new String[0]);
				
				Pair<SNPClass[][], String[][][]> annotationPair = annotateVariants(variantsInRegion, region, genes);
				snpClasses[regionCtr] = annotationPair.getLeft();
				String[][][] genenames = annotationPair.getRight(); // [variant][coding/promoter][genename]
				
				for (int datasetNr = 0; datasetNr < data.length; datasetNr++) {
					dataForPlotting[datasetNr][regionCtr] = new double[variantToInt.size()];
					dataForPlotting2[datasetNr][regionCtr] = new double[variantToInt.size()];
					
					// fill with nans
					for (int q = 0; q < variantToInt.size(); q++) {
						dataForPlotting[datasetNr][regionCtr][q] = Double.NaN;
						dataForPlotting2[datasetNr][regionCtr][q] = Double.NaN;
					}
					
					AssociationResult[] credibleset = crediblesets[datasetNr][regionId];
					
					boolean abovethresh = false;

//					ArrayList<AssociationResult> dsResults = resultsPerDs.get(i);
					AssociationResult[] allDsResults = data[datasetNr][regionId];
					for (int s = 0; s < allDsResults.length; s++) {
						
						AssociationResult r = allDsResults[s];
//						String variant = dsResults.get(s).getSnp().getName();
//						Integer variantId = variantToInt.get(variant);
//						double posterior = dsResults.get(s).getPosterior();
//						if (dsResults.get(s).getPval() < significanceThreshold) {
//							abovethresh = true;
//						}
						String variant = r.getSnp().getName();
						Integer variantId = variantToInt.get(variant);
						double posterior = r.getPosterior();
						if (r.getPval() < significanceThreshold) {
							abovethresh = true;
						}
						if (variantId != null) {
							if (variantId > dataForPlotting[datasetNr][regionCtr].length - 1) {
								System.out.println();
								System.out.println("WEIRDNESSSSSSSSSSSSSSSSS: " + variant);
								System.out.println();
							} else {
								if (r.getPval() < significanceThreshold) {
									dataForPlotting[datasetNr][regionCtr][variantId] = posterior;
									dataForPlotting2[datasetNr][regionCtr][variantId] = 1;
								} else {
									dataForPlotting[datasetNr][regionCtr][variantId] = Double.NaN;
									dataForPlotting2[datasetNr][regionCtr][variantId] = Double.NaN;
								}
							}
						}
					}
					
					if (abovethresh && credibleset.length <= nrVariantsInCredibleSet) {
						datasetHasCredibleSetInRegion[datasetNr][regionCtr] = true;
					}
					
				}
				
				DecimalFormat decimalFormat = new DecimalFormat("#");
				decimalFormat.setGroupingUsed(true);
				decimalFormat.setGroupingSize(3);
				
				String locusName = region.getChromosome().toString() + ":" + decimalFormat.format(region.getStart()) + "-" + decimalFormat.format(region.getStop());
				// locusName = region.toString();
				locusName += "; " + locusToGene.get(region.toString());
				groupnames.add(locusName);
				
				groups.add(new Triple<Integer, Integer, String>(regionCtr, regionCtr + 1, locusName));
				System.out.println("region\t" + region.toString() + "\t" + locusName);
				regionCtr++;
			}
		}
		
		
		Grid grid = new Grid(600, 600, 1, 1, 100, 0);
		CircularHeatmapPanel panel = new CircularHeatmapPanel(1, 1);
		panel.setRange(new Range(0, 0, 1, 1));
		panel.setData(datasetNames, groupnames.toArray(new String[0]), dataForPlotting);
		panel.setBinAnnotations(snpClasses);
		panel.setGroupAnnotations(datasetHasCredibleSetInRegion);
		panel.setGroups(groups);
		grid.addPanel(panel);
		grid.draw(outfile);

//		grid = new Grid(1000, 1000, 1, 1, 0, 0);
//		panel = new CircularHeatmapPanel(1, 1);
//		panel.setData(datasetNames, groupnames.toArray(new String[0]), dataForPlotting2);
//		panel.setBinAnnotations(snpClasses);
//		panel.setGroups(groups);
//		grid.addPanel(panel);
//		grid.draw(outfile + "-bin.pdf");
	
	}
	
	private Pair<SNPClass[][], String[][][]> annotateVariants(ArrayList<SNPFeature> variantsInRegion, Feature region, TreeSet<Gene> genes) {
		
		
		ArrayList<Gene> overlappingGenes = new ArrayList<>();
		for (Gene g : genes) {
			if (g.overlaps(region)) {
//				System.out.println(g.toString());
				overlappingGenes.add(g);
			}
		}
		SNPClass[][] output = new SNPClass[variantsInRegion.size()][3];
		
		String[][][] genenames = new String[variantsInRegion.size()][2][];
		
		// gene overlap
		for (int i = 0; i < variantsInRegion.size(); i++) {
			HashSet<String> exonGenesForVariant = new HashSet<>();
			HashSet<String> promoterGenesForVariant = new HashSet<>();
			boolean isGenic = false;
			boolean isUTR = false;
			boolean isExonic = false;
			boolean isPromoter = false;
			
			for (Gene g : overlappingGenes) {
				
				
				SNPFeature variant = variantsInRegion.get(i);
				if (variant.overlaps(g)) {
					// get exons
					isGenic = true;
					exonGenesForVariant.add(g.getGeneSymbol());
					ArrayList<Transcript> transcripts = g.getTranscripts();
					for (Transcript t : transcripts) {
						ArrayList<UTR> utrs = t.getUTRs();
						
						// determine if the exon overlaps an UTR
						boolean overlapsUTR = false;
						if (utrs != null) {
							for (UTR u : utrs) {
								if (u.overlaps(variant)) {
									isUTR = true;
									overlapsUTR = true;
								}
							}
						}
						
						// might not overlap the UTR of another transcript..
						if (!overlapsUTR) {
							ArrayList<Exon> exons = t.getExons();
							for (Exon e : exons) {
								if (e.overlaps(variant)) {
									isExonic = true;
								}
							}
						}
					}
				}
				
				// check whether it is in the promotor
				Feature promotor = null;
				if (g.getStrand().equals(Strand.POS)) {
					promotor = new Feature(g.getChromosome(), g.getStart() - promotordistance, g.getStop());
				} else {
					promotor = new Feature(g.getChromosome(), g.getStart(), g.getStop() + promotordistance);
				}
				if (promotor.overlaps(variant)) {
					isPromoter = true;
					promoterGenesForVariant.add(g.getGeneSymbol());
				}
				
			}
			
			// give priority to coding, exonic annotation, if for example the variant overlaps two genes.
//			if (output[i][0] == null || !output[i][0].equals(SNPClass.EXONIC)) {
			if (isGenic) {
				if (isExonic) {
					output[i][0] = SNPClass.EXONIC;
				} else if (isUTR) {
					output[i][0] = SNPClass.UTR;
				} else {
					output[i][0] = SNPClass.INTRONIC;
				}
			} else {
				output[i][0] = SNPClass.NONCODING;
				if (isPromoter) {
					output[i][1] = SNPClass.PROMOTER;
				}
			}
			
			genenames[i][0] = exonGenesForVariant.toArray(new String[0]);
			genenames[i][1] = promoterGenesForVariant.toArray(new String[0]);
			
		}
		
		// assign indel status
		for (int i = 0; i < variantsInRegion.size(); i++) {
			SNPFeature variant = variantsInRegion.get(i);
			String[] alleles = variant.getAlleles();
			boolean indel = false;
			for (String s : alleles) {
				if (s.length() > 1) {
					indel = true;
				}
			}
			if (indel) {
				output[i][2] = SNPClass.INDEL;
			}
		}
		
		return new Pair<>(output, genenames);
	}
	
	private ArrayList<AssociationResult> getTopVariants(AssociationResult[][][] data, int i, int d, int len) {
		ArrayList<AssociationResultPValuePair> pairs = new ArrayList<AssociationResultPValuePair>();
		
		for (AssociationResult r : data[i][d]) {
			AssociationResultPValuePair p = new AssociationResultPValuePair(r, r.getPosterior(), false);
			if (!Double.isInfinite(p.getP()) && !Double.isNaN(p.getP())) {
				pairs.add(p);
			}
		}
		
		ArrayList<AssociationResult> credibleSet = new ArrayList<AssociationResult>();
		if (!pairs.isEmpty()) {
			Collections.sort(pairs);
			ArrayList<AssociationResult> output = new ArrayList<>();
			int ctr = 0;
			while (ctr < len) {
				output.add(pairs.get(ctr).getAssociationResult());
				ctr++;
			}
			return output;
		}
		
		return null;
	}
}

class AnnotationData {
	
	ArrayList<String> names = new ArrayList<>();
	ArrayList<String> files = new ArrayList<>();
	ArrayList<String> groups = new ArrayList<String>();
	ArrayList<String> uniqueGroupNames = new ArrayList<String>();
	ArrayList<ArrayList<Feature>> annotations = new ArrayList<>();
	HashMap<String, ArrayList<Integer>> groupToFile = new HashMap<>();
	
	public AnnotationData(String filename, ArrayList<Feature> regions) throws IOException {
		TextFile tf1 = new TextFile(filename, TextFile.R);
		String[] elems = tf1.readLineElems(TextFile.tab); // file\tname\tgroup
		BedFileReader reader = new BedFileReader();
		while (elems != null) {
			if (elems.length == 2) {
				// file\tname
				names.add(elems[1]);
				groups.add("None");
			} else if (elems.length == 3) {
				// file\tname\tgroup
				names.add(elems[1]);
				groups.add(elems[2]);
			} else {
				groups.add("None");
				names.add(elems[0]);
			}
			if (elems.length >= 1 && elems[0].trim().length() > 0) {
				files.add(elems[0]);
				ArrayList<Feature> tmp = reader.readAsList(elems[0], regions);
				System.out.println(tmp.size() + " annotations overlap in: " + elems[0]);
				annotations.add(tmp);
			}
			
			
			elems = tf1.readLineElems(TextFile.tab); // file\tname\tgroup
		}
		tf1.close();
		System.out.println(files.size() + " files in " + filename);
		
		for (int i = 0; i < files.size(); i++) {
			ArrayList<Integer> ids = groupToFile.get(groups.get(i));
			if (ids == null) {
				ids = new ArrayList<>();
			}
			ids.add(i);
			groupToFile.put(groups.get(i), ids);
		}
		System.out.println(groupToFile.size() + " unique annotation groups.");
		
		HashSet<String> grouphash = new HashSet<>();
		grouphash.addAll(groups);
		ArrayList<String> groupsId = new ArrayList<>();
		groupsId.addAll(grouphash);
		Collections.sort(groupsId);
		uniqueGroupNames = groupsId;
	}
	
	public ArrayList<String> getUniqueGroups() {
		return uniqueGroupNames;
	}
	
	public int countOverlappingAnnotations(Feature locus, int group) {
		String name = uniqueGroupNames.get(group);
		ArrayList<Integer> ids = groupToFile.get(name);
		
		
		int ctr = 0;
		for (int g = 0; g < ids.size(); g++) {
			Integer id = ids.get(g);
			ArrayList<Feature> annotation = annotations.get(id);
			
			boolean overlap = overlap(annotation, locus);
			if (overlap) {
				ctr++;
			}
		}
		
		return ctr;
		
	}
	
	
	public boolean overlap(ArrayList<Feature> annotation, Feature query) {
		for (Feature f : annotation) {
			if (f.overlaps(query)) {
				return true;
			}
		}
		return false;
	}
	
	public int getNrAnnotationsInGroup(int group) {
		String name = uniqueGroupNames.get(group);
		ArrayList<Integer> ids = groupToFile.get(name);
		return ids.size();
	}
}
