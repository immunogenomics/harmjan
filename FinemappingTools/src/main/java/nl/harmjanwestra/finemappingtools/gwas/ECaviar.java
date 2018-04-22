package nl.harmjanwestra.finemappingtools.gwas;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.utilities.annotation.ensembl.EnsemblStructures;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.math.stats.ZScores;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFTabix;
import nl.harmjanwestra.utilities.vcf.VCFVariant;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class ECaviar {
	
	
	public static void main(String[] args) {
		
		String disk = "d:";
		
		String[] eqtlfiles = new String[]{
				disk + "\\Data\\eQTLs\\ImmVar\\Raj\\tableS12_meta_cd4T_cis_fdr05-upd.tab",
				disk + "\\Data\\eQTLs\\Milani\\CD4-cis-eQTLs-ProbeLevelFDR0.5.txt.gz",
				disk + "\\Data\\eQTLs\\Milani\\CD8-cis-eQTLs-ProbeLevelFDR0.5.txt.gz",
				disk + "\\Data\\eQTLs\\BiosEQTLs\\eQTLsFDR0.05-ProbeLevel.txt.gz",
				disk + "\\Data\\eQTLs\\BluePrint\\mono_gene_nor_combat_peer_10_all_summary-fdr005.tab.gz",
				disk + "\\Data\\eQTLs\\BluePrint\\mono_K27AC_log2rpm_peer_10_all_summary-fdr005.tab.gz",
				disk + "\\Data\\eQTLs\\BluePrint\\mono_K4ME1_log2rpm_peer_10_all_summary-fdr005.tab.gz",
				disk + "\\Data\\eQTLs\\BluePrint\\mono_meth_M_peer_10_all_summary-fdr005.tab.gz",
				disk + "\\Data\\eQTLs\\BluePrint\\mono_psi_peer_10_all_summary-fdr005.tab.gz",
				disk + "\\Data\\eQTLs\\BluePrint\\neut_gene_nor_combat_peer_10_all_summary-fdr005.tab.gz",
				disk + "\\Data\\eQTLs\\BluePrint\\neut_K27AC_log2rpm_peer_10_all_summary-fdr005.tab.gz",
				disk + "\\Data\\eQTLs\\BluePrint\\neut_K4ME1_log2rpm_peer_10_all_summary-fdr005.tab.gz",
				disk + "\\Data\\eQTLs\\BluePrint\\neut_meth_M_peer_10_all_summary-fdr005.tab.gz",
				disk + "\\Data\\eQTLs\\BluePrint\\neut_psi_peer_10_all_summary-fdr005.tab.gz",
				disk + "\\Data\\eQTLs\\BluePrint\\tcel_gene_nor_combat_peer_10_all_summary-fdr005.tab.gz",
				disk + "\\Data\\eQTLs\\BluePrint\\tcel_K27AC_log2rpm_peer_10_all_summary-fdr005.tab.gz",
				disk + "\\Data\\eQTLs\\BluePrint\\tcel_K4ME1_log2rpm_peer_10_all_summary-fdr005.tab.gz",
				disk + "\\Data\\eQTLs\\BluePrint\\tcel_meth_M_peer_10_all_summary-fdr005.tab.gz",
				disk + "\\Data\\eQTLs\\BluePrint\\tcel_psi_peer_10_all_summary-fdr005.tab.gz",
				disk + "\\Data\\eQTLs\\Sun-pQTL\\table1.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Adipose_Subcutaneous_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Adipose_Visceral_Omentum_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Adrenal_Gland_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Artery_Aorta_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Artery_Coronary_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Artery_Tibial_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Brain_Anterior_cingulate_cortex_BA24_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Brain_Caudate_basal_ganglia_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Brain_Cerebellar_Hemisphere_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Brain_Cerebellum_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Brain_Cortex_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Brain_Frontal_Cortex_BA9_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Brain_Hippocampus_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Brain_Hypothalamus_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Brain_Nucleus_accumbens_basal_ganglia_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Brain_Putamen_basal_ganglia_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Breast_Mammary_Tissue_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Cells_EBV-transformed_lymphocytes_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Cells_Transformed_fibroblasts_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Colon_Sigmoid_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Colon_Transverse_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Esophagus_Gastroesophageal_Junction_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Esophagus_Mucosa_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Esophagus_Muscularis_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Heart_Atrial_Appendage_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Heart_Left_Ventricle_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Liver_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Lung_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Muscle_Skeletal_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Nerve_Tibial_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Ovary_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Pancreas_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Pituitary_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Prostate_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Skin_Not_Sun_Exposed_Suprapubic_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Skin_Sun_Exposed_Lower_leg_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Small_Intestine_Terminal_Ileum_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Spleen_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Stomach_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Testis_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Thyroid_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Uterus_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Vagina_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz",
				disk + "\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v6p_eQTL/Whole_Blood_Analysis.v6p.signif_snpgene_pairs.txt.gz.tab.gz"
		};
		
		String[] eqtlfilenames = new String[]{
				"Raj-Tcell",
				"Milani-CD4",
				"Milani-CD8",
				"Bios-WholeBlood",
				"Blueprint-Monocyte-eQTL",
				"Blueprint-Monocyte-hQTL-K27AC",
				"Blueprint-Monocyte-hQTL-K4ME1",
				"Blueprint-Monocyte-mQTL",
				"Blueprint-Monocyte-sQTL",
				"Blueprint-Neutrophil-eQTL",
				"Blueprint-Neutrophil-hQTL-K27AC",
				"Blueprint-Neutrophil-hQTL-K4ME1",
				"Blueprint-Neutrophil-mQTL",
				"Blueprint-Neutrophil-sQTL",
				"Blueprint-TCell-eQTL",
				"Blueprint-TCell-hQTL-K27AC",
				"Blueprint-TCell-hQTL-K4ME1",
				"Blueprint-TCell-mQTL",
				"Blueprint-TCell-sQTL",
				"Sun-pQTL",
				"GTEx-Adipose_Subcutaneous",
				"GTEx-Adipose_Visceral_Omentum",
				"GTEx-Adrenal_Gland",
				"GTEx-Artery_Aorta",
				"GTEx-Artery_Coronary",
				"GTEx-Artery_Tibial",
				"GTEx-Brain_Anterior_cingulate_cortex_BA24",
				"GTEx-Brain_Caudate_basal_ganglia",
				"GTEx-Brain_Cerebellar_Hemisphere",
				"GTEx-Brain_Cerebellum",
				"GTEx-Brain_Cortex",
				"GTEx-Brain_Frontal_Cortex_BA9",
				"GTEx-Brain_Hippocampus",
				"GTEx-Brain_Hypothalamus",
				"GTEx-Brain_Nucleus_accumbens_basal_ganglia",
				"GTEx-Brain_Putamen_basal_ganglia",
				"GTEx-Breast_Mammary_Tissue",
				"GTEx-Cells_EBV-transformed_lymphocytes",
				"GTEx-Cells_Transformed_fibroblasts",
				"GTEx-Colon_Sigmoid",
				"GTEx-Colon_Transverse",
				"GTEx-Esophagus_Gastroesophageal_Junction",
				"GTEx-Esophagus_Mucosa",
				"GTEx-Esophagus_Muscularis",
				"GTEx-Heart_Atrial_Appendage",
				"GTEx-Heart_Left_Ventricle",
				"GTEx-Liver",
				"GTEx-Lung",
				"GTEx-Muscle_Skeletal",
				"GTEx-Nerve_Tibial",
				"GTEx-Ovary",
				"GTEx-Pancreas",
				"GTEx-Pituitary",
				"GTEx-Prostate",
				"GTEx-Skin_Not_Sun_Exposed_Suprapubic",
				"GTEx-Skin_Sun_Exposed_Lower_leg",
				"GTEx-Small_Intestine_Terminal_Ileum",
				"GTEx-Spleen",
				"GTEx-Stomach",
				"GTEx-Testis",
				"GTEx-Thyroid",
				"GTEx-Uterus",
				"GTEx-Vagina",
				"GTEx-Whole_Blood"
		};


//		QTLOverlap o = new QTLOverlap();
//		String ensembl = "D:\\Data\\Ref\\Ensembl\\GrCH37-b86-Structures.txt.gz";
//		double ldthresh = 0.8;
//		int ciswindow = 1000000;
//		String variantfile = disk + "\\Sync\\Dropbox\\FineMap\\2018-01-Rebuttal\\tables\\listofsnpswithposterior0.2.txt";
//		String tabix = disk + "\\Data\\Ref\\1kg\\ALL.chrCHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";
//		String samplefile = disk + "\\Data\\Ref\\1kg-europeanpopulations.txt.gz";
//		String output = disk + "\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\qtloverlap\\output.txt";
//		try {
//			o.run(variantfile, eqtlfiles, eqtlfilenames, tabix, samplefile, ciswindow, ldthresh, ensembl, output);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
	
	}
	
	
	public void run(String variantfile,
					String[] qfiles,
					String[] qfilenames,
					String tabixprefix,
					String samplefile,
					int ciswindow,
					double ldthresh,
					String ens,
					String output) throws IOException {
		
		// read the snps
		
		
		ArrayList<SNPFeature> snps = new ArrayList<SNPFeature>();
		
		EnsemblStructures eg = new EnsemblStructures(ens);
		Collection<Gene> genes = eg.getGenes();
		HashMap<String, Gene> strToGene = new HashMap<>();
		for (Gene g : genes) {
			strToGene.put(g.getName(), g);
		}
		
		Collections.sort(snps, new FeatureComparator(true));
		ArrayList<Feature> regions = new ArrayList<Feature>();
		for (Feature f : snps) {
			Feature f2 = new Feature(f.getChromosome(), f.getStart() - ciswindow, f.getStop() + ciswindow);
			regions.add(f2);
		}
		
		// read qtls
		EQTL[][][] eqtls = loadEQTLs(qfiles, regions);
		DetermineLD ldcal = new DetermineLD();
		
		
		String[][] outputlns = new String[qfiles.length][snps.size()];
		
		// write header;
		ExecutorService ex = Executors.newFixedThreadPool(6);
		boolean[] done = new boolean[snps.size()];
		for (int s = 0; s < snps.size(); s++) {
//			ex.submit(new QTLOverLapTask(done, regions, s, snps, tabixprefix, samplefile, outputlns, eqtls, qfiles, ldthresh, strToGene));
		}
		
		while (!alldone(done)) {
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		ex.shutdown();
		
		
	}
	
	private boolean alldone(boolean[] d) {
		boolean alldone = true;
		int done = 0;
		for (boolean b : d) {
			
			if (!b) {
				alldone = false;
			} else {
				done++;
			}
		}
		System.out.println(done + " tasks done.");
		return alldone;
	}
	
	private class QTLOverLapTask implements Runnable {
		
		private final boolean[] done;
		private final ArrayList<AssociationResult> associationResults;
		ArrayList<Feature> regions;
		int s;
		String tabixprefix;
		String samplefile;
		String outputprefix;
		EQTL[][][] eqtls;
		String[] qfiles;
		double ldthresh;
		HashMap<String, Gene> strToGene;
		
		public QTLOverLapTask(boolean[] done, ArrayList<Feature> regions,
							  int s,
							  ArrayList<SNPFeature> snps,
							  String tabixprefix,
							  String samplefile,
							  String outputprefix,
							  EQTL[][][] eqtls,
							  ArrayList<AssociationResult> associationResults,
							  String[] qfiles,
							  double ldthresh,
							  HashMap<String, Gene> strToGene) {
			this.done = done;
			this.regions = regions;
			this.s = s;
			this.associationResults = associationResults;
			this.tabixprefix = tabixprefix;
			this.samplefile = samplefile;
			this.outputprefix = outputprefix;
			this.eqtls = eqtls;
			this.qfiles = qfiles;
			this.ldthresh = ldthresh;
			this.strToGene = strToGene;
		}
		
		private ArrayList<AssociationResult> getRegionAssocResults(ArrayList<AssociationResult> assocResults, Feature region) {
			ArrayList<AssociationResult> output = new ArrayList<>();
			for (int r = 0; r < assocResults.size(); r++) {
				AssociationResult result = assocResults.get(r);
				if (result.getSnp().overlaps(region)) {
					output.add(result);
				}
			}
			return output;
		}
		
		@Override
		public void run() {
			try {
				Feature region = regions.get(s);
				String tabixFile = tabixprefix.replaceAll("CHR", "" + region.getChromosome().getNumber());
				
				VCFTabix t = new VCFTabix(tabixFile);
				boolean[] filter = null;
				if (samplefile != null) {
					filter = t.getSampleFilter(samplefile);
				}
				ArrayList<VCFVariant> allvariants = t.getAllVariants(region, filter);
				
				for (int q = 0; q < qfiles.length; q++) {
					EQTL[] regioneqtls = eqtls[q][s];
					ArrayList<AssociationResult> regionAssocs = getRegionAssocResults(associationResults, region);
					
					if (regioneqtls == null || regioneqtls.length == 0 || regionAssocs.size() == 0) {
						// do nothing
					} else {
						
						// get all associations per gene
						HashMap<String, ArrayList<EQTL>> genes = new HashMap<String, ArrayList<EQTL>>();
						for (EQTL e : regioneqtls) {
							ArrayList<EQTL> es = genes.get(e.genename);
							if (es == null) {
								es = new ArrayList<>();
							}
							es.add(e);
							genes.put(e.getGenename(), es);
						}
						
						
						for (String gene : genes.keySet()) {
							// sort
							ArrayList<EQTL> es = genes.get(gene);
							Collections.sort(es, new FeatureComparator(true));
							
							// get the snps that overlap with the reference panel...
							ArrayList<EQTL> tmpqtlswithsnpspresent = new ArrayList<>();
							ArrayList<VCFVariant> tmpvcfVariants = new ArrayList<>();
							for (EQTL e : es) {
								VCFVariant v = getVariant(e.getSnp(), allvariants);
								if (v != null) {
									tmpqtlswithsnpspresent.add(e);
									tmpvcfVariants.add(v);
								}
							}
							
							// filter the associationresults for those variants
							ArrayList<AssociationResult> assocswithsnpspresent = new ArrayList<>();
							ArrayList<EQTL> qtlswithsnpspresent = new ArrayList<>();
							ArrayList<VCFVariant> vcfVariants = new ArrayList<>();
							for (int v = 0; v < vcfVariants.size(); v++) {
								VCFVariant var = vcfVariants.get(v);
								AssociationResult r = getAssoc(var, regionAssocs);
								if (r != null) {
									qtlswithsnpspresent.add(tmpqtlswithsnpspresent.get(v));
									vcfVariants.add(var);
									assocswithsnpspresent.add(r);
									
									// TODO: check whether allele codings are identical..
									// otherwise, FLIP Z!
								}
								
							}
							
							// write Z-file
							TextFile outz1 = new TextFile(outputprefix + "_" + region.toString() + "_" + gene + "-eqtl.z", TextFile.W);
							TextFile outz2 = new TextFile(outputprefix + "_" + region.toString() + "_" + gene + "-gwas.z", TextFile.W);
							for (int i = 0; i < vcfVariants.size(); i++) {
								EQTL e = qtlswithsnpspresent.get(i);
								AssociationResult r = assocswithsnpspresent.get(i);
								VCFVariant v = vcfVariants.get(i);
								outz1.writeln(v.asSNPFeature().toString() + "\t" + e.getZscore());
								outz1.writeln(v.asSNPFeature().toString() + "\t" + e.getZscore());
							}
							outz1.close();
							outz2.close();
							
							double[][] corMat = new double[vcfVariants.size()][vcfVariants.size()];
							for (int i = 0; i < vcfVariants.size(); i++) {
								corMat[i][i] = 1d;
								double[] x = toArray(vcfVariants.get(i));
								for (int j = i + 1; j < vcfVariants.size(); j++) {
									double[] y = toArray(vcfVariants.get(j));
									double r = JSci.maths.ArrayMath.correlation(x, y);
									if (r <= -1d) {
										r = -1d;
									}
									if (r >= 1) {
										r = 1d;
									}
									corMat[i][j] = r;
									corMat[j][i] = r;
								}
							}
							
							
							TextFile ldout = new TextFile(outputprefix + "_" + region.toString() + "_" + gene + ".ld", TextFile.W);
							System.out.println("LD info: " + ldout);
							for (int i = 0; i < vcfVariants.size(); i++) {
								StringBuilder ln = new StringBuilder();
								for (int j = 0; j < vcfVariants.size(); j++) {
									if (ln.length() == 0) {
										ln.append(corMat[i][j]);
									} else {
										ln.append(" ").append(corMat[i][j]);
									}
								}
								ldout.writeln(ln.toString());
							}
							ldout.close();
							
						}
						
						
					}
				}
			} catch (IOException e) {
				done[s] = true;
				e.printStackTrace();
			}
			done[s] = true;
		}
		
		private AssociationResult getAssoc(VCFVariant var, ArrayList<AssociationResult> regionAssocs) {
			for (AssociationResult r : regionAssocs) {
				if (r.getSnp().overlaps(var.asSNPFeature())) {
					return r;
				}
			}
			return null;
		}
		
		private double[] toArray(VCFVariant vcfVariant) {
			DoubleMatrix2D mat = vcfVariant.getDosagesAsMatrix2D();
			double[] output = mat.viewColumn(0).toArray();
			return output;
		}
		
	}
	
	
	public VCFVariant getVariant(Feature snpFeature, ArrayList<VCFVariant> allvariants) {
		for (VCFVariant v : allvariants) {
			if (v.asFeature().overlaps(snpFeature)) {
				return v;
			}
		}
		return null;
		
	}
	
	private EQTL[][][] loadEQTLs(String[] eqtlfilenames, ArrayList<Feature> regions) throws IOException {
		EQTL[][][] output = new EQTL[eqtlfilenames.length][regions.size()][];
		
		for (int d = 0; d < eqtlfilenames.length; d++) {
			System.out.println("Parsing: " + eqtlfilenames[d]);
			
			
			ArrayList<EQTL> allEQTLs = new ArrayList<>();
			String filename = eqtlfilenames[d];
			boolean gtex = false;
			TextFile tfh = new TextFile(filename, TextFile.R);
			String fh = tfh.readLine();
			if (fh.startsWith("variant_id")) {
				gtex = true;
			}
			tfh.close();
			if (gtex) {
				TextFile tf = new TextFile(eqtlfilenames[d], TextFile.R);
				tf.readLine();
//				variant_id      gene_id tss_distance    pval_nominal    slope   slope_se        slope_fpkm      slope_fpkm_se   pval_nominal_threshold  min_pval_nominal        pval_beta
				String[] elems = tf.readLineElems(TextFile.tab);
				while (elems != null) {
					
					
					// 1_984847_G_C_b37
					String[] snpelems = elems[0].split("_");
					
					Chromosome chr = Chromosome.parseChr(snpelems[0]);
					Integer pos = Integer.parseInt(snpelems[1]);
					String snpname = elems[2];
					String gene = elems[1];
					double pval = Double.parseDouble(elems[2]);
					
					SNPFeature feature = new SNPFeature(chr, pos, pos);
					feature.setName(snpname);
					feature.setAlleles(new String[]{snpelems[2], snpelems[3]});
					
					EQTL eqtl = new EQTL();
					eqtl.setSnp(feature);
					eqtl.setGenename(gene);
					eqtl.setName(gene);
					eqtl.setPval(pval);
					double z = ZScores.zToP(pval);
					Double slpe = Double.parseDouble(elems[4]);
					if (slpe < 0) {
						z *= -1;
					}
					eqtl.setZscore(z);
					if (eqtloverlap(eqtl, regions)) {
						allEQTLs.add(eqtl);
					}
					
					elems = tf.readLineElems(TextFile.tab);
				}
				tf.close();
			} else if (filename.endsWith("tab") || filename.endsWith("tab.gz")) {
				
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
					
					if (eqtloverlap(eqtl, regions)) {
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
					double z = Double.parseDouble(elems[10]);
					String[] alleles = elems[8].split("/");
					String directionalleles = elems[9];
					if (!directionalleles.equals(alleles[0])) {
						alleles = new String[]{alleles[1], alleles[0]};
						z *= -1;
					}
					eqtl.setZscore(z);
					eqtl.getSnp().setAlleles(alleles);
					
					
					boolean sig = true;
					if (fdrcol != -1) {
						Double fdr = Double.parseDouble(elems[fdrcol]);
						if (fdr > 0.05) {
							sig = false;
						}
					}
					if (sig && eqtloverlap(eqtl, regions)) {
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
			for (int r = 0; r < regions.size(); r++) {
				ArrayList<EQTL> overlapping = new ArrayList<>();
				Feature eqtlregion = regions.get(r);
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
	
	
}
