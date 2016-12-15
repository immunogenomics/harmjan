package nl.harmjanwestra.finemapping;

import com.itextpdf.text.DocumentException;
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
import nl.harmjanwestra.utilities.graphics.panels.CircularHeatmapPanel;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFTabix;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
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


			String bedregions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
//		String bedregions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/RA-significantloci-75e7.bed";
			String genenames = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/AllLoci-GenesPerLocus.txt";
			String geneAnnotation = "/Data/Ref/Ensembl/GrCH37-b86-Structures.txt.gz"; //"/Data/Ref/Annotation/UCSC/genes.gtf.gz";
			String outfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/MergedCredibleSets/mergedCredibleSets.txt";
			String outeqtlfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/MergedCredibleSets/mergedCredibleSets-eqtls.txt";
			String outoverlapfile = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/MergedCredibleSets/mergedCredibleSets-overlap.txt";
			String outoverlapplot = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/MergedCredibleSets/annotationplots/mergedCredibleSets-";
			String outplot = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/MergedCredibleSets/mergedCredibleSets-plot-promoter" + promotordistance + ".pdf";
			String[] assocfiles = new String[]{
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/RA-assoc0.3-COSMO-merged-posterior.txt.gz",
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
					"/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/2016-09-06-SummaryStats/NormalHWEP1e4/META-assoc0.3-COSMO-merged-posterior.txt.gz"
			};
			if (windows) {
				bedregions = "D:\\Cloud\\Dropbox\\2016-03-RAT1D-Finemappng\\Data\\LocusDefinitions\\AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
				genenames = "D:\\Cloud\\Dropbox\\2016-03-RAT1D-Finemappng\\Data\\AllLoci-GenesPerLocus.txt";
				geneAnnotation = "D:\\Cloud\\Dropbox\\2016-03-RAT1D-Finemappng\\Data\\genes.gtf.gz";

				outfile = "D:\\Cloud\\Dropbox\\2016-03-RAT1D-Finemappng\\Data\\2016-09-06-SummaryStats\\NormalHWEP1e4\\MergedCredibleSets\\mergedCredibleSets.txt";
				outplot = "D:\\Cloud\\Dropbox\\2016-03-RAT1D-Finemappng\\Data\\2016-09-06-SummaryStats\\NormalHWEP1e4\\MergedCredibleSets\\mergedCredibleSets-plot-2k5promoter.pdf";
				assocfiles = new String[]{
						"D:\\Cloud\\Dropbox\\2016-03-RAT1D-Finemappng\\Data\\2016-09-06-SummaryStats\\NormalHWEP1e4\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
						"D:\\Cloud\\Dropbox\\2016-03-RAT1D-Finemappng\\Data\\2016-09-06-SummaryStats\\NormalHWEP1e4\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
						"D:\\Cloud\\Dropbox\\2016-03-RAT1D-Finemappng\\Data\\2016-09-06-SummaryStats\\NormalHWEP1e4\\META-assoc0.3-COSMO-merged-posterior.txt.gz"
				};
			}


			String[] datasetnames = new String[]{
					"RA",
					"T1D",
					"Combined"
			};
			double threshold = 7.5E-7;
			int nrVariantsInCredibleSet = 10;
			double maxPosteriorCredibleSet = 0.9;
			boolean includeAllLoci = true;

			c.mergeCredibleSets(bedregions, assocfiles, datasetnames, genenames, outfile, maxPosteriorCredibleSet, threshold, nrVariantsInCredibleSet, geneAnnotation, includeAllLoci);

//			c.makeCircularPlot(bedregions, assocfiles, datasetnames, genenames, outplot, threshold, maxPosteriorCredibleSet, nrVariantsInCredibleSet, geneAnnotation);
//			System.exit(-1);

			String[] eqtlfiles = new String[]{
					"/Data/eQTLs/ImmVar/Raj/tableS12_meta_cd4T_cis_fdr05-upd.tab",
					"/Data/eQTLs/BiosEQTLs/eQTLsFDR0.05-ProbeLevel.txt.gz"
			};
			String[] eqtlfilenames = new String[]{
					"Raj",
					"Bios"
			};
//			String dbsnpvcf = "/Data/Ref/dbSNP/human_9606_b147_GRCh37p13/00-All.vcf.gz";
//			c.rewriteEQTLFile("/Data/eQTLs/ImmVar/Raj/tableS12_meta_cd4T_cis_fdr05.tsv", dbsnpvcf,
//					"/Data/eQTLs/ImmVar/Raj/tableS12_meta_cd4T_cis_fdr05-upd.tab");

			String tabixprefix = "/Data/Ref/beagle_1kg/1kg.phase3.v5a.chrCHR.vcf.gz";
			String tabixfilter = "/Data/Ref/1kg-europeanpopulations.txt.gz";
			c.eQTLOverlap(bedregions, eqtlfiles, eqtlfilenames, tabixprefix, tabixfilter, assocfiles, datasetnames, genenames, outeqtlfile, maxPosteriorCredibleSet, nrVariantsInCredibleSet, geneAnnotation);

			String[] bedfiles = new String[]{
					"/Data/Enhancers/Roadmap/dnase-groups.txt",
					"/Data/Enhancers/Roadmap/h3k4me3-groups.txt",
					"/Data/Enhancers/ChromHMM/ChromHMMEnhancers-groups.txt",
					"/Data/Enhancers/ChromHMM/ChromHMMPromotors-groups.txt"
			};
			String[] bedfilenames = new String[]{
					"DNASE", "H3K4me3", "ChromHMM-Enhancers", "ChromHMM-Promotors"
			};

//			c.bedOverlapOld(bedregions,
//					bedfiles,
//					bedfilenames,
//					assocfiles, datasetnames, genenames,
//					outoverlapfile, maxPosteriorCredibleSet, nrVariantsInCredibleSet, geneAnnotation);

			c.bedOverlap(bedregions,
					bedfiles,
					bedfilenames,
					assocfiles, datasetnames, genenames,
					outoverlapfile, outoverlapplot,
					maxPosteriorCredibleSet, nrVariantsInCredibleSet, geneAnnotation);
//
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
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
							String geneAnnotation) throws IOException, DocumentException {

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
					Grid grid = new Grid(1000, 100, 1, 1, 100, 100);
					AnnotationTrackPanel p = new AnnotationTrackPanel(1, 1);

					int[] highlight = new int[crediblesets[ds][regionId].length];
					for (int i = 0; i < highlight.length; i++) {
						highlight[i] = crediblesets[ds][regionId][i].getSnp().getStart();
					}
					String[][] groupnames = null;
					p.setData(overlap, highlight, region, groupnames);
					p.setMarginX(10);
					p.setMarginY(10);
					grid.addPanel(p);
					grid.draw(outplot + datasetnames[ds] + "-" + region.toString() + ".pdf");
				}

//
//				qctr++;
//				System.out.println(ctr + " regions processed.. ");

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
		String header2 = "\t\t";

		String header1 = "region\tgene";
		for (int i = 0; i < data.length; i++) {
			header2 += datasetNames[i]
					+ "\t\t\t\t\t";
			header1 += "\tNrVariantsInCredibleSet" +
					"\tVariants" +
					"\tPosterior";
			for (int e = 0; e < eqtlfilenames.length; e++) {
				header1 += "\teSNP-" + eqtlfilenames[e]
						+ "\teGene-" + eqtlfilenames[e]
						+ "\tP(eQTL)-" + eqtlfilenames[e]
						+ "\tDistance-" + eqtlfilenames[e]
						+ "\tLD-" + eqtlfilenames[e];
				header2 += "\t\t\t\t";
			}
		}
		out.writeln(header2);
		out.writeln(header1);

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
				eqtlregion.setStart(eqtlregion.getStart() - 1000000);
				eqtlregion.setStop(eqtlregion.getStop() + 1000000);
				ArrayList<VCFVariant> all1kgvariants = tabix.getAllVariants(region, filter);

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


								String lnblock = eqtllineblock(all1kgvariants,
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


								String lnblock = eqtllineblock(all1kgvariants,
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
											+ "\t";
									for (int e = 0; e < eqtlfilenames.length; e++) {
										ln += "\t\t\t\t\t";
									}
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
			if (filename.endsWith("tab")) {

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
					allEQTLs.add(eqtl);
					elems = tf.readLineElems(TextFile.tab);
				}
				tf.close();


			} else {

				// BIOS eQTL file format
				TextFile tf = new TextFile(eqtlfilenames[d], TextFile.R);
				tf.readLine();
				String[] elems = tf.readLineElems(TextFile.tab);
				int ctr = 0;
				while (elems != null) {

					String snp = elems[1];
					Double pval = Double.parseDouble(elems[0]);
					Chromosome chr = Chromosome.parseChr(elems[2]);
					Integer pos = Integer.parseInt(elems[3]);
					SNPFeature feature = new SNPFeature(chr, pos, pos);
					feature.setName(snp);

					EQTL eqtl = new EQTL();
					eqtl.setSnp(feature);
					eqtl.setGenename(elems[16]);
					eqtl.setName(elems[16]);
					eqtl.setPval(pval);

					for (int r = 0; r < eqtlregions.size(); r++) {

						Feature eqtlregion = eqtlregions.get(r);
						if (eqtl.getSnp().overlaps(eqtlregion)) {
							allEQTLs.add(eqtl);
						}
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
					"\tAFCases" +
					"\tAFControls" +
					"\tOR" +
					"\tPval" +
					"\tPosterior" +
					"\tSignificant" +
					"\tAnnotation";
		}

		out.writeln(header2);
		out.writeln(header1);
		for (int regionId = 0; regionId < regions.size(); regionId++) {
			Feature region = regions.get(regionId);
			if (regionsWithCredibleSetsHash.contains(region)) {
				// region nrCrediblesetVariants posteriorsumtop5 topvariants alleles or pval posterior
				ArrayList<ArrayList<AssociationResult>> resultsPerDs = new ArrayList<>();
				SNPClass[][][] variantAnnotations = new SNPClass[data.length][][];
				String[][][][] variantGeneAnnotations = new String[data.length][][][];
				for (int i = 0; i < data.length; i++) {
					ArrayList<AssociationResult> topResults = getTopVariants(data, i, regionId, len);
					resultsPerDs.add(topResults);

					ArrayList<SNPFeature> variants = new ArrayList<>();

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


								AssociationResult r = resultsPerDs.get(datasetId).get(snpId); //data[datasetId][regionId][snpId];
								ln += "\t" + regionSignificant[datasetId][regionId]
										+ "\t" + crediblesets[datasetId][regionId].length
										+ "\t" + sumsperregion[datasetId]
										+ "\t" + r.getSnp().toString()
										+ "\t" + Strings.concat(r.getSnp().getAlleles(), Strings.comma)
										+ "\t" + r.getSnp().getImputationQualityScore()
										+ "\t" + r.getSnp().getAFCases()
										+ "\t" + r.getSnp().getAFControls()
										+ "\t" + Strings.concat(r.getORs(), Strings.semicolon)
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
											+ "\t" + r.getSnp().getAFCases()
											+ "\t" + r.getSnp().getAFControls()
											+ "\t" + Strings.concat(r.getORs(), Strings.semicolon)
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


		}
		outi.close();
		outa.close();
		outg.close();
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
								  String[] assocFiles, String[] datasetNames, String genenamefile, String outfile, double threshold, double maxPosteriorCredibleSet, int nrVariantsInCredibleSet, String annot) throws IOException, DocumentException {


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

			for (int d = 0; d < regions.size(); d++) {
				boolean hasSet = false;
				for (int i = 0; i < assocFiles.length; i++) {
					ArrayList<AssociationResult> allDatasetData = filter(allResults.get(i), regions.get(d));
					data[i][d] = allDatasetData.toArray(new AssociationResult[0]);
					ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(allDatasetData, maxPosteriorCredibleSet);

					boolean abovethresh = false;
					for (AssociationResult r : credibleSet) {
						if (r.getPval() < threshold) {
							abovethresh = true;
						}
					}

					crediblesets[i][d] = credibleSet.toArray(new AssociationResult[0]);
					if (credibleSet.size() <= nrVariantsInCredibleSet && abovethresh) {
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

		int len = nrVariantsInCredibleSet;

		double[][][] dataForPlotting = new double[data.length][regionsWithCredibleSets.size()][];
		double[][][] dataForPlotting2 = new double[data.length][regionsWithCredibleSets.size()][];
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
				for (int i = 0; i < data.length; i++) {
					ArrayList<AssociationResult> topResults = getTopVariants(data, i, regionId, len);
					resultsPerDs.add(topResults);

					double sum = 0;
					for (int s = 0; s < topResults.size(); s++) {
						double posterior = topResults.get(s).getPosterior();
						if (sum < maxPosteriorCredibleSet) {
							String variant = topResults.get(s).getSnp().getName();
							boolean isSignificant = (topResults.get(s).getPval() < threshold);
							if (!variantToInt.containsKey(variant) && isSignificant) {
								variantToInt.put(variant, variantToInt.size());
								variantsInRegion.add(topResults.get(s).getSnp());
								variantsNamesInRegion.add(variant);
							}
						}
						sum += posterior;
					}
				}

				snpNames[regionCtr] = variantsNamesInRegion.toArray(new String[0]);

				Pair<SNPClass[][], String[][][]> annotationPair = annotateVariants(variantsInRegion, region, genes);
				snpClasses[regionCtr] = annotationPair.getLeft();
				String[][][] genenames = annotationPair.getRight(); // [variant][coding/promoter][genename]

				for (int i = 0; i < data.length; i++) {
					dataForPlotting[i][regionCtr] = new double[variantToInt.size()];
					dataForPlotting2[i][regionCtr] = new double[variantToInt.size()];

					// fill with nans
					for (int q = 0; q < variantToInt.size(); q++) {
						dataForPlotting[i][regionCtr][q] = Double.NaN;
						dataForPlotting2[i][regionCtr][q] = Double.NaN;
					}

					ArrayList<AssociationResult> dsResults = resultsPerDs.get(i);
					for (int s = 0; s < dsResults.size(); s++) {
						String variant = dsResults.get(s).getSnp().getName();
						Integer variantId = variantToInt.get(variant);
						double posterior = dsResults.get(s).getPosterior();

						if (variantId != null) {
							if (variantId > dataForPlotting[i][regionCtr].length - 1) {
								System.out.println();
								System.out.println("WEIRDNESSSSSSSSSSSSSSSSS: " + variant);
								System.out.println();
							} else {
								if (dsResults.get(s).getPval() < threshold) {
									dataForPlotting[i][regionCtr][variantId] = posterior;
									dataForPlotting2[i][regionCtr][variantId] = 1;
								} else {
									dataForPlotting[i][regionCtr][variantId] = Double.NaN;
									dataForPlotting2[i][regionCtr][variantId] = Double.NaN;
								}
							}

						}
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


		Grid grid = new Grid(1000, 1000, 1, 1, 100, 0);
		CircularHeatmapPanel panel = new CircularHeatmapPanel(1, 1);
		panel.setRange(new Range(0, 0, 1, 1));
		panel.setData(datasetNames, groupnames.toArray(new String[0]), dataForPlotting);
		panel.setAnnotations(snpClasses);
		panel.setGroups(groups);
		grid.addPanel(panel);
		grid.draw(outfile);

		grid = new Grid(1000, 1000, 1, 1, 100, 0);
		panel = new CircularHeatmapPanel(1, 1);
		panel.setData(datasetNames, groupnames.toArray(new String[0]), dataForPlotting2);
		panel.setAnnotations(snpClasses);
		panel.setGroups(groups);
		grid.addPanel(panel);
		grid.draw(outfile + "-bin.pdf");

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
