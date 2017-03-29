package nl.harmjanwestra.finemapping.annotation;

import nl.harmjanwestra.utilities.annotation.Annotation;
import nl.harmjanwestra.utilities.annotation.ensembl.EnsemblStructures;
import nl.harmjanwestra.utilities.annotation.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

/**
 * Created by hwestra on 11/14/16.
 */
public class CodingAndIndelEnrichment {


	public static void main(String[] args) {
		CodingAndIndelEnrichment c = new CodingAndIndelEnrichment();

		try {
			String annot = "/Data/Ref/Ensembl/GrCH37-b86-Structures.txt.gz";

			System.out.println((double) 0 / 100);
			System.exit(-1);
			String annotationfiles = "/Data/Enhancers/Roadmap/dnase-groups.txt";
			annotationfiles = "/Data/Enhancers/ChromHMM/ChromHMMPromotorsEnhancers-groups.txt";
			String bedregions = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/T1D-assoc0.3-COSMO-significantregions-75e7.bed";
			String assocfile = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/T1D-assoc0.3-COSMO-merged-posterior-significantloci-75e7.txt.gz";
//			c.run(annot, bedregions, annotationfiles, assocfile, true);
			c.run(annot, bedregions, annotationfiles, assocfile, true);


			System.out.println();
			bedregions = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/RA-assoc0.3-COSMO-significantregions-75e7.bed";
			assocfile = "/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/RA-assoc0.3-COSMO-merged-posterior-significantloci-75e7.txt.gz";
//			c.run(annot, bedregions, annotationfiles, assocfile, true);
			c.run(annot, bedregions, annotationfiles, assocfile, true);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String annot, String bedregions, String annotationfiles, String assocFile, boolean splitpergroup) throws IOException {

		Annotation geneAnnotation = null;
		if (annot.endsWith(".gtf.gz") || annot.endsWith(".gtf")) {
			geneAnnotation = new GTFAnnotation(annot);
		} else {
			geneAnnotation = new EnsemblStructures(annot);
		}

		//
		Collection<Gene> allgenes = geneAnnotation.getGenes();

		BedFileReader reader = new BedFileReader();
		ArrayList<Feature> regions = reader.readAsList(bedregions);


		BedFileGroup annotation = null;
		if (annotationfiles != null) {
			annotation = new BedFileGroup();
			annotation.loadGroupsFromFile(annotationfiles, regions);
		}


		// region variant1 pval1 posterior1 variant2 pval2 posterior2 variant3 pval3 posterior3
		AssociationFile f = new AssociationFile();

		AssociationResult[][] data = new AssociationResult[regions.size()][];
		AssociationResult[][] crediblesets = new AssociationResult[regions.size()][];

		ApproximateBayesPosterior abp = new ApproximateBayesPosterior();
		ArrayList<Feature> regionsWithCredibleSets = new ArrayList<>();
		double maxPosteriorCredibleSet = 0.9;
		int maxNrVariantsInCredibleSet = 10;
		for (int d = 0; d < regions.size(); d++) {
			boolean hasSet = false;

			ArrayList<AssociationResult> allDatasetData = f.read(assocFile, regions.get(d));
			data[d] = allDatasetData.toArray(new AssociationResult[0]);
			ArrayList<AssociationResult> credibleSet = abp.createCredibleSet(allDatasetData, maxPosteriorCredibleSet);
			crediblesets[d] = credibleSet.toArray(new AssociationResult[0]);
			if (credibleSet.size() <= maxNrVariantsInCredibleSet) {
				hasSet = true;
			}

			if (hasSet) {
				regionsWithCredibleSets.add(regions.get(d));
			}
		}

		// now we have the variants, see if they overlap
		if (annotationfiles != null) {

			for (int i = 0; i < annotation.size(); i++) {
				double sigmaposteriorIndel = 0;
				double sigmaposteriorCoding = 0;
				double sigmaPosteriorOther = 0;

				double fractionIndel = 0;
				double fractionCoding = 0;
				double fractionOther = 0;

				ArrayList<Feature> functionalAnnotation = null;
				if (splitpergroup) {
					functionalAnnotation = annotation.getJointFeaturesForGroup(i, null);
				} else {
					functionalAnnotation = annotation.getJointFeaturesForGroup(null, null);
				}

				for (int d = 0; d < regions.size(); d++) {
					int nrCoding = 0;
					int nrTotal = 0;
					int nrIndel = 0;
					int nrOther = 0;

					Feature region = regions.get(d);
					ArrayList<Gene> genes = new ArrayList<>();
					for (Gene g : allgenes) {
						if (g.overlaps(region)) {
							genes.add(g);
						}
					}

					for (AssociationResult r : data[d]) {
						// check if variant overlaps coding region
						SNPFeature snp = r.getSnp();
						boolean overlaps = false;
						if (snp.isIndel() && snp.overlaps(functionalAnnotation)) {
							sigmaposteriorIndel += r.getPosterior();
							overlaps = true;
							nrIndel++;
						}

						boolean coding = getIsCoding(snp, genes);
						if (coding && (overlaps || snp.overlaps(functionalAnnotation))) {
							sigmaposteriorCoding += r.getPosterior();
							overlaps = true;
							nrCoding++;
						}

						if (!snp.isIndel() && !coding) {
							if (overlaps || snp.overlaps(functionalAnnotation)) {
								sigmaPosteriorOther += r.getPosterior();
								overlaps = true;
								nrOther++;
							}
						}
//						if (overlaps) {
						nrTotal++;
//						}
					}
					if (nrTotal != 0) {
						if (nrCoding > 0) {
							fractionCoding += ((double) nrCoding / nrTotal);
						} else {
							System.out.println("Zero coding variants.. in region " + region.toString());
						}
						if (nrIndel > 0) {
							fractionIndel += ((double) nrIndel / nrTotal);
						} else {
							System.out.println("Zero indel variants.. in region " + region.toString());
						}
						if (nrOther > 0) {
							fractionOther += ((double) nrOther / nrTotal);
						} else {
							System.out.println("Zero other variants.. in region " + region.toString());
						}
					}
				}

				System.out.println(annotation.getGroupName(i) + "\tcoding\tsumposterior:\t" + sigmaposteriorCoding + "\tfraction:\t" + fractionCoding + "\tenrich:\t" + (sigmaposteriorCoding / fractionCoding));
				System.out.println(annotation.getGroupName(i) + "\tindel\tsumposterior:\t" + sigmaposteriorIndel + "\tfraction:\t" + fractionIndel + "\tenrich:\t" + (sigmaposteriorIndel / fractionIndel));
				System.out.println(annotation.getGroupName(i) + "\tother\tsumposterior:\t" + sigmaPosteriorOther + "\tfraction:\t" + fractionOther + "\tenrich:\t" + (sigmaPosteriorOther / fractionOther));
			}
		} else {
			double sigmaposteriorIndel = 0;
			double sigmaposteriorCoding = 0;
			double sigmaPosteriorOther = 0;

			double fractionIndel = 0;
			double fractionCoding = 0;
			double fractionOther = 0;


			for (int d = 0; d < regions.size(); d++) {
				int nrCoding = 0;
				int nrTotal = 0;
				int nrIndel = 0;
				int nrOther = 0;

				Feature region = regions.get(d);
				ArrayList<Gene> genes = new ArrayList<>();
				for (Gene g : allgenes) {
					if (g.overlaps(region)) {
						genes.add(g);
					}
				}


				for (AssociationResult r : data[d]) {
					// check if variant overlaps coding region
					SNPFeature snp = r.getSnp();
					if (snp.isIndel()) {
						sigmaposteriorIndel += r.getPosterior();
						nrIndel++;

					}

					boolean coding = getIsCoding(snp, genes);
					if (coding) {
						sigmaposteriorCoding += r.getPosterior();
						nrCoding++;

					}

					if (!snp.isIndel() && !coding) {
						sigmaPosteriorOther += r.getPosterior();
						nrOther++;
					}
					nrTotal++;
				}
				fractionCoding += ((double) nrCoding / nrTotal);
				fractionIndel += ((double) nrIndel / nrTotal);
				fractionOther += ((double) nrOther / nrTotal);
			}

			System.out.println("coding\tsumposterior:\t" + sigmaposteriorCoding + "\tfraction:\t" + fractionCoding + "\tenrich:\t" + (sigmaposteriorCoding / fractionCoding));
			System.out.println("indel\tsumposterior:\t" + sigmaposteriorIndel + "\tfraction:\t" + fractionIndel + "\tenrich:\t" + (sigmaposteriorIndel / fractionIndel));
			System.out.println("other\tsumposterior:\t" + sigmaPosteriorOther + "\tfraction:\t" + fractionOther + "\tenrich:\t" + (sigmaPosteriorOther / fractionOther));
		}
	}

	private boolean getIsCoding(SNPFeature snp, ArrayList<Gene> genes) {
		for (Gene g : genes) {
			if (g.overlaps(snp)) {
				ArrayList<Transcript> transcripts = g.getTranscripts();
				for (Transcript t : transcripts) {
					ArrayList<Exon> exons = t.getExons();
					ArrayList<UTR> utrs = t.getUTRs();
					boolean overlapsUTR = false;
					if (utrs != null) {
						for (UTR u : utrs) {
							if (u.overlaps(snp)) {
								overlapsUTR = true;
							}
						}
					}

					if (!overlapsUTR) {
						for (Exon e : exons) {
							if (e.overlaps(snp)) {
								return true;
							}
						}
					}
				}
			}
		}
		return false;
	}

}
