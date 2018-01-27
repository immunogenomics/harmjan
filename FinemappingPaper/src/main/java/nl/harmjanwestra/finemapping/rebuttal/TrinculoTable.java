package nl.harmjanwestra.finemapping.rebuttal;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.association.AssociationResultPValuePair;
import nl.harmjanwestra.utilities.association.approximatebayesposterior.ApproximateBayesPosterior;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class TrinculoTable {
	
	public static void main(String[] args) {
		TrinculoTable t = new TrinculoTable();
		String prefix = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\trinculo\\input\\META-COSMO-0.3-chr";
		String regionfile = "c:/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
		String origassocfile = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\META-assoc0.3-COSMO-merged-posterior.txt.gz";
		
		String pcafile = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\multinomial\\METAPCA-assoc0.3-COSMO-merged.txt.gz";
		try {
			String output = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\trinculo\\output\\table-multinomcomp.txt";
			t.compareToMultiNom(prefix, regionfile, origassocfile, output);
			output = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\trinculo\\output\\table-pcacomp.txt";
			t.compareToPCA(regionfile, origassocfile, pcafile, output);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void compareToPCA(String regionfile, String origassocfile, String pcaassocfile, String output) throws IOException {
		BedFileReader b = new BedFileReader();
		ArrayList<Feature> regions = b.readAsList(regionfile);
		ArrayList<TrinculoAssoc> allassocs = new ArrayList<>();
		
		AssociationFile f = new AssociationFile();
		ArrayList<AssociationResult> origAssocs = f.read(origassocfile);
		
		
		TextFile out = new TextFile(output, TextFile.W);
		FinemapTable t = new FinemapTable();
		ApproximateBayesPosterior bp = new ApproximateBayesPosterior();
		ArrayList<AssociationResult> pcaAssocs = f.read(pcaassocfile);
		
		String header = "Region\t" +
				"#AssociationsTrinculo\t" +
				"PosteriorCorrelationPearson\t" +
				"PosteriorCorrelationSpearman\t" +
				"CredibleSetSize\t" +
				"CredibleSetSizePCA\t" +
				"TopVariant\t" +
				"TopVariantP\t" +
				"TopVariantPosterior\t" +
				"TopVariantSignificant\t" +
				"TopPCAVariant\t" +
				"TopPCAVariantP\t" +
				"TopPCAVariantPosterior\t" +
				"TopPCAVariantSignificant";
		out.writeln(header);
		
		for (Feature region : regions) {
			// META-COSMO-0.3-chr1-region-Chr1_113863087-114527968.out.assoc.finemap.gz
			ArrayList<AssociationResult> regionassoc = t.getRegionAssocs(pcaAssocs, region);
			bp.calculatePosterior(regionassoc);
			
			ArrayList<AssociationResultPValuePair> pairs = new ArrayList<AssociationResultPValuePair>();
			for (AssociationResult r : regionassoc) {
				pairs.add(new AssociationResultPValuePair(r, r.getPosterior(), false));
			}
			Collections.sort(pairs);
			ArrayList<AssociationResult> tmp = new ArrayList<>();
			for (AssociationResultPValuePair p : pairs) {
				tmp.add(p.getAssociationResult());
			}
			regionassoc = tmp;
			
			ArrayList<AssociationResult> origregionassoc = t.getRegionAssocs(origAssocs, region);
			HashMap<String, AssociationResult> assocmap = new HashMap<String, AssociationResult>();
			for (AssociationResult r : origregionassoc) {
				assocmap.put(r.getSnp().toString(), r);
			}
			
			bp.calculatePosterior(origregionassoc);
			ArrayList<AssociationResult> credibleSetHJW = bp.createCredibleSet(origregionassoc, 0.95);
			ArrayList<AssociationResult> credibleSetPCA = bp.createCredibleSet(regionassoc, 0.95);
			
			
			double[] posteriorTrinculo = new double[regionassoc.size()];
			double[] posteriorHJW = new double[regionassoc.size()];
			
			
			double sum = 0;
			AssociationResult topassoc = null;
			double maxp = 0;
			for (int r = 0; r < regionassoc.size(); r++) {
				AssociationResult assoc = regionassoc.get(r);
				
				AssociationResult origassoc = assocmap.get(assoc.getSnp().toString());
				if (origassoc != null) {
					posteriorTrinculo[r] = assoc.getPosterior();
					posteriorHJW[r] = origassoc.getPosterior();
				}
				
				
				// sum the posterior until we hit 95%
//				if (sum < 0.95) {
//					crediblesetTrinculo.add(assoc);
//					sum += regionassoc.get(r).pcausal;
//				}
				
				if (assoc.getLog10Pval() > maxp) {
					maxp = assoc.getLog10Pval();
					topassoc = assoc;
				}
				
			}
			
			AssociationResult tophjwassoc = null;
			maxp = 0;
			for (AssociationResult r : origregionassoc) {
				if (r.getLog10Pval() > maxp) {
					maxp = r.getLog10Pval();
					tophjwassoc = r;
				}
			}
			
			
			SpearmansCorrelation corrs = new SpearmansCorrelation();
			PearsonsCorrelation corrp = new PearsonsCorrelation();
			
			double spearmanPost = corrs.correlation(posteriorTrinculo, posteriorHJW);
			double pearsonPost = corrp.correlation(posteriorTrinculo, posteriorHJW);
			
			
			double log10p = -Math.log10(7.5E-7);
			
			String ln = region.toString()
					+ "\t" + regionassoc.size()
					+ "\t" + pearsonPost
					+ "\t" + spearmanPost
					+ "\t" + credibleSetHJW.size()
					+ "\t" + credibleSetPCA.size()
					+ "\t" + tophjwassoc.getSnp().toString()
					+ "\t" + tophjwassoc.getLog10Pval()
					+ "\t" + tophjwassoc.getPosterior()
					+ "\t" + (tophjwassoc.getLog10Pval() > log10p)
					+ "\t" + topassoc.getSnp().toString()
					+ "\t" + topassoc.getLog10Pval()
					+ "\t" + topassoc.getPosterior()
					+ "\t" + (topassoc.getLog10Pval() > log10p);
			
			
			out.writeln(ln);
			
		}
		out.close();
		
		
	}
	
	public void compareToMultiNom(String prefix, String regionfile, String origassocfile, String output) throws IOException {
		
		BedFileReader b = new BedFileReader();
		ArrayList<Feature> regions = b.readAsList(regionfile);
		ArrayList<TrinculoAssoc> allassocs = new ArrayList<>();
		
		AssociationFile f = new AssociationFile();
		ArrayList<AssociationResult> origAssocs = f.read(origassocfile);
		
		TextFile out = new TextFile(output, TextFile.W);
		FinemapTable t = new FinemapTable();
		ApproximateBayesPosterior bp = new ApproximateBayesPosterior();
		
		String header = "Region\t" +
				"#AssociationsTrinculo\t" +
				"DiseaseORCorrelationPearson\t" +
				"DiseaseORCorrelationSpearman\t" +
				"PosteriorCorrelationPearson\t" +
				"PosteriorCorrelationSpearman\t" +
				"CredibleSetSizeBinomial\t" +
				"CredibleSetSizeMultinomial\t" +
				"TopBinomialVariant\t" +
				"TopBinomialVariantP\t" +
				"TopBinomialVariantPosterior\t" +
				"TopBinomialVariantSignificant\t" +
				"TopMultinomialVariant\t" +
				"TopMultinomialVariantBF\t" +
				"TopMultinomialVariantPosterior\t" +
				"TopMultinomialVariantSignificant";
		out.writeln(header);
		
		
		for (Feature region : regions) {
			// META-COSMO-0.3-chr1-region-Chr1_113863087-114527968.out.assoc.finemap.gz
			ArrayList<TrinculoAssoc> regionassoc = parseAssoc(prefix + region.getChromosome().getNumber() + "-region-" + region.toString() + ".out.assoc.finemap.gz");
			Collections.sort(regionassoc);
			
			ArrayList<AssociationResult> origregionassoc = t.getRegionAssocs(origAssocs, region);
			HashMap<String, AssociationResult> assocmap = new HashMap<String, AssociationResult>();
			for (AssociationResult r : origregionassoc) {
				assocmap.put(r.getSnp().toString(), r);
			}
			
			bp.calculatePosterior(origregionassoc);
			ArrayList<AssociationResult> credibleSetHJW = bp.createCredibleSet(origregionassoc, 0.95);
			
			
			double[] xOR = new double[regionassoc.size()];
			double[] yOR = new double[regionassoc.size()];
			
			double[] posteriorTrinculo = new double[regionassoc.size()];
			double[] posteriorHJW = new double[regionassoc.size()];
			
			ArrayList<TrinculoAssoc> crediblesetTrinculo = new ArrayList<>();
			double sum = 0;
			TrinculoAssoc topassoc = null;
			double maxp = 0;
			for (int r = 0; r < regionassoc.size(); r++) {
				TrinculoAssoc assoc = regionassoc.get(r);
				xOR[r] = assoc.ord1;
				yOR[r] = assoc.ord2;
				
				AssociationResult origassoc = assocmap.get(assoc.rsid);
				if (origassoc != null) {
					posteriorTrinculo[r] = assoc.pcausal;
					posteriorHJW[r] = origassoc.getPosterior();
				}
				
				
				// sum the posterior until we hit 95%
				if (sum < 0.95) {
					crediblesetTrinculo.add(assoc);
					sum += regionassoc.get(r).pcausal;
				}
				
				if (assoc.logbf > maxp) {
					maxp = assoc.logbf;
					topassoc = assoc;
				}
				
			}
			
			AssociationResult tophjwassoc = null;
			maxp = 0;
			for (AssociationResult r : origregionassoc) {
				if (r.getLog10Pval() > maxp) {
					maxp = r.getLog10Pval();
					tophjwassoc = r;
				}
			}
			
			
			SpearmansCorrelation corrs = new SpearmansCorrelation();
			PearsonsCorrelation corrp = new PearsonsCorrelation();
			double spearmanOR = corrs.correlation(xOR, yOR);
			double pearsonOR = corrp.correlation(xOR, yOR);
			double spearmanPost = corrs.correlation(posteriorTrinculo, posteriorHJW);
			double pearsonPost = corrp.correlation(posteriorTrinculo, posteriorHJW);
			
			
			double log10p = -Math.log10(7.5E-7);
			
			String ln = region.toString()
					+ "\t" + regionassoc.size()
					+ "\t" + pearsonOR
					+ "\t" + spearmanOR
					+ "\t" + pearsonPost
					+ "\t" + spearmanPost
					+ "\t" + credibleSetHJW.size()
					+ "\t" + crediblesetTrinculo.size()
					+ "\t" + tophjwassoc.getSnp().toString()
					+ "\t" + tophjwassoc.getLog10Pval()
					+ "\t" + tophjwassoc.getPosterior()
					+ "\t" + (tophjwassoc.getLog10Pval() > log10p)
					+ "\t" + topassoc.rsid
					+ "\t" + topassoc.logbf
					+ "\t" + topassoc.pcausal
					+ "\t" + (topassoc.logbf > log10p);
			
			
			out.writeln(ln);
			
			allassocs.addAll(regionassoc);
		}
		out.close();
		
		
	}
	
	public ArrayList<TrinculoAssoc> parseAssoc(String assoc) throws IOException {
		System.out.println("Parsing " + assoc);
		TextFile tf = new TextFile(assoc, TextFile.R);
		
		String[] header = tf.readLineElems(TextFile.tab);
		
		ArrayList<TrinculoAssoc> assocs = new ArrayList<>();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			
			TrinculoAssoc a = new TrinculoAssoc();
			a.rsid = elems[0];
			a.ord1 = Double.parseDouble(elems[1]);
			a.ord2 = Double.parseDouble(elems[2]);
			a.lp0 = Double.parseDouble(elems[3]);
			a.lp1 = Double.parseDouble(elems[4]);
			a.logbf = Double.parseDouble(elems[5]);
			a.pcausal = Double.parseDouble(elems[6]);
			assocs.add(a);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		return assocs;
	}
	
	public class TrinculoAssoc implements Comparable<TrinculoAssoc> {
		// RSID    OR_Disease1     OR_Disease2     L+P0    L+P1    logBF   P_CAUSAL
		String rsid;
		double ord1;
		double ord2;
		double lp0;
		double lp1;
		double logbf;
		double pcausal;
		
		@Override
		public boolean equals(Object o) {
			if (this == o) return true;
			if (o == null || getClass() != o.getClass()) return false;
			
			TrinculoAssoc that = (TrinculoAssoc) o;
			
			if (Double.compare(that.ord1, ord1) != 0) return false;
			if (Double.compare(that.ord2, ord2) != 0) return false;
			if (Double.compare(that.lp0, lp0) != 0) return false;
			if (Double.compare(that.lp1, lp1) != 0) return false;
			if (Double.compare(that.logbf, logbf) != 0) return false;
			if (Double.compare(that.pcausal, pcausal) != 0) return false;
			return rsid != null ? rsid.equals(that.rsid) : that.rsid == null;
		}
		
		@Override
		public int hashCode() {
			int result;
			long temp;
			result = rsid != null ? rsid.hashCode() : 0;
			temp = Double.doubleToLongBits(ord1);
			result = 31 * result + (int) (temp ^ (temp >>> 32));
			temp = Double.doubleToLongBits(ord2);
			result = 31 * result + (int) (temp ^ (temp >>> 32));
			temp = Double.doubleToLongBits(lp0);
			result = 31 * result + (int) (temp ^ (temp >>> 32));
			temp = Double.doubleToLongBits(lp1);
			result = 31 * result + (int) (temp ^ (temp >>> 32));
			temp = Double.doubleToLongBits(logbf);
			result = 31 * result + (int) (temp ^ (temp >>> 32));
			temp = Double.doubleToLongBits(pcausal);
			result = 31 * result + (int) (temp ^ (temp >>> 32));
			return result;
		}
		
		@Override
		public int compareTo(TrinculoAssoc o) {
			if (this.equals(o)) {
				return 0;
			} else if (this.pcausal > o.pcausal) {
				return -1;
			} else {
				return 1;
			}
			
		}
	}
	
}
