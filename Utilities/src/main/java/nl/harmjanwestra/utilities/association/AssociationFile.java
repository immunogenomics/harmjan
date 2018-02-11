package nl.harmjanwestra.utilities.association;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by hwestra on 10/27/15.
 */
public class AssociationFile {
	
	
	private String model = null;
	
	public ArrayList<AssociationResult> readTabFile(String pvaluefile, Feature region) throws IOException {
		System.out.println("Reading tab path: " + pvaluefile);
		HashSet<String> variantHash = new HashSet<String>();
		TextFile textfile = new TextFile(pvaluefile, TextFile.R);
		
		ArrayList<AssociationResult> output = new ArrayList<AssociationResult>();
		
		// skip header
		String[] headerelems = textfile.readLineElems(TextFile.tab);
		
		// Marker	Chr	Position	PValue	OR(MinAllele)	LowerOR	UpperOR	Alleles(Maj>Min)
		int markercol = 0;
		int chrcol = 0;
		int poscol = 0;
		int pvalcol = 0;
		int orcol = 0;
		int lowerorcol = 0;
		int upperorcol = 0;
		int allelescol = 0;
		
		for (int c = 0; c < headerelems.length; c++) {
			String elem = headerelems[c].toLowerCase();
			if (elem.equals("marker")) {
				markercol = c;
			} else if (elem.equals("chr")) {
				chrcol = c;
			} else if (elem.equals("position")) {
				poscol = c;
			} else if (elem.equals("pvalue")) {
				pvalcol = c;
			} else if (elem.equals("OR(MinAllele)".toLowerCase())) {
				orcol = c;
			} else if (elem.equals("LowerOR".toLowerCase())) {
				lowerorcol = c;
			} else if (elem.equals("UpperOR".toLowerCase())) {
				upperorcol = c;
			} else if (elem.equals("Alleles(Maj>Min)".toLowerCase())) {
				allelescol = c;
			}
		}
		
		String[] elems = textfile.readLineElems(TextFile.tab);
		int pvalctr = 0;
		while (elems != null) {
			// Marker	Chr	Position	PValue	Odds Ratio
			
			Chromosome chr = Chromosome.parseChr(elems[1]);
			if (region.getChromosome().equals(chr)) {
				String variant = elems[chrcol] + "_" + elems[poscol] + "_" + elems[markercol];
				SNPFeature f2 = new SNPFeature();
				f2.setChromosome(chr);
				f2.setStart(Integer.parseInt(elems[poscol]));
				f2.setStop(Integer.parseInt(elems[poscol]));
				if (f2.overlaps(region)) {
					try {
						
						String alleles = elems[allelescol];
						String[] allelesArr = alleles.split(">");
						f2.setAlleles(allelesArr);
						if (allelesArr.length < 2) {
							f2.setMinorAllele(allelesArr[0]);
						} else {
							f2.setMinorAllele(allelesArr[1]);
						}
						
						Double or = Double.parseDouble(elems[orcol]);
						
						Double pval = 1d;
						try {
							pval = Double.parseDouble(elems[pvalcol]);
						} catch (NumberFormatException e) {
						
						}
						variantHash.add(variant);
						
						f2.setName(elems[markercol]);
						AssociationResult r = new AssociationResult();
						r.setSnp(f2);
						r.setPval(pval);
						double[][] ors = new double[1][1];
						ors[0][0] = or;
						r.setOR(ors);
						
						output.add(r);
					} catch (NumberFormatException e) {
						e.printStackTrace();
					}
					pvalctr++;
				}
				
			}
			elems = textfile.readLineElems(TextFile.tab);
		}
		textfile.close();
		
		System.out.println(pvalctr + " pvalues for " + pvalctr + " positions from path: " + pvaluefile);
		
		return output;
		
	}
	
	public String getModel() {
		return model;
	}
	
	public ArrayList<AssociationResult> read(String file) throws IOException {
		Feature region = null;
		return read(file, region);
	}
	
	
	public ArrayList<AssociationResult> read(String file, Feature region) throws IOException {
		
		if (file.endsWith("tab") || file.endsWith("tab.gz")) {
			return readTabFile(file, region);
		}
		
		if (file.endsWith("OKADA.gz")) {
			return readOkada(file, region);
		}
		
		if (region == null) {
			ArrayList<Feature> f = null;
			return read(file, f);
		}
		ArrayList<Feature> regions = new ArrayList<>();
		regions.add(region);
		return read(file, regions);
		
		
	}
	
	private ArrayList<AssociationResult> readOkada(String pvaluefile, Feature region) throws IOException {
		System.out.println("Reading tab path: " + pvaluefile);
		HashSet<String> variantHash = new HashSet<String>();
		TextFile textfile = new TextFile(pvaluefile, TextFile.R);
		
		ArrayList<AssociationResult> output = new ArrayList<AssociationResult>();
		
		// skip header
		String[] headerelems = textfile.readLineElems(TextFile.tab);
		
		// Marker	Chr	Position	PValue	OR(MinAllele)	LowerOR	UpperOR	Alleles(Maj>Min)
		int markercol = 0;
		int chrcol = 0;
		int poscol = 0;
		int pvalcol = 0;
		int orcol = 0;
		int lowerorcol = 0;
		int upperorcol = 0;
		int a1col = 0;
		int a2col = 0;
		
		for (int c = 0; c < headerelems.length; c++) {
			String elem = headerelems[c].toLowerCase();
			if (elem.equals("SNPID".toLowerCase())) {
				markercol = c;
			} else if (elem.equals("Chr".toLowerCase())) {
				chrcol = c;
			} else if (elem.equals("Position(hg19)".toLowerCase())) {
				poscol = c;
			} else if (elem.equals("P-val".toLowerCase())) {
				pvalcol = c;
			} else if (elem.equals("OR(A1)".toLowerCase())) {
				orcol = c;
			} else if (elem.equals("OR_95%CIlow".toLowerCase())) {
				lowerorcol = c;
			} else if (elem.equals("OR_95%CIup".toLowerCase())) {
				upperorcol = c;
			} else if (elem.equals("A1".toLowerCase())) {
				a1col = c;
			} else if (elem.equals("A2".toLowerCase())) {
				a2col = c;
			}
		}
		
		String[] elems = textfile.readLineElems(TextFile.tab);
		int pvalctr = 0;
		while (elems != null) {
			// Marker	Chr	Position	PValue	Odds Ratio
			
			Chromosome chr = Chromosome.parseChr(elems[1]);
			if (region.getChromosome().equals(chr)) {
				String variant = elems[chrcol] + "_" + elems[poscol] + "_" + elems[markercol];
				SNPFeature f2 = new SNPFeature();
				f2.setChromosome(chr);
				f2.setStart(Integer.parseInt(elems[poscol]));
				f2.setStop(Integer.parseInt(elems[poscol]));
				if (f2.overlaps(region)) {
					try {
						
						String a1 = elems[a1col];
						String a2 = elems[a2col];
						String[] allelesArr = new String[]{a1, a2};
						f2.setAlleles(allelesArr);
						if (allelesArr.length < 2) {
							f2.setMinorAllele(allelesArr[0]);
						} else {
							f2.setMinorAllele(allelesArr[1]);
						}
						
						Double or = Double.parseDouble(elems[orcol]);
						
						Double pval = 1d;
						try {
							pval = Double.parseDouble(elems[pvalcol]);
						} catch (NumberFormatException e) {
						
						}
						variantHash.add(variant);
						
						f2.setName(elems[markercol]);
						AssociationResult r = new AssociationResult();
						r.setSnp(f2);
						r.setPval(pval);
						double[][] ors = new double[1][1];
						ors[0][0] = or;
						r.setOR(ors);
						
						output.add(r);
					} catch (NumberFormatException e) {
						e.printStackTrace();
					}
					pvalctr++;
				}
				
			}
			elems = textfile.readLineElems(TextFile.tab);
		}
		textfile.close();
		
		System.out.println(pvalctr + " pvalues for " + pvalctr + " positions from path: " + pvaluefile);
		
		return output;
	}
	
	public ArrayList<AssociationResult> read(String file, ArrayList<Feature> regions) throws IOException {
		
		System.out.println("Reading assoc path: " + file);
		TextFile tf = new TextFile(file, TextFile.R);
		String ln = tf.readLine();
		
		model = null;
		
		// Chr     Pos     Id      CombinedId      N       MAF     DevianceNull    DevianceGeno    Df      Beta(Genotype)  SE(Genotype)    OR      OR-Hi   OR-Lo   Pval    -Log10(pval)
		int chrcol = -1;
		int poscol = -1;
		int idcol = -1;
		int allelesCol = -1;
		int minorAlleleCol = -1;
		int combinedIdCol = -1;
		int ncol = -1;
		int crcasescol = -1;
		int crcontrolscol = -1;
		int missingnesspcol = -1;
		int mafcol = -1;
		int afcasescol = -1;
		int afcontrolscol = -1;
		int impqualscorecol = -1;
		int deviancenullcol = -1;
		int deviancegenocol = -1;
		int dfnullcol = -1;
		int dfaltcol = -1;
		int dfcol = -1;
		
		int betacol = -1;
		int secol = -1;
		int orcol = -1;
		int orhicol = -1;
		int orlocol = -1;
		int pvalcol = -1;
		int log10pvalcol = -1;
		int posteriorcol = -1;
		int bfcol = -1;
		int regioncol = -1;
		int hwepcol = -1;
		
		
		ArrayList<AssociationResult> results = new ArrayList<AssociationResult>();
		int nr = 0;
		while (ln != null) {
			if (ln.startsWith("#Chromosome") || ln.startsWith("Chr\tPos")) {
// skip header
//				System.out.println("Found header");
				
				String[] elems = Strings.tab.split(ln);
				for (int i = 0; i < elems.length; i++) {
					String e = elems[i];
					if (e.equals("#Chromosome") || e.equals("Chr")) {
						chrcol = i;
					} else if (e.equals("Pos")) {
						poscol = i;
					} else if (e.equals("Id")) {
						idcol = i;
					} else if (e.equals("CombinedId")) {
						combinedIdCol = i;
					} else if (e.equals("Alleles")) {
						allelesCol = i;
					} else if (e.equals("MinorAllele")) {
						minorAlleleCol = i;
					} else if (e.equals("N")) {
						ncol = i;
					} else if (e.equals("CallRateCases")) {
						crcasescol = i;
					} else if (e.equals("CallRateControls")) {
						crcontrolscol = i;
					} else if (e.equals("DifferentialMissingnessP")) {
						missingnesspcol = i;
					} else if (e.equals("MAF")) {
						mafcol = i;
					} else if (e.equals("AFCases")) {
						afcasescol = i;
					} else if (e.equals("AFControls")) {
						afcontrolscol = i;
					} else if (e.equals("HWEP")) {
						hwepcol = i;
					} else if (e.equals("DevianceNull")) {
						deviancenullcol = i;
					} else if (e.equals("DevianceGeno")) {
						deviancegenocol = i;
					} else if (e.equals("DfNull")) {
						dfnullcol = i;
					} else if (e.equals("DfAlt")) {
						dfaltcol = i;
					} else if (e.equals("DiffDf")) {
						dfcol = i;
					} else if (e.equals("Beta(Genotype)")) {
						betacol = i;
					} else if (e.equals("ImputationQualScore")) {
						impqualscorecol = i;
					} else if (e.equals("SE(Genotype)")) {
						secol = i;
					} else if (e.equals("OR")) {
						orcol = i;
					} else if (e.equals("OR-Hi")) {
						orhicol = i;
					} else if (e.equals("OR-Lo")) {
						orlocol = i;
					} else if (e.equals("Pval")) {
						pvalcol = i;
					} else if (e.equals("-Log10(pval)")) {
						log10pvalcol = i;
					} else if (e.equals("Region")) {
						regioncol = i;
					} else if (e.equals("BF")) {
						bfcol = i;
					} else if (e.equals("Posterior")) {
						posteriorcol = i;
					}
				}
			} else if (ln.startsWith("#")) {
				// get the model
				model = ln;
			} else {
				String[] elems = Strings.tab.split(ln);
				if (elems.length > 4) {
					// VariantID	N	MAF	DevianceNull	DfNull	DevianceGeno	DfAlt	Beta(Genotype)	SE(Genotype)	OR	OR-Hi	OR-Lo	Pval	-Log10(pval)
					Chromosome chr = Chromosome.NA;
					int pos = -1;
					String id = null;
					int n = 0;
					double maf = 0d;
					double crcases = 1d;
					double crcontrols = 1d;
					double missingnessp = 1;
					double hwep = 0d;
					double deviancenull = 0d;
					double deviancegeno = 0d;
					
					double afcases = 0d;
					double afcontrols = 0d;
					int dfnull = 0;
					int dfalt = 0;
					int df = 0;
					double[][] beta = null;
					double[][] se = null;
					
					double pval = 1d;
					
					double bf = 0d;
					double posterior = 0d;
					Feature assocregion = null;
					double impqualscore = Double.NaN;
					String[] alleles = null;
					String minorAllele = null;
					
					if (chrcol != -1) {
						chr = Chromosome.parseChr(elems[chrcol]);
					}
					if (poscol != -1) {
						try {
							pos = Integer.parseInt(elems[poscol]);
						} catch (NumberFormatException e) {
						
						}
					}
					if (idcol != -1) {
						id = elems[idcol];
					}
					
					if (allelesCol != -1) {
						String alleleStr = elems[allelesCol];
						String[] alleletmp = alleleStr.split(",");
						alleles = new String[alleletmp.length];
						for (int i = 0; i < alleletmp.length; i++) {
							alleles[i] = new String(alleletmp[i]).intern();
						}
					}
					if (minorAlleleCol != -1) {
						minorAllele = new String(elems[minorAlleleCol]).intern();
					}
					if (ncol != -1) {
						try {
							n = Integer.parseInt(elems[ncol]);
						} catch (NumberFormatException e) {
						
						}
					}
					if (mafcol != -1) {
						try {
							maf = Double.parseDouble(elems[mafcol]);
						} catch (NumberFormatException e) {
						
						}
					}
					
					if (crcasescol != -1) {
						try {
							crcases = Double.parseDouble(elems[crcasescol]);
						} catch (NumberFormatException e) {
						
						}
					}
					
					if (crcontrolscol != -1) {
						try {
							crcontrols = Double.parseDouble(elems[crcontrolscol]);
						} catch (NumberFormatException e) {
						
						}
					}
					
					if (missingnesspcol != -1) {
						try {
							missingnessp = Double.parseDouble(elems[missingnesspcol]);
						} catch (NumberFormatException e) {
						
						}
					}
					
					if (afcasescol != -1) {
						try {
							afcases = Double.parseDouble(elems[afcasescol]);
						} catch (NumberFormatException e) {
						
						}
					}
					
					if (afcontrolscol != -1) {
						try {
							afcontrols = Double.parseDouble(elems[afcontrolscol]);
						} catch (NumberFormatException e) {
						
						}
					}
					
					if (hwepcol != -1) {
						try {
							hwep = Double.parseDouble(elems[hwepcol]);
						} catch (NumberFormatException e) {
						
						}
					}
					if (deviancenullcol != -1) {
						try {
							deviancenull = Double.parseDouble(elems[deviancenullcol]);
						} catch (NumberFormatException e) {
						
						}
					}
					if (deviancegenocol != -1) {
						try {
							deviancegeno = Double.parseDouble(elems[deviancegenocol]);
						} catch (NumberFormatException e) {
						
						}
					}
					if (dfnullcol != -1) {
						try {
							dfnull = Integer.parseInt(elems[dfnullcol]);
						} catch (NumberFormatException e) {
						
						}
						
					}
					if (dfaltcol != -1) {
						try {
							dfalt = Integer.parseInt(elems[dfaltcol]);
						} catch (NumberFormatException e) {
						
						}
						
					}
					
					if (dfcol != -1) {
						try {
							df = Integer.parseInt(elems[dfcol]);
						} catch (NumberFormatException e) {
						
						}
						
					}
					
					if (betacol != -1) {
						String betaStr = elems[betacol];
						String[] betaStrElems = betaStr.split(";");
						
						beta = new double[betaStrElems.length][];
						for (int i = 0; i < betaStrElems.length; i++) {
							String[] betaAlleleElems = betaStrElems[i].split(",");
							beta[i] = new double[betaAlleleElems.length];
							for (int j = 0; j < betaAlleleElems.length; j++) {
								try {
									beta[i][j] = Double.parseDouble(betaAlleleElems[j]);
								} catch (NumberFormatException e) {
								
								}
							}
						}
					}
					
					if (secol != -1) {
						String seStr = elems[secol];
						String[] seStrElems = seStr.split(";");
						se = new double[seStrElems.length][];
						for (int i = 0; i < seStrElems.length; i++) {
							String[] seAlleleElems = seStrElems[i].split(",");
							se[i] = new double[seAlleleElems.length];
							for (int j = 0; j < seAlleleElems.length; j++) {
								try {
									se[i][j] = Double.parseDouble(seAlleleElems[j]);
								} catch (NumberFormatException e) {
								
								}
							}
						}
					}
					
					
					if (pvalcol != -1 && pvalcol < elems.length) {
						try {
							pval = Double.parseDouble(elems[pvalcol]);
						} catch (NumberFormatException e) {
						
						}
					}
					
					
					if (impqualscorecol != -1) {
						try {
							impqualscore = Double.parseDouble(elems[impqualscorecol]);
						} catch (NumberFormatException e) {
						
						}
					}
					
					if (regioncol != -1) {
						String regionStr = elems[regioncol];
						assocregion = Feature.parseFeature(regionStr);
					}
					
					if (bfcol != -1) {
						try {
							bf = Double.parseDouble(elems[bfcol]);
						} catch (NumberFormatException e) {
						
						}
					}
					if (posteriorcol != -1) {
						try {
							posterior = Double.parseDouble(elems[posteriorcol]);
						} catch (NumberFormatException e) {
						
						}
					}
					
					
					SNPFeature snp = new SNPFeature(chr, pos, pos);
					snp.setName(id);
					snp.setImputationQualityScore(impqualscore);
					snp.setAlleles(alleles);
					snp.setMinorAllele(minorAllele);
					snp.setMissingnessP(missingnessp);
					snp.setCrCases(crcases);
					snp.setCrControls(crcontrols);
					snp.setMaf(maf);
					snp.setAFCases(afcases);
					snp.setAFControls(afcontrols);
					
					if (regions == null || snp.overlaps(regions)) {
						AssociationResult result = new AssociationResult();
						result.setSnp(snp);
						
						result.setN(n);
						
						result.setDevianceNull(deviancenull);
						result.setDevianceGeno(deviancegeno);
						result.setDfnull(dfnull);
						result.setDfalt(dfalt);
						result.setDf(df);
						
						result.setBeta(beta);
						result.setSe(se);
						result.setPval(pval);
						result.setBf(bf);
						snp.setHwep(hwep);
						result.setPosterior(posterior);
						result.setRegion(assocregion);
						
						results.add(result);
						
					}
				}
			}
			
			ln = tf.readLine();
		}
		tf.close();
		
		return results;
	}
	
	public String getHeader() {
		
		String str = "";
		if (model != null) {
			str = model + "\n";
		}
		str += "#Chromosome" +
				"\tPos" +
				"\tId" +
				"\tCombinedId" +
				"\tAlleles" +
				"\tMinorAllele" +
				"\tImputationQualScore" +
				"\tN" +
				"\tCallRateCases" +
				"\tCallRateControls" +
				"\tDifferentialMissingnessP" +
				"\tMAF" +
				"\tAFCases" +
				"\tAFControls" +
				"\tHWEP" +
				"\tDevianceNull" +
				"\tDevianceGeno" +
				"\tDfNull" +
				"\tDfAlt" +
				"\tDiffDf" +
				"\tBeta(Genotype)" +
				"\tSE(Genotype)" +
				"\tOR" +
				"\tOR-Hi" +
				"\tOR-Lo" +
				"\tPval" +
				"\tLog10(p)";
		
		return str;
	}
	
}
