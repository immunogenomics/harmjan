package nl.harmjanwestra.finemapping.assoc;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class AssocVariantFilter {
	
	public static void main(String[] args) {

//		String oldFile = "D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\missp\\META-assoc0.3-COSMO-merged-posterior.txt.gz";
//		String newFile = "D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\DataRev2\\normalRerunOfAssoc\\META-assoc0.3-COSMO-merged-posterior.txt.gz";
//		String oldFile = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\multinomial\\outputwithoutcohortcovariate\\META-assoc0.3-COSMO-merged-posterior.txt.gz";
//		String newFile = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\multinomial\\METAPCA-assoc0.3-COSMO-merged-posterior.txt.gz";

//		String blacklist = "D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\DataRev2\\normalRerunOfAssoc\\META-blacklist.txt";
		AssocVariantFilter f = new AssocVariantFilter();
//		try {
//			f.compareVariantLists(oldFile, newFile, blacklist);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		
		
		String[] d = new String[]{
//				"RA",
//				"T1D",
				"META"
		};
		
		
		
		for (String ds : d) {
			
			String blacklist = "d:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\DataRev2\\normalRerunOfAssoc\\" + ds + "-blacklist.txt";
			String[] rafiles = new String[]{
					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\multinomial\\outputwithoutcohortcovariate\\META-assoc0.3-COSMO-merged-posterior.txt.gz"
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\missp\\" + ds + "-assoc0.3-COSMO-merged-posterior.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\conditional\\output\\" + ds + "-assoc0.3-COSMO-gwas-0-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\conditional\\output\\" + ds + "-assoc0.3-COSMO-gwas-1-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\conditional\\output\\" + ds + "-assoc0.3-COSMO-gwas-2-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\conditional\\output\\" + ds + "-assoc0.3-COSMO-gwas-3-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\conditional\\output\\" + ds + "-assoc0.3-COSMO-gwas-4-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\conditional\\output\\" + ds + "-assoc0.3-COSMO-gwas-5-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\conditional\\outputTNFAIPFix\\" + ds + "-assoc0.3-COSMO-gwas-1-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\conditional\\outputTNFAIPFix\\" + ds + "-assoc0.3-COSMO-gwas-2-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\conditional\\outputTNFAIPFix\\" + ds + "-assoc0.3-COSMO-gwas-3-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\conditional\\outputTNFAIPFix\\" + ds + "-assoc0.3-COSMO-gwas-4-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\conditional\\outputTNFAIPFix\\" + ds + "-assoc0.3-COSMO-gwas-5-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\exhaustive\\data\\" + ds + "-assoc0.3-COSMO-chr2-pairwise.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\exhaustive\\data\\" + ds + "-assoc0.3-COSMO-chr6-pairwise.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\exhaustive\\data\\" + ds + "-assoc0.3-COSMO-chr10-pairwise.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\exhaustive\\data\\" + ds + "-assoc0.3-COSMO-chr11-pairwise.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\exhaustive\\data\\" + ds + "-assoc0.3-COSMO-chr12-pairwise.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\exhaustive\\data\\" + ds + "-assoc0.3-COSMO-chr19-pairwise.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\exhaustive\\data\\" + ds + "-assoc0.3-COSMO-chr21-pairwise.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\exhaustive\\data\\" + ds + "-assoc0.3-COSMO-chr6-TNFAIP3-pairwise.txt.gz",
			
			};
			
			String[] raoutfiles = new String[]{
					"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\multinomial\\METAPCA-assoc0.3-COSMO-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\conditional\\output\\" + ds + "-assoc0.3-COSMO-gwas-0-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\conditional\\output\\" + ds + "-assoc0.3-COSMO-gwas-1-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\conditional\\output\\" + ds + "-assoc0.3-COSMO-gwas-2-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\conditional\\output\\" + ds + "-assoc0.3-COSMO-gwas-3-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\conditional\\output\\" + ds + "-assoc0.3-COSMO-gwas-4-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\conditional\\output\\" + ds + "-assoc0.3-COSMO-gwas-5-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\conditional\\outputTNFAIPFix\\" + ds + "-assoc0.3-COSMO-gwas-1-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\conditional\\outputTNFAIPFix\\" + ds + "-assoc0.3-COSMO-gwas-2-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\conditional\\outputTNFAIPFix\\" + ds + "-assoc0.3-COSMO-gwas-3-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\conditional\\outputTNFAIPFix\\" + ds + "-assoc0.3-COSMO-gwas-4-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\conditional\\outputTNFAIPFix\\" + ds + "-assoc0.3-COSMO-gwas-5-merged.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\exhaustive\\data\\" + ds + "-assoc0.3-COSMO-chr2-pairwise.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\exhaustive\\data\\" + ds + "-assoc0.3-COSMO-chr6-pairwise.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\exhaustive\\data\\" + ds + "-assoc0.3-COSMO-chr10-pairwise.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\exhaustive\\data\\" + ds + "-assoc0.3-COSMO-chr11-pairwise.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\exhaustive\\data\\" + ds + "-assoc0.3-COSMO-chr12-pairwise.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\exhaustive\\data\\" + ds + "-assoc0.3-COSMO-chr19-pairwise.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\exhaustive\\data\\" + ds + "-assoc0.3-COSMO-chr21-pairwise.txt.gz",
//					"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\exhaustive\\data\\" + ds + "-assoc0.3-COSMO-chr6-TNFAIP3-pairwise.txt.gz",
			};
			try {
				for (int q = 0; q < rafiles.length; q++) {
					String infile = rafiles[q];
					String outfile = raoutfiles[q];
					if (infile.contains("pairwise")) {
						f.filterForBlackList(infile, outfile, blacklist, true);
					} else {
						f.filterForBlackList(infile, outfile, blacklist, false);
					}
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		}
		
		
	}
	
	public void filterForBlackList(String infile, String outfile, String blacklist, boolean ispairwise) throws IOException {
		System.out.println(infile);
		HashSet<String> variantlist = new HashSet<>();
		TextFile bf = new TextFile(blacklist, TextFile.R);
		variantlist.addAll(bf.readAsArrayList());
		bf.close();
		TextFile in = new TextFile(infile, TextFile.R);
		TextFile out = new TextFile(outfile, TextFile.W);
		out.writeln(in.readLine());
		String[] elems = in.readLineElems(TextFile.tab);
		while (elems != null) {
			if (ispairwise) {
				String var1 = elems[3];
				String var2 = elems[7];
				if (!variantlist.contains(var1) && !variantlist.contains(var2)) {
					out.writeln(Strings.concat(elems, Strings.tab));
				}
			} else {
				String var = elems[3];
				if (!variantlist.contains(var)) {
					out.writeln(Strings.concat(elems, Strings.tab));
				}
			}
			elems = in.readLineElems(TextFile.tab);
		}
		out.close();
		in.close();
	}
	
	public void compareVariantLists(String oldfile, String newfile, String blfile) throws IOException {
		
		
		// find out if new results give us new variants
		AssociationFile f = new AssociationFile();
		
		ArrayList<AssociationResult> oldresults = f.read(oldfile);
		HashSet<String> oldResultVariantNames = new HashSet<String>();
		for (AssociationResult r : oldresults) {
			oldResultVariantNames.add(r.getSnp().toString());
		}
		
		ArrayList<AssociationResult> newresults = f.read(newfile);
		HashSet<String> newResultVariantNames = new HashSet<String>();
		int ctr = 0;
		for (AssociationResult r : newresults) {
			String var = r.getSnp().toString();
			newResultVariantNames.add(var);
			if (!oldResultVariantNames.contains(var)) {
				System.out.println("New variant: " + var + "\t" + r.getLog10Pval() + "\t" + r.getPosterior());
				ctr++;
			}
		}
		
		
		TextFile blacklist = new TextFile(blfile, TextFile.W);
		System.out.println(ctr + " / " + newresults.size());
		System.out.println();
		int ctr2 = 0;
		for (AssociationResult r : oldresults) {
			String var = r.getSnp().toString();
			if (!newResultVariantNames.contains(var)) {
				System.out.println("Old variant: " + var + "\t" + r.getLog10Pval() + "\t" + r.getPosterior());
				blacklist.writeln(var);
				ctr2++;
			}
		}
		blacklist.close();
		System.out.println(ctr2 + " / " + oldresults.size());
		
		
	}
	
}
