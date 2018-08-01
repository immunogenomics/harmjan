package nl.harmjanwestra.finemapping.comparisons;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class CompareDatasetPosteriors {
	public static void main(String[] args) {
		
		String[] files = new String[]{
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\multinomial\\outputwithcohortcovariate\\META-ORIG-assoc0.3-COSMO-merged-posterior.txt.gz",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\multinomial\\outputwithcohortcovariate\\META-assoc0.3-COSMO-merged-posterior.txt.gz",
				"C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\multinomial\\outputwithcohortcovariate\\MULTI-assoc0.3-COSMO-merged-posterior.txt.gz"
		};
		String[] filedesc = new String[]{
				"META-Covariates",
				"META-Multinomial"
		};
		String out = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4\\multinomial\\outputwithcohortcovariate\\posteriorcomp.txt";
		CompareDatasetPosteriors c = new CompareDatasetPosteriors();
		String regionAnnot = "C:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\AllLoci-GenesPerLocus.txt";
		try {
			c.run(files, filedesc, regionAnnot, out);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String[] files, String[] filedesc, String regionAnnot, String out) throws IOException {
		
		HashMap<String, String> regionToGenes = new HashMap<>();
		if (regionAnnot != null) {
			TextFile tf = new TextFile(regionAnnot, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				if (elems.length >= 2) {
					regionToGenes.put(elems[0], elems[1]);
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
		}
		
		AssociationFile f = new AssociationFile();
		
		HashMap<String, Integer> snpIndex = new HashMap<String, Integer>();
		
		ArrayList<ArrayList<AssociationResult>> results = new ArrayList<>();
		ArrayList<String> snps = new ArrayList<>();
		ArrayList<String> regions = new ArrayList<>();
		
		int ctr = 0;
		for (int i = 0; i < files.length; i++) {
			ArrayList<AssociationResult> r1 = f.read(files[i]);
			results.add(r1);
			for (AssociationResult r : r1) {
				String snp = r.getSnp().toString();
				if (!snpIndex.containsKey(snp)) {
					snpIndex.put(snp, ctr);
					regions.add(r.getRegion().toString());
					snps.add(snp);
					ctr++;
				}
			}
		}
		
		double[][] mat = new double[snpIndex.size()][files.length];
		for (int i = 0; i < mat.length; i++) {
			for (int j = 0; j < mat[i].length; j++) {
				mat[i][j] = Double.NaN;
			}
		}
		
		for (int i = 0; i < files.length; i++) {
			ArrayList<AssociationResult> r1 = results.get(i);
			for (AssociationResult r : r1) {
				String snp = r.getSnp().toString();
				Integer index = snpIndex.get(snp);
				if (index != null && r.getPval() < 7.5E-7) {
					mat[index][i] = r.getPosterior();
				}
			}
		}
		
		TextFile outf = new TextFile(out, TextFile.W);
		String header = "-\tRegion\tGenes";
		for (int i = 0; i < filedesc.length; i++) {
			header += "\t" + filedesc[i];
		}
		
		outf.writeln(header);
		for (int i = 0; i < mat.length; i++) {
			String region = regions.get(i);
			String annot = regionToGenes.get(region);
			String ln = snps.get(i) + "\t" + region + "\t" + annot + "\t" + Strings.concat(mat[i], Strings.tab);
			outf.writeln(ln);
		}
		outf.close();
		
	}
	
}
