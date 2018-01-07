package nl.harmjanwestra.finemapping.rebuttal1;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class R1Q1CompareSignificantVariants {
	
	public static void main(String[] args) {
		
		String[] datasets = new String[]{
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\META-assoc0.3-COSMO-merged-posterior.txt.gz",
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
				"D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
		};
		
		
		String[] dsnames = new String[]{
				"META",
				"RA",
				"T1D"
		};
		
		
		String outfile = "D:\\Sync\\OneDrive\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\PValueComparisonToMeta.txt";
		double threshold = 7.5E-7;
		R1Q1CompareSignificantVariants c = new R1Q1CompareSignificantVariants();
		try {
//			c.compare(datasets, dsnames, threshold, outfile);
			
			String crediblesetsizes = "D:\\Sync\\Dropbox\\FineMap\\2018-01-Rebuttal\\Tables\\Rev1-Q1-SizesOfCredibleSets.txt";
			outfile = "D:\\Sync\\Dropbox\\FineMap\\2018-01-Rebuttal\\Tables\\Rev1-Q1-SizesOfCredibleSets-";
			c.conmpareCredibleSetSizes(crediblesetsizes, outfile);
			
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private void conmpareCredibleSetSizes(String crediblesetsizes, String out) throws IOException {
		
		TextFile tf = new TextFile(crediblesetsizes, TextFile.R);
		String ln = tf.readLine();
		ln = tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		
		int raregioncol = 2;
		int ranrvarscol = 3;
		
		int t1dregioncol = 4;
		int t1dnrvarscol = 5;
		
		int metaregioncol = 6;
		int metanrvarscol = 7;
		
		int nrSigRa = 0;
		int nrSigT1D = 0;
		int nrSmallerEqualRA = 0;
		int nrSmallerEqualT1D = 0;
		int nrSmallerEqualBoth = 0;
		int nrSmallerRA = 0;
		int nrSmallerT1D = 0;
		int nrSmallerBoth = 0;
		int nrSigBoth = 0;
		
		int nroverallSigRA = 0;
		int nroverallSmallerRA = 0;
		
		int nroverallSigT1D = 0;
		int nroverallSmallerT1D = 0;
		
		TextFile tfoutra = new TextFile(out + "ra.txt", TextFile.W);
		TextFile tfoutt1d = new TextFile(out + "t1d.txt", TextFile.W);
		TextFile tfoutshared = new TextFile(out + "all.txt", TextFile.W);
		while (elems != null) {
			
			
			String region = elems[0];
			String genes = elems[1];
			Boolean ra1 = Boolean.parseBoolean(elems[raregioncol]);
			Boolean t1d1 = Boolean.parseBoolean(elems[t1dregioncol]);
			Boolean meta1 = Boolean.parseBoolean(elems[metaregioncol]);
			
			Integer rasize = Integer.parseInt(elems[ranrvarscol]);
			Integer t1dsize = Integer.parseInt(elems[t1dnrvarscol]);
			Integer metasize = Integer.parseInt(elems[metanrvarscol]);
			
			
			if (ra1) {
				nroverallSigRA++;
				if (metasize <= rasize) {
					nroverallSmallerRA++;
				}
			}
			
			if (t1d1) {
				nroverallSigT1D++;
				if (metasize <= t1dsize) {
					nroverallSmallerT1D++;
				}
			}
			
			if (ra1 && meta1) {
				tfoutra.writeln(region + "\t" + genes + "\t" + rasize + "\t" + metasize + "\t" + (rasize - metasize));
				if (metasize <= rasize) {
					nrSmallerRA++;
				}
				if (metasize < rasize) {
					nrSmallerEqualRA++;
				}
				nrSigRa++;
			}
			
			if (t1d1 && meta1) {
				tfoutt1d.writeln(region + "\t" + genes + "\t" + t1dsize + "\t" + metasize + "\t" + (t1dsize - metasize));
				if (metasize <= t1dsize) {
					nrSmallerT1D++;
				}
				if (metasize < t1dsize) {
					nrSmallerEqualT1D++;
				}
				nrSigT1D++;
			}
			
			if (t1d1 && meta1 && ra1) {
				tfoutshared.writeln(region + "\t" + genes + "\t" + t1dsize + "\t" + rasize + "\t" + metasize + "\t" + (t1dsize - metasize) + "\t" + (rasize - metasize));
				if (metasize <= t1dsize && metasize <= rasize) {
					nrSmallerBoth++;
				}
				if (metasize < t1dsize && metasize < rasize) {
					nrSmallerEqualBoth++;
				}
				nrSigBoth++;
			}
			
			elems = tf.readLineElems(TextFile.tab);
		}
		tfoutra.writeln(nrSmallerRA + "\t" + nrSmallerEqualRA + "\t" + ((double) nrSmallerRA / nrSigRa) + "\t" + ((double) nrSmallerEqualRA / nrSigRa));
		tfoutra.writeln(nroverallSmallerRA + "\t" + nroverallSigRA + "\t" + ((double) nroverallSmallerRA / nroverallSigRA));
		tfoutt1d.writeln(nrSmallerT1D + "\t" + nrSmallerEqualT1D + "\t" + ((double) nrSmallerT1D / nrSigT1D) + "\t" + ((double) nrSmallerEqualT1D / nrSigT1D));
		tfoutt1d.writeln(nroverallSmallerT1D + "\t" + nroverallSigT1D + "\t" + ((double) nroverallSmallerT1D / nroverallSigT1D));
		
		tfoutshared.writeln(nrSmallerBoth + "\t" + nrSmallerEqualBoth + "\t" + ((double) nrSmallerBoth / nrSigBoth) + "\t" + ((double) nrSmallerEqualBoth / nrSigBoth));
		tfoutra.close();
		tfoutt1d.close();
		tfoutshared.close();
		tf.close();
	}
	
	private void compare(String[] datasets, String[] dsnames, double threshold, String outfile) throws IOException {
		
		double log10threshold = -Math.log10(threshold);
		
		AssociationFile f = new AssociationFile();
		ArrayList<ArrayList<AssociationResult>> t = new ArrayList<>();
		HashSet<String> vars = new HashSet<>();
		ArrayList<HashMap<String, AssociationResult>> resultmap = new ArrayList<>();
		for (int d = 0; d < datasets.length; d++) {
			ArrayList<AssociationResult> res = f.read(datasets[d]);
			t.add(res);
			HashMap<String, AssociationResult> map = new HashMap<>();
			for (AssociationResult r : res) {
				vars.add(r.getSnp().toString());
				map.put(r.getSnp().toString(), r);
			}
			resultmap.add(map);
		}
		
		TextFile out = new TextFile(outfile, TextFile.W);
		String header = "Variant\tSignificantInDs";
		for (int d = 0; d < dsnames.length; d++) {
			header += "\tP" + dsnames[d] + "\tPosterior-" + dsnames[d];
		}
		out.writeln(header);
		
		for (String var : vars) {
			String ln = var;
			
			String significantin = "";
			String dsln = "";
			for (int d = 0; d < dsnames.length; d++) {
				
				HashMap<String, AssociationResult> map = resultmap.get(d);
				AssociationResult r = map.get(var);
				if (r == null) {
					dsln += "\t-\t-";
				} else {
					if (r.getPval() < threshold) {
						if (significantin.length() > 0) {
							significantin += "," + dsnames[d];
						} else {
							significantin += "" + dsnames[d];
						}
					}
					dsln += "\t" + r.getPval() + "\t" + r.getPosterior();
				}
			}
			
			if (significantin.length() == 0) {
				significantin = "-";
			}
			ln += "\t" + significantin + dsln;
			
			out.writeln(ln);
		}
		out.close();
	}
}
