package nl.harmjanwestra.finemapping.annotation;

import htsjdk.tribble.readers.TabixReader;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.vcf.VCFTabix;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 2/21/17.
 */
public class IndelRSIdRewriter {

	public static void main(String[] args) {

		String tabixprefix = "/Data/Ref/beagle_1kg/1kg.phase3.v5a.chrCHR.vcf.gz";
		String samplefilter = "/Data/Ref/1kg-europeanpopulations.txt.gz";

//		String[] tables = new String[]{
//				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/tables/UpdateRSIds/SupplementaryTable7-2016-06-19-RA-LocusComparisonWithOkada-SignificantLoci.txt",
//				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/tables/UpdateRSIds/SupplementaryTable8-2016-06-19-T1D-LocusComparisonWithOnengutCC-SignificantLoci.txt",
//				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/tables/UpdateRSIds/SupplementaryTable10-2016-06-19-RA-LocusComparisonWithMeta-SignificantLoci.txt",
//				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/tables/UpdateRSIds/SupplementaryTable11-2016-06-19-T1D-LocusComparisonWithMeta-SignificantLoci.txt",
//				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/tables/UpdateRSIds/SupplementaryTable12-Crediblesets.txt",
//				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/tables/UpdateRSIds/SupplementaryTable13-Conditional.txt",
//				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/tables/UpdateRSIds/SupplementaryTable14-Crediblesets-FunctionalOverlap.txt",
//				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/tables/UpdateRSIds/SupplementaryTable16-ATACoverlap.txt",
//				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/tables/UpdateRSIds/SupplementaryTable17-Crediblesets-eQTLs.txt"
//		};
//
//
//		int[] headers = new int[]{
//				1,
//				1,
//				1,
//				1,
//				2,
//				2,
//				2,
//				2,
//				2
//		};
//
//		int[][] columns = new int[][]{
//				new int[]{1, 10},
//				new int[]{1, 10},
//				new int[]{1, 10},
//				new int[]{1, 10},
//				new int[]{5, 18, 31},
//				new int[]{3, 10, 17},
//				new int[]{3, 14, 25},
//				new int[]{3, 14, 25},
//				new int[]{3, 16, 29}
//		};

//		IndelRSIdRewriter v = new IndelRSIdRewriter();
//		for (int i = 0; i < tables.length; i++) {
//			String table = tables[i];
//			int[] column = columns[i];
//
//			try {
//				v.run(table, table + "-rewrite.txt", tabixprefix, samplefilter, headers[i], column);
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//		}

		String[] tables = new String[]{
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/wrongrsids/RA-assoc0.3-COSMO-merged-posterior-significantloci-75e7.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/wrongrsids/RA-assoc0.3-COSMO-merged-posterior.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/wrongrsids/RA-assoc0.3-COSMO-merged.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/wrongrsids/T1D-assoc0.3-COSMO-merged-posterior-significantloci-75e7.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/wrongrsids/T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/wrongrsids/T1D-assoc0.3-COSMO-merged.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/wrongrsids/META-assoc0.3-COSMO-merged-posterior-significantloci-75e7.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/wrongrsids/META-assoc0.3-COSMO-merged-posterior.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/wrongrsids/META-assoc0.3-COSMO-merged.txt.gz",
		};

		String[] tablesOut = new String[]{
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/RA-assoc0.3-COSMO-merged-posterior-significantloci-75e7.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/RA-assoc0.3-COSMO-merged-posterior.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/RA-assoc0.3-COSMO-merged.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/T1D-assoc0.3-COSMO-merged-posterior-significantloci-75e7.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/T1D-assoc0.3-COSMO-merged.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/META-assoc0.3-COSMO-merged-posterior-significantloci-75e7.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/META-assoc0.3-COSMO-merged-posterior.txt.gz",
				"/Sync/OneDrive/Postdoc/2016-03-RAT1D-Finemapping/Data/2017-03-25-SummaryStats/normal/META-assoc0.3-COSMO-merged.txt.gz",
		};

		IndelRSIdRewriter v = new IndelRSIdRewriter();
		for (int i = 0; i < tables.length; i++) {

			try {
				v.runOnAssocFiles(tables[i], tablesOut[i], tabixprefix, samplefilter);
			} catch (IOException e) {
				e.printStackTrace();
			}

		}

	}

	public void runOnAssocFiles(String infile, String outfile, String tabixprefix, String samplefilter) throws IOException {
		int idCol = 2;
		int combinedIdCol = 3;

		TextFile input = new TextFile(infile, TextFile.R);
		TextFile output = new TextFile(outfile, TextFile.W);

		// header
		output.writeln(input.readLine());
		String[] elems = input.readLineElems(TextFile.tab);
		int lnnr = 0;

		while (elems != null) {

			Chromosome chrobj = Chromosome.parseChr(elems[0]);
			Integer pos = Integer.parseInt(elems[1]);


			// get variants in reference
			Feature snp = new Feature(chrobj, pos, pos);
			ArrayList<VCFVariant> select = getSNP(snp, tabixprefix, samplefilter, pos);


			if (select.isEmpty()) {
				// ??
				System.out.println("Could not find variant at position " + pos + ". Keeping id: " + elems[idCol]);
			} else if (select.size() == 1) {
				//System.out.println("ln:\t" + lnnr + "\tcol:\t" + idCol + "\tSingle variant: " + pos + "\t\t-->\t\t" + select.get(0).getId());
				elems[idCol] = select.get(0).getId();
				elems[combinedIdCol] = chrobj.toString() + "_" + pos + "_" + select.get(0).getId();

			} else {
				// iterate variants
//						System.out.println("Multiple variants: " + posStr);
//						for (VCFVariant v : select) {
//							System.out.println(v.getChr() + "\t" + v.getPos() + "\t" + v.getId());
//						}
//						System.exit(-1);
			}

			if (lnnr % 1000 == 0) {
				System.out.println(lnnr + " lines parsed..");
			}

			output.writeln(Strings.concat(elems, Strings.tab));
			elems = input.readLineElems(TextFile.tab);
			lnnr++;
		}


		output.close();
		input.close();

	}

	public void run(String infile, String outfile, String tabixprefix, String samplefilter, int nrheaderlines, int[] columns) throws IOException {


		System.out.println("Processing " + infile);
		TextFile tf = new TextFile(infile, TextFile.R);
		TextFile tfout = new TextFile(outfile, TextFile.W);
		int ctr = 0;
		int lnnr = 0;
		while (ctr < nrheaderlines) {
			tfout.writeln(tf.readLine());
			ctr++;
			lnnr++;
		}

		String[] elems = tf.readLineElems(TextFile.tab);

		while (elems != null) {
			for (int i : columns) {
				if (i > elems.length) {
					System.err.println("Warning: column " + i + " is out of bounds at line " + lnnr + " for " + infile);
				} else {
					String posStr = elems[i];
					if (posStr.length() > 0) {
						String[] strElems = posStr.split("_");
						Chromosome chrobj = Chromosome.parseChr(strElems[0]);
						Integer pos = Integer.parseInt(strElems[1]);


						// get variants in reference
						Feature snp = new Feature(chrobj, pos, pos);
						ArrayList<VCFVariant> select = getSNP(snp, tabixprefix, samplefilter, pos);

						if (select.isEmpty()) {
							// ??
							System.out.println("Could not find variant at position " + posStr);
						} else if (select.size() == 1) {
							System.out.println("ln:\t" + lnnr + "\tcol:\t" + i + "\tSingle variant: " + posStr + "\t\t-->\t\t" + select.get(0).getId());
							elems[i] = chrobj.toString() + "_" + pos + "_" + select.get(0).getId();
						} else {
							// iterate variants
//						System.out.println("Multiple variants: " + posStr);
//						for (VCFVariant v : select) {
//							System.out.println(v.getChr() + "\t" + v.getPos() + "\t" + v.getId());
//						}
//						System.exit(-1);
						}
					}
				}
			}
			String ln = Strings.concat(elems, Strings.tab);
			tfout.writeln(ln);
			elems = tf.readLineElems(TextFile.tab);
			lnnr++;
		}
		tf.close();
		tfout.close();
		System.out.println("---");
		System.out.println();
	}

	private ArrayList<VCFVariant> getSNP(Feature snp, String tabixprefix, String samplefilter, Integer pos) throws IOException {

		String tabixfile = tabixprefix.replace("CHR", "" + snp.getChromosome().getNumber());
//		System.out.println(tabixfile);
		VCFTabix reader = new VCFTabix(tabixfile);
		boolean[] snpSampleFilter = reader.getSampleFilter(samplefilter);


		Feature snpF = new Feature(snp);
		snp.setStart(snpF.getStart() - 1);
		snp.setStop(snpF.getStop() + 1);

		TabixReader.Iterator inputSNPiter = reader.query(snpF);
		String snpStr = inputSNPiter.next();
		ArrayList<VCFVariant> select = new ArrayList<>();

		while (snpStr != null) {
			VCFVariant variant = new VCFVariant(snpStr, VCFVariant.PARSE.ALL, snpSampleFilter);
			if (variant.asFeature().overlaps(snpF)) {
				select.add(variant);
			}

			snpStr = inputSNPiter.next();
		}

		reader.close();

		ArrayList<VCFVariant> select2 = new ArrayList<>();
		for (VCFVariant v : select) {
			if (pos.equals(v.getPos())) {
				select2.add(v);
			}
		}
		return select2;
	}

}
