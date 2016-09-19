package nl.harmjanwestra.finemapping;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.graphics.ChromosomePlot;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;

/**
 * Created by hwestra on 8/26/15.
 */
public class ChrPlot {

	public static void main(String[] args) {
//		ChrPlot c = new ChrPlot();
//		try {
//			c.reannotateVariants();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		System.exit(-1);

		String cytoband = "/Data/Ref/Annotation/UCSC/cytoBand.txt";
		ChromosomePlot plot = null;
		try {
			plot = new ChromosomePlot("/Data/tmp/chrplot-tb.pdf", 1200, 1200);
			plot.setMargin(200);
//			String[] gfffiles = new String[]{
//					"/Data/ImmunoBase//Hs_GRCh38-AA-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-AS-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-ATD-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-CEL-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-CRO-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-IBD-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-JIA-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-MS-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-NAR-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-PBC-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-PSC-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-PSO-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-RA-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-SJO-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-SLE-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-SSC-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-T1D-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-UC-assoc_genesGFF",
//					"/Data/ImmunoBase//Hs_GRCh38-VIT-assoc_genesGFF"};
			String[] gfffiles = new String[]{
					"/Data/ImmunoBase/IB//Hs_GRCh38-AA-assoc_genesGFF",
					"/Data/ImmunoBase/IB//Hs_GRCh38-IBD-assoc_genesGFF",
					"/Data/ImmunoBase/IB//Hs_GRCh38-NAR-assoc_genesGFF",
					"/Data/ImmunoBase/IB//Hs_GRCh38-PSC-assoc_genesGFF",
					"/Data/ImmunoBase/IB//Hs_GRCh38-SJO-assoc_genesGFF",
					"/Data/ImmunoBase/IB//Hs_GRCh38-SSC-assoc_genesGFF",
					"/Data/ImmunoBase/IB//Hs_GRCh38-VIT-assoc_genesGFF"
			};

//			gfffiles = new String[]{
//					"/Data/tmp/tb-edit.txt"};
//			String sequencedRegionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/2015-03-30-allRegionsWith50PercentOfSamplesAbove20x.bed";
//			plot.plot(cytoband, gfffiles, false, null);



			plot = new ChromosomePlot("/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ChromosomePlots/ra-t1d.pdf", 1200, 1200);
			plot.setMargin(200);
			gfffiles = new String[]{
					"/Data/Ref/ImmunoBase/ImmunoBase/Hs_GRCh37-RA-assoc_genesGFF",
					"/Data/Ref/ImmunoBase/ImmunoBase/Hs_GRCh37-T1D-assoc_genesGFF"
			};
			plot.plot(cytoband, gfffiles, true, null);

			String regions = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/LocusDefinitions/OverlapWithSequencing/ICLociOverlappingWithSequencingRegions.bed";
			plot = new ChromosomePlot("/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/ChromosomePlots/ra-t1d-SequencedRegions.pdf", 1200, 1200);
			plot.setMargin(200);
			gfffiles = new String[]{
					"/Data/Ref/ImmunoBase/ImmunoBase/Hs_GRCh37-RA-assoc_genesGFF",
					"/Data/Ref/ImmunoBase/ImmunoBase/Hs_GRCh37-T1D-assoc_genesGFF"
			};



			plot.plot(cytoband, gfffiles, true, regions);

//			plot.heatmap(gfffiles, true);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}


	}

	public void reannotateVariants() throws IOException {
		String filein = "/Data/tmp/tb.txt";
		String fileout = "/Data/tmp/tb-edit.txt";
		String dbsnpvcf = "/Data/Projects/2014-FR-Reseq/GATKResources/dbsnp_138.b37.vcf.gz";


		HashMap<String, String> variants = new HashMap<String, String>();
		TextFile tf = new TextFile(filein, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length == 10) {
				variants.put(elems[9], null);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		System.out.println(variants.size() + " variants loaded");

		TextFile vcf = new TextFile(dbsnpvcf, TextFile.R);

		System.out.println("parsing DBSNP VCF: " + dbsnpvcf);
		String[] lineElems = vcf.readLineElems(TextFile.tab);
		int ln = 0;
		while (lineElems != null) {

			if (lineElems[0].startsWith("#")) {
				// header
			} else {

				String rs = lineElems[2];
				if (variants.containsKey(rs)) {
					variants.put(rs, lineElems[1]);
				}
			}

			ln++;

			if (ln % 2500000 == 0) {
				System.out.println(ln + " positions parsed...");
			}
			lineElems = vcf.readLineElems(TextFile.tab);
		}
		vcf.close();

		TextFile tfout = new TextFile(fileout, TextFile.W);
		tf = new TextFile(filein, TextFile.R);
		elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length == 10) {
				String variant = elems[9];
				String pos = variants.get(variant);

				if (pos != null) {
					elems[3] = "" + (Integer.parseInt(pos) - 1);
					elems[4] = pos;
				}
				tfout.writeln(Strings.concat(elems, Strings.tab));
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		tfout.close();


	}
}
