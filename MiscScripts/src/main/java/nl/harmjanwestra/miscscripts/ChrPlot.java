package nl.harmjanwestra.miscscripts;

import com.lowagie.text.DocumentException;
import nl.harmjanwestra.utilities.graphics.ChromosomePlot;

import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * Created by hwestra on 8/26/15.
 */
public class ChrPlot {

	public static void main(String[] args) {
		String cytoband = "/Data/Annotation/UCSC/cytoBand.txt";
		ChromosomePlot plot = null;
		try {
			plot = new ChromosomePlot("/Data/tmp/chrplot.pdf", 3400, 2400);
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

			gfffiles = new String[]{
					"/Data/ImmunoBase/OD/Hs_GRCh38-AS-assoc_genesGFF",
					"/Data/ImmunoBase/OD/Hs_GRCh38-ATD-assoc_genesGFF",
					"/Data/ImmunoBase/OD/Hs_GRCh38-CEL-assoc_genesGFF",
					"/Data/ImmunoBase/OD/Hs_GRCh38-CRO-assoc_genesGFF",
					"/Data/ImmunoBase/OD/Hs_GRCh38-JIA-assoc_genesGFF",
					"/Data/ImmunoBase/OD/Hs_GRCh38-MS-assoc_genesGFF",
					"/Data/ImmunoBase/OD/Hs_GRCh38-PBC-assoc_genesGFF",
					"/Data/ImmunoBase/OD/Hs_GRCh38-PSO-assoc_genesGFF",
					"/Data/ImmunoBase/OD/Hs_GRCh38-RA-assoc_genesGFF",
					"/Data/ImmunoBase/OD/Hs_GRCh38-SLE-assoc_genesGFF",
					"/Data/ImmunoBase/OD/Hs_GRCh38-T1D-assoc_genesGFF",
					"/Data/ImmunoBase/OD/Hs_GRCh38-UC-assoc_genesGFF"
			};
			String sequencedRegionFile = "/Data/Projects/2014-FR-Reseq/2015-finalRun/2015-04-07-Analysis/2015-03-30-allRegionsWith50PercentOfSamplesAbove20x.bed";
			plot.plot(cytoband, gfffiles, true, null);
//			plot.heatmap(gfffiles, true);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
