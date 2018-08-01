package nl.harmjanwestra.ngs;

import nl.harmjanwestra.utilities.plink.PedAndMapFunctions;

import java.io.IOException;

/**
 * Created by hwestra on 11/30/15.
 */
public class tmp {
	public static void main(String[] args) {

		try {
			PedAndMapFunctions p = new PedAndMapFunctions();
			p.filterFAM("/Data/Projects/2014-FR-Reseq/2015-finalRun/AllSequencedImmunoChipIDs.txt",
					"/Data/ImmunoChip/2015-11-28-hg19/Merged/T1D.fam",
					"/Data/ImmunoChip/2015-11-28-hg19/Merged/T1D.fam.sequencedIndividuals.txt");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
