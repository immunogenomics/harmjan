package nl.harmjanwestra.finemapping.genotypesandimputation;

import nl.harmjanwestra.utilities.enums.DiseaseStatus;
import nl.harmjanwestra.utilities.individuals.Individual;
import nl.harmjanwestra.utilities.plink.PlinkFamFile;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Created by hwestra on 6/28/16.
 */
public class SupTable1SampleCounter {


	public static void main(String[] args) {

		String[] famfilesRA = new String[]{
				"/Data/ImmunoChip/2015-11-28-hg19/IndividualDatasets/RA-ES.fam",
				"/Data/ImmunoChip/2015-11-28-hg19/IndividualDatasets/RA-NL.fam",
				"/Data/ImmunoChip/2015-11-28-hg19/IndividualDatasets/RA-SEE.fam",
				"/Data/ImmunoChip/2015-11-28-hg19/IndividualDatasets/RA-SEU.fam",
				"/Data/ImmunoChip/2015-11-28-hg19/IndividualDatasets/RA-UK.fam",
				"/Data/ImmunoChip/2015-11-28-hg19/IndividualDatasets/RA-US.fam"
		};

		String ralist = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/SampleLists/RA-samplelist.txt";

		String[] famfilesT1D = new String[]{
				"/Data/ImmunoChip/2015-11-28-hg19/IndividualDatasets/T1D-UK.fam",
				"/Data/ImmunoChip/2015-11-28-hg19/IndividualDatasets/T1D-EUR.fam"
		};

		String t1dlist = "/Sync/Dropbox/2016-03-RAT1D-Finemappng/Data/SampleLists/T1D-samplelist.txt";

		SupTable1SampleCounter s = new SupTable1SampleCounter();
		try {
			for (int i = 0; i < famfilesRA.length; i++) {
				s.run(famfilesRA[i], ralist);
				System.out.println();
			}

			for (int i = 0; i < famfilesT1D.length; i++) {
				s.run(famfilesT1D[i], t1dlist);
				System.out.println();
			}

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String famfile, String listfile) throws IOException {


		PlinkFamFile pf = new PlinkFamFile(famfile);

		TextFile tf = new TextFile(listfile, TextFile.R);
		String ln = tf.readLine();
		HashSet<String> testedSampleList = new HashSet<String>();
		while (ln != null) {
			testedSampleList.add(ln);
			ln = tf.readLine();
		}
		tf.close();
		System.out.println(testedSampleList.size() + " total tested samples");

		ArrayList<Individual> allFamSamples = pf.getSamples();
		int nrTestedSamples = 0;
		int nrTestedSamplesCases = 0;
		int nrTestedSamplesControls = 0;
		int nrPseudocontrolsTested = 0;
		int nrExcludedCases = 0;
		int nrExcludedControls = 0;
		for (int i = 0; i < allFamSamples.size(); i++) {
			if (testedSampleList.contains(allFamSamples.get(i).getName())) {
				nrTestedSamples++;
				if (allFamSamples.get(i).getDiseaseStatus().equals(DiseaseStatus.CASE)) {
					nrTestedSamplesCases++;
				} else {
					nrTestedSamplesControls++;
				}
			} else {

				if (allFamSamples.get(i).getDiseaseStatus().equals(DiseaseStatus.CASE)) {
					nrExcludedCases++;
				} else {
					nrExcludedControls++;
				}


//				System.out.println(testedSamples.get(i).getName());
//				System.exit(-1);
			}
			if (testedSampleList.contains(allFamSamples.get(i).getName() + "-pseudo")) {
				nrPseudocontrolsTested++;
			}
		}

		System.out.println(nrTestedSamples + " samples in tested files: " + nrTestedSamplesCases + " \t " + nrTestedSamplesControls);
		System.out.println(nrPseudocontrolsTested + " pseudo controls");
		System.out.println(nrExcludedCases + " excluded cases");
		System.out.println(nrExcludedControls + " excluded controls");

		if (pf.getSamples().size() > pf.getFamilies().size()) {
			int nrTrios = 0;
			int nrunrelatedCases = 0;
			int nrunrelatedControls = 0;
			int nrCasesWhereParentIsCase = 0;
			int nrDuos = 0;
			for (Individual i : allFamSamples) {
				if (i.getFamily() == null || i.getFamily().getName().equals(i.getName())) {
					if (i.getDiseaseStatus().equals(DiseaseStatus.CONTROL)) {
						nrunrelatedControls++;
					} else {
						nrunrelatedCases++;
					}
				}
				boolean parentiscase = false;
				if (i.getFather() != null) {
					if (i.getFather().getDiseaseStatus().equals(DiseaseStatus.CASE)) {
						parentiscase = true;
					}
				}

				if (i.getMother() != null) {
					if (i.getMother().getDiseaseStatus().equals(DiseaseStatus.CASE)) {
						parentiscase = true;
					}
				}

				if (parentiscase) {
					nrCasesWhereParentIsCase++;
				}

				if (i.getMother() != null && i.getFather() != null) {
					nrTrios++;
				}

				if ((i.getMother() != null && i.getFather() == null) || (i.getMother() == null && i.getFather() != null)) {
					nrDuos++;
				}
			}
			System.out.println(nrunrelatedCases + " unrelated cases");
			System.out.println(nrunrelatedControls + " unrelated controls");
			System.out.println(nrCasesWhereParentIsCase + " cases with parent case");
			System.out.println(nrTrios + " full trios ");
			System.out.println(nrDuos + " duos");
		}
	}
}
