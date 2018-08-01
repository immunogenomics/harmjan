package nl.harmjanwestra.finemapping.genotypesandimputation;

import nl.harmjanwestra.utilities.enums.DiseaseStatus;
import nl.harmjanwestra.utilities.individuals.Family;
import nl.harmjanwestra.utilities.individuals.Individual;
import nl.harmjanwestra.utilities.plink.PlinkFamFile;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 6/29/16.
 */
public class PseudocontrolCheck {

	public static void main(String[] args) {
		PseudocontrolCheck c = new PseudocontrolCheck();
		try {
			c.filter();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

//	public void picksamples() throws IOException {
//
//
//		// if mom/dad == case --> pseudocontrol is included, but case is excluded
//		// if mom and dad == absent --> case is included, without pseudo control
//		// how many children controls are included that are not pseudocontrols
//
//		// load original fam path
//		// load new fam path
//
//
//		String newfamFile = "/Data/tmp/2016-06-24/T1D-recode-maf0005-ICRegions-samplenamefix-pseudo.vcf.gz.fam";
//		String origfamFile = "/Data/ImmunoChip/2015-11-28-hg19/IndividualDatasets/T1D-EUR.fam";
//		PlinkFamFile origfam = new PlinkFamFile(origfamFile);
//		String finalSampleListFile = "/Data/tmp/2016-06-29/testoutgwas-hj.txtsamplelist.txt";
//
//
//		HashSet<String> sampleTestHash = new HashSet<String>();
//		TextFile tf = new TextFile(finalSampleListFile, TextFile.R);
//		ArrayList<String> samplesDuringTesting = tf.readAsArrayList();
//		tf.close();
//
//		sampleTestHash.addAll(samplesDuringTesting);
//
//
//		PlinkFamFile newfam = new PlinkFamFile(newfamFile);
//
//		HashMap<String, Individual> newFamHash = hashInds(newfam);
//		HashMap<String, Individual> origFamHash = hashInds(origfam);
//
//
//		System.out.println();
//		// count the number of individuals in the original fam path that have a sibling in the tested data
//		ArrayList<Family> origFams = origfam.getFamilies();
//		System.out.println(origFams.size() + " families in dataset");
//		for (Family f : origFams) {
//			ArrayList<Individual> individuals = f.getIndividuals();
//		}
//
//
//		ArrayList<Individual> origInds = origfam.getSamples();
//		HashSet<String> visitedIndividuals = new HashSet<String>();
//
//		int nrFathersStillPresent = 0;
//		int nrMothersStillPresent = 0;
//		int nrFatherPseudosStillPresent = 0;
//		int nrMotherPseudosStillPresent = 0;
//
//		int nrSiblings = 0;
//		int nrSiblingsPseudo = 0;
//		int nrSiblingsCase = 0;
//
//
//		for (Individual i : origInds) {
//
//			// try to figure out if there any direct siblings being tested
//			HashSet<String> siblings = new HashSet<>();
//
//			if (i.getMother() != null) {
//				if (sampleTestHash.contains(i.getMother().getName())) {
//					nrMothersStillPresent++;
//				}
//				if (sampleTestHash.contains(i.getMother().getName() + "-pseudo")) {
//					nrMotherPseudosStillPresent++;
//				}
//				ArrayList<Individual> sibs = i.getMother().getChildren();
////				System.out.println(sibs.size());
//				for (Individual s : sibs) {
//					siblings.add(s.getName());
//				}
//
//			}
//			if (i.getFather() != null) {
//				if (sampleTestHash.contains(i.getFather().getName())) {
//					nrFathersStillPresent++;
//				}
//				if (sampleTestHash.contains(i.getFather().getName() + "-pseudo")) {
//					nrFatherPseudosStillPresent++;
//				}
//
//				ArrayList<Individual> sibs = i.getFather().getChildren();
//				for (Individual s : sibs) {
//					siblings.add(s.getName());
//				}
//			}
//
//			for (String sibling : siblings) {
//
//				if (!visitedIndividuals.contains(sibling)) {
//					if (sampleTestHash.contains(sibling)) {
//						nrSiblings++;
//					}
//					if (sampleTestHash.contains(sibling + "-pseudo")) {
//						nrSiblingsPseudo++;
//					}
//
//					if (origFamHash.get(sibling).getDiseaseStatus().equals(DiseaseStatus.CASE)) {
//						nrSiblingsCase++;
//					}
//
//					visitedIndividuals.add(sibling);
//				}
//
//
//			}
//
//		}
//
//		System.out.println(nrFathersStillPresent + " fathers tested ");
//		System.out.println(nrFatherPseudosStillPresent + " father pseudos tested ");
//
//		System.out.println(nrMothersStillPresent + " mothers tested ");
//		System.out.println(nrMotherPseudosStillPresent + " mother pseudos tested ");
//
//		System.out.println(nrSiblings + " siblings tested ");
//		System.out.println(nrSiblingsPseudo + " sibling pseudos tested ");
//
//		System.out.println(nrSiblingsCase + " siblings are case");
//
//
//		// get a sense of the actual number of individuals we should get
//		System.out.println();
//
//
//	}

	public HashMap<String, Individual> hashInds(PlinkFamFile f) {

		HashMap<String, Individual> hash = new HashMap<String, Individual>();
		ArrayList<Individual> samples = f.getSamples();
		for (Individual i : samples) {
			hash.put(i.getName(), i);
		}
		return hash;
	}

	public void filter() throws IOException {

		String origfamFile = "/Data/ImmunoChip/2015-11-28-hg19/IndividualDatasets/T1D-EUR.fam";
		PlinkFamFile origfam = new PlinkFamFile(origfamFile);
		ArrayList<Individual> allIndividuals = origfam.getSamples();

		System.out.println(allIndividuals.size() + "\tindividuals total");
		boolean unrelated = false;

		HashSet<Individual> individualsToExclude = new HashSet<Individual>();

		// exclude all parents
		for (Individual ind : allIndividuals) {
			if (ind.getMother() != null) {
				individualsToExclude.add(ind.getMother());
			}
			if (ind.getFather() != null) {
				individualsToExclude.add(ind.getFather());
			}
			if (ind.getDiseaseStatus().equals(DiseaseStatus.UNKNOWN)) {
				individualsToExclude.add(ind);
			}
		}

		System.out.println(individualsToExclude.size() + "\tindividuals removed after accounting for parents");
		System.out.println();

		ArrayList<Individual> cases = new ArrayList<>();
		ArrayList<Individual> controls = new ArrayList<>();
		for (Individual ind : allIndividuals) {
			if (!individualsToExclude.contains(ind)) {
				if (ind.getDiseaseStatus().equals(DiseaseStatus.CASE)) {
					cases.add(ind);
				} else if (ind.getDiseaseStatus().equals(DiseaseStatus.CONTROL)) {
					controls.add(ind);
				}
			}
		}
		System.out.println(cases.size() + "\tcases remain after removing parents");
		System.out.println(controls.size() + "\tcontrols remain after removing parents");
		System.out.println((controls.size() + cases.size()) + "\ttotal remain");
		System.out.println();

		for (Individual ind : cases) {
			if (!individualsToExclude.contains(ind)) {
				HashSet<Individual> sibs = new HashSet<>();
				if (ind.getMother() != null) {
					sibs.addAll(ind.getMother().getChildren());
				}
				if (ind.getFather() != null) {
					sibs.addAll(ind.getFather().getChildren());
				}
				Family fam = ind.getFamily();
				sibs.addAll(fam.getIndividuals());

				ArrayList<Individual> sibsArr = new ArrayList<>();
				sibsArr.addAll(sibs);
				for (Individual sib : sibsArr) {
					if (!sib.equals(ind)) {
						individualsToExclude.add(sib);
					}
				}
			}
		}

		cases = new ArrayList<>();
		controls = new ArrayList<>();
		for (Individual ind : allIndividuals) {
			if (!individualsToExclude.contains(ind)) {
				if (ind.getDiseaseStatus().equals(DiseaseStatus.CASE)) {
					cases.add(ind);
				} else if (ind.getDiseaseStatus().equals(DiseaseStatus.CONTROL)) {
					controls.add(ind);
				}
			}
		}
		System.out.println(cases.size() + "\tcases remain after removing cases' sibs");
		System.out.println(controls.size() + "\tcontrols remain after removing  cases' sibs");
		System.out.println((controls.size() + cases.size()) + "\ttotal remain");
		System.out.println();


		for (Individual ind : controls) {
			if (!individualsToExclude.contains(ind)) {
				HashSet<Individual> sibs = new HashSet<>();
				if (ind.getMother() != null) {
					sibs.addAll(ind.getMother().getChildren());
				}
				if (ind.getFather() != null) {
					sibs.addAll(ind.getFather().getChildren());
				}

				Family fam = ind.getFamily();
				sibs.addAll(fam.getIndividuals());

				ArrayList<Individual> sibsArr = new ArrayList<>();
				sibsArr.addAll(sibs);
				for (Individual sib : sibsArr) {
					if (!sib.equals(ind)) {
						individualsToExclude.add(sib);
					}
				}
			}
		}

		cases = new ArrayList<>();
		controls = new ArrayList<>();
		for (Individual ind : allIndividuals) {
			if (!individualsToExclude.contains(ind)) {
				if (ind.getDiseaseStatus().equals(DiseaseStatus.CASE)) {
					cases.add(ind);
				} else if (ind.getDiseaseStatus().equals(DiseaseStatus.CONTROL)) {
					controls.add(ind);
				}
			}
		}
		System.out.println(cases.size() + "\tcases remain after removing controls' sibs");
		System.out.println(controls.size() + "\tcontrols remain after removing controls' sibs");
		System.out.println((controls.size() + cases.size()) + "\ttotal remain");

		// add pseudos
		HashSet<Individual> pseudos = new HashSet<>();
		for (Individual ind : cases) {
			if (!individualsToExclude.contains(ind)) {
				if (ind.getMother() != null && ind.getFather() != null) {
					pseudos.add(ind);
				}
			}
		}

		ArrayList<String> finalSamples = new ArrayList<>();
		for (Individual ind : cases) {
			finalSamples.add(ind.getName());
		}
		for (Individual ind : controls) {
			finalSamples.add(ind.getName());
		}
		for (Individual ind : pseudos) {
			finalSamples.add(ind.getName() + "-pseudo");
		}


		System.out.println((pseudos.size() + controls.size()) + "\tcontrols and pseudo controls");
		int mothermissing = 0;
		int fathermissing = 0;
		int fatherormothermissing = 0;
		int bothmissing = 0;
		int bothpresent = 0;
		for (Individual ind : cases) {
			if (!pseudos.contains(ind)) {
				if (ind.getMother() == null) {
					mothermissing++;
				}
				if (ind.getFather() == null) {
					fathermissing++;
				}


				if (
						(ind.getFather() != null && ind.getMother() == null)
								||
								(ind.getFather() == null && ind.getMother() != null)
						) {
					fatherormothermissing++;
				}
				if (ind.getFather() == null && ind.getMother() == null) {
					bothmissing++;
				}
			}
			if (ind.getFather() != null && ind.getMother() != null) {
				bothpresent++;
			}
		}
		System.out.println();


		System.out.println(mothermissing + "\tcases with mother missing");
		System.out.println(fathermissing + "\tcases with father missing");
		System.out.println(bothmissing + "\tcases with both missing");
		System.out.println(bothpresent + "\tcases with both present");
		System.out.println(fatherormothermissing + "\tfather or mother missing");

		System.out.println(finalSamples.size() + " final samples");


		String generatedFamFile = "/Data/tmp/2016-06-29/T1D-recode-maf0005-ICRegions-samplenamefix-pseudo.vcf.gz.fam";
		PlinkFamFile generatedfam = new PlinkFamFile(generatedFamFile);

		HashSet<String> sampleset = new HashSet<>();
		sampleset.addAll(finalSamples);

		String generatedFamFileOut = "/Data/tmp/2016-06-29/T1D-recode-maf0005-ICRegions-samplenamefix-pseudo.vcf.gz-filtered.fam";
		TextFile famIn = new TextFile(generatedFamFile, TextFile.R);
		TextFile famOut = new TextFile(generatedFamFileOut, TextFile.W);

		int written = 0;
		HashSet<String> writtenHash = new HashSet<String>();
		String[] elems = elems = famIn.readLineElems(Strings.whitespace);
		while (elems != null) {
			if (elems.length > 2) {
				if (sampleset.contains(elems[1])) {
					famOut.writeln(Strings.concat(elems, Strings.tab));
					written++;
					writtenHash.add(elems[1]);
				}
			}
			elems = famIn.readLineElems(Strings.whitespace);
		}
		famIn.close();
		famOut.close();

		for (String sample : finalSamples) {
			if (!writtenHash.contains(sample)) {
				System.out.println(sample);
				System.exit(-1);
			}
		}
		System.out.println(written + " written");

	}

}
