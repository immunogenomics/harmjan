package nl.harmjanwestra.vcfutils;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.utilities.vcf.VCFGenotypeData;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Triple;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.*;
import java.util.regex.Pattern;

/**
 * Created by hwestra on 7/24/15.
 */
public class PseudoControls {

	public static void main(String[] args) {
		PseudoControls c = new PseudoControls();
		try {
//			for (int i = 0; i < 24; i++) {
//
//				String vcfin = "/Data/ImmunoChip/T1D/merged/merged-Chr" + i + ".vcf.gz";
//				String vcfout = "/Data/ImmunoChip/T1D/merged/merged-Chr" + i + "-pseudoCC.vcf.gz";
//				String famfileout = "/Data/ImmunoChip/T1D/binary/T1D-Chr" + i + ".fam";
//
//				if (i == 24) {
//					vcfin = "/Data/ImmunoChip/T1D/merged/merged-ChrX.vcf.gz";
//					vcfout = "/Data/ImmunoChip/T1D/merged/merged-ChrX-pseudoCC.vcf.gz";
//					famfileout = "/Data/ImmunoChip/T1D/binary/T1D-ChrX.fam";
//				}
//
//				String famfile = "/Data/ImmunoChip/T1D/binary/T1D.fam";
//
//
//				if (Gpio.exists(vcfin)) {
//					c.make(vcfin, vcfout, famfile, famfileout);
//					///	c.check(vcfout, famfile);
//				}
//				if (i == 24) {
//					// make males homozygous
//					//	c.makeMalesHomozygous(vcfin, famfile);
//				}
//				System.out.println();
//			}

			String covariates = "/Data/ImmunoChip/T1D/binary/covarmerged.txtmergedCovariates.txt";
			String famfile = "/Data/tmp/2016-03-11/T1D-recode-regionsfiltered-allelesfiltered-samplenamefix-pseudo.vcf.gz.fam";
			String output = "/Data/tmp/2016-03-11/2016-03-11-T1D-covarmerged.txtmergedCovariates-withPseudos.txt";
			c.makeNewCovariateFile(covariates, famfile, output);
			c.makeNewDiseaseStatusFile(famfile, famfile + "-diseaseStatus.txt");

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void makeNewDiseaseStatusFile(String famFile, String diseaseStatusFile) throws IOException {

		TextFile tf = new TextFile(famFile, TextFile.R);
		TextFile out = new TextFile(diseaseStatusFile, TextFile.W);
		String[] elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {

			if (elems.length > 1) {
				out.writeln(elems[1] + "\t" + elems[elems.length - 1]);

			}
			elems = tf.readLineElems(Strings.whitespace);
		}
		out.close();

		tf.close();
	}

	private void makeNewCovariateFile(String originalCovariates, String famFileWithPseudoControls, String output) throws IOException {

		// read covariate matrix

		// samples on rows
		// covariates on columns
		nl.harmjanwestra.utilities.legacy.genetica.math.matrix.DoubleMatrixDataset<String, String> ds = new nl.harmjanwestra.utilities.legacy.genetica.math.matrix.DoubleMatrixDataset<String, String>(originalCovariates);

		// get a list of pseudocontrols
		ArrayList<String> famFileSamples = new ArrayList<String>();
		ArrayList<String> pseudos = new ArrayList<String>();
		TextFile tf = new TextFile(famFileWithPseudoControls, TextFile.R);
		String[] elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {
			if (elems.length > 1) {
				String sample = elems[1];
				if (sample.endsWith("pseudo")) {
					pseudos.add(sample);
				}
				famFileSamples.add(sample);
			}
			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();

		System.out.println(pseudos.size() + " pseudos found in text path: " + famFileWithPseudoControls);

		// set the covariates of the pseudo controls to match those of the original kid
		ArrayList<String> samplesToAdd = new ArrayList<String>();
		ArrayList<double[]> dataToAdd = new ArrayList<double[]>();
		for (int i = 0; i < pseudos.size(); i++) {
			String pseudo = pseudos.get(i);
			String origsample = pseudo.replaceAll("-pseudo", "");
			Integer index = ds.hashRows.get(origsample);
			if (index != null) {
				double[] data = ds.rawData[index];
				dataToAdd.add(data);
				samplesToAdd.add(pseudo);
			}
		}

		System.out.println(samplesToAdd.size() + " pseudo samples will be added to covariate path");
		double[][] dataout = new double[ds.nrRows + samplesToAdd.size()][ds.nrCols];
		for (int r = 0; r < ds.nrRows; r++) {
			for (int c = 0; c < ds.nrCols; c++) {
				dataout[r][c] = ds.rawData[r][c];
			}
		}

		List<String> rows = ds.rowObjects;
		for (int r = ds.nrRows; r < dataout.length; r++) {
			int actual = r - ds.nrRows;
			dataout[r] = dataToAdd.get(actual);
			rows.add(samplesToAdd.get(actual));
		}

		DoubleMatrixDataset<String, String> dsout = new DoubleMatrixDataset<String, String>();
		dsout.setMatrix(dataout);
		try {
			dsout.setRowObjects(rows);
			dsout.setColObjects(ds.colObjects);
			dsout.save(output);
		} catch (Exception e) {
			e.printStackTrace();
		}


	}

	private void makeMalesHomozygous(String vcfIn, String famfile) throws IOException {


		System.out.println("Checking output....");
		System.out.println("VCF in: " + vcfIn);
		System.out.println("FAM in: " + famfile);

		VCFGenotypeData data = new VCFGenotypeData(vcfIn);
		ArrayList<String> vcfSamples = data.getSamples();

		System.out.println(vcfSamples.size() + "vcf samples");
		// load fam path
		// get the trios (only trios where kid is a case)
		ArrayList<Pair<String, Triple<String, String, String>>> famData = getTrios(famfile); // format: familyname<kid, mom, dad>
		HashMap<String, Boolean> genderPerSample = getGenderFromFamFile(famfile);

		HashMap<String, Integer> sampleMap = new HashMap<String, Integer>();
		for (int i = 0; i < vcfSamples.size(); i++) {
			sampleMap.put(vcfSamples.get(i), i);
		}

		Integer[][] momAndDad = new Integer[2][vcfSamples.size()];
		int nrTrios = 0;
		// boolean[] excludeSample = new boolean[vcfSamples.size()];
		HashMap<String, Integer> kidToPseudoControl = new HashMap<String, Integer>();
		ArrayList<String> pseudoCCNames = new ArrayList<String>();
		for (Pair<String, Triple<String, String, String>> family : famData) {
			String famId = family.getLeft();
			Triple<String, String, String> trio = family.getRight();
			String kid = trio.getLeft();
			String mom = trio.getMiddle();
			String dad = trio.getRight();

			Integer kidId = sampleMap.get(kid);
			if (kidId != null) {

				Integer momId = sampleMap.get(mom);
				Integer dadId = sampleMap.get(dad);
				if (momId != null && dadId != null) {
					if (!kid.endsWith("-pseudo")) {
						kidToPseudoControl.put(kid, nrTrios + vcfSamples.size());
						pseudoCCNames.add(kid + "-pseudo");

						momAndDad[0][kidId] = momId;
						momAndDad[1][kidId] = dadId;
						nrTrios++;
					}
				}
			}
		}

		System.out.println("Nr trios: " + nrTrios);

		System.out.println("Nr trios: " + nrTrios);

		// link samples to pseudocontrols
		Integer[] sampleToPseudo = new Integer[vcfSamples.size()];
		for (int i = 0; i < vcfSamples.size(); i++) {
			String sample = vcfSamples.get(i);
			if (sample.endsWith("-pseudo")) {
				String origSample = sample.replaceAll("-pseudo", "");
				Integer origSampleId = sampleMap.get(origSample);
				sampleToPseudo[origSampleId] = i;
			}
		}


		data.close();

		// index pseudocontrols
		data = new VCFGenotypeData(vcfIn);
		Integer[] kidToPseudoIdArr = new Integer[vcfSamples.size()];
		for (int i = 0; i < vcfSamples.size(); i++) {
			String kidName = vcfSamples.get(i);
			Integer pseudoId = kidToPseudoControl.get(kidName);

			kidToPseudoIdArr[i] = pseudoId;
		}

		while (data.hasNext()) {
			VCFVariant current = data.next();
			DoubleMatrix2D alleles = current.getGenotypeAllelesAsMatrix2D();


			for (int kidId = 0; kidId < alleles.rows(); kidId++) {

				Boolean gender = genderPerSample.get(vcfSamples.get(kidId));

				if (gender != null && gender) {
					// person is male
					if (alleles.getQuick(kidId, 0) != alleles.get(kidId, 1)) {
						// if the person is homozygous
						String b = (byte) alleles.getQuick(kidId, 0) + "/" + (byte) alleles.get(kidId, 1);
						System.err.println("heterozygous male found for variant: " + current.getId()
								+ "\t" + kidId
								+ "\t" + vcfSamples.get(kidId)
								+ "\t" + b);
					}
				}

			}
		}

	}

	private HashMap<String, Boolean> getGenderFromFamFile(String famfile) throws IOException {
		System.out.println("Loading trios from FAM path: " + famfile);
		HashMap<String, Boolean> output = new HashMap<String, Boolean>();
		TextFile tf = new TextFile(famfile, TextFile.R);
		String[] elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {

			/*
			Family ID
     Individual ID
     Paternal ID
     Maternal ID
     Sex (1=male; 2=female; other=unknown)
     Phenotype
			 */

			String genderStr = elems[4];
			Boolean gender = null;
			if (genderStr.equals("1")) {
				gender = true;
			} else if (genderStr.equals("2")) {
				gender = false;
			}

			output.put(elems[1], gender);


			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();
		System.out.println(output.size() + " trios found in FAM path");
		return output;
	}


	private void check(String vcfIn, String famfile) throws IOException {
		System.out.println("Checking output....");
		System.out.println("VCF in: " + vcfIn);
		System.out.println("FAM in: " + famfile);

		VCFGenotypeData data = new VCFGenotypeData(vcfIn);
		ArrayList<String> vcfSamples = data.getSamples();

		System.out.println(vcfSamples.size() + "vcf samples");
		// load fam path
		// get the trios (only trios where kid is a case)
		ArrayList<Pair<String, Triple<String, String, String>>> famData = getTrios(famfile); // format: familyname<kid, mom, dad>

		HashMap<String, Integer> sampleMap = new HashMap<String, Integer>();
		for (int i = 0; i < vcfSamples.size(); i++) {
			sampleMap.put(vcfSamples.get(i), i);
		}

		Integer[][] momAndDad = new Integer[2][vcfSamples.size()];
		int nrTrios = 0;
		// boolean[] excludeSample = new boolean[vcfSamples.size()];

		HashMap<String, Integer> kidToPseudoControl = new HashMap<String, Integer>();
		ArrayList<String> pseudoCCNames = new ArrayList<String>();
		for (Pair<String, Triple<String, String, String>> family : famData) {
			String famId = family.getLeft();
			Triple<String, String, String> trio = family.getRight();
			String kid = trio.getLeft();
			String mom = trio.getMiddle();
			String dad = trio.getRight();

			Integer kidId = sampleMap.get(kid);
			if (kidId != null) {

				Integer momId = sampleMap.get(mom);
				Integer dadId = sampleMap.get(dad);
				if (momId != null && dadId != null) {
					if (!kid.endsWith("-pseudo")) {
						kidToPseudoControl.put(kid, nrTrios + vcfSamples.size());
						pseudoCCNames.add(kid + "-pseudo");

						momAndDad[0][kidId] = momId;
						momAndDad[1][kidId] = dadId;
						nrTrios++;
					}
				}
			}
		}

		System.out.println("Nr trios: " + nrTrios);

		// link samples to pseudocontrols
		Integer[] sampleToPseudo = new Integer[vcfSamples.size()];
		for (int i = 0; i < vcfSamples.size(); i++) {
			String sample = vcfSamples.get(i);
			if (sample.endsWith("-pseudo")) {
				String origSample = sample.replaceAll("-pseudo", "");
				Integer origSampleId = sampleMap.get(origSample);
				sampleToPseudo[origSampleId] = i;
			}
		}


		data.close();

		// index pseudocontrols
		data = new VCFGenotypeData(vcfIn);
		Integer[] kidToPseudoIdArr = new Integer[vcfSamples.size()];
		for (int i = 0; i < vcfSamples.size(); i++) {
			String kidName = vcfSamples.get(i);
			Integer pseudoId = kidToPseudoControl.get(kidName);

			kidToPseudoIdArr[i] = pseudoId;
		}

		int finalNrSamples = vcfSamples.size() + nrTrios;
		while (data.hasNext()) {

			VCFVariant current = data.next();
			DoubleMatrix2D alleles = current.getGenotypeAllelesAsMatrix2D();
			byte[][] finalAlleles = new byte[2][finalNrSamples];


			for (int kidId = 0; kidId < alleles.rows(); kidId++) {

				Integer pseudoId = kidToPseudoIdArr[kidId];
				Integer momId = momAndDad[0][kidId];
				Integer dadId = momAndDad[1][kidId];

				// copy original alleles
				finalAlleles[0][kidId] = (byte) alleles.getQuick(kidId, 0);
				finalAlleles[1][kidId] = (byte) alleles.getQuick(kidId, 1);


				if (pseudoId != null) {
					// make the pseudo control
					makePseudoControl(kidId, momId, dadId, pseudoId, alleles, finalAlleles, current);

					// check generated pseudocontrol against the one in the VCF path
					int pseudoSample = sampleToPseudo[kidId];

					if (pseudoSample != -1) {
						byte[] pseudoAllelesInFile = getAlleles(alleles, pseudoSample);
						if (!(pseudoAllelesInFile[0] + "-" + pseudoAllelesInFile[1]).equals(finalAlleles[0][pseudoId] + "-" + finalAlleles[1][pseudoId])) {
							byte[] allelesMom = getAlleles(alleles, momId);
							byte[] allelesDad = getAlleles(alleles, dadId);
							byte[] allelesKid = getAlleles(alleles, kidId);

							System.err.println("Error in pseudo alleles written to disk: " +
									current.getId() + "\tdad - " + allelesDad[0] + "/" + allelesDad[1]
									+ "\tmom - " + allelesMom[0] + "/" + allelesMom[1]
									+ "\tkid - " + kidId + " - " + allelesKid[0] + "/" + allelesKid[1]
									+ "\tPseudo - " + finalAlleles[0][pseudoId] + "/" + finalAlleles[1][pseudoId]
									+ "\tPseudoVCF - " + pseudoSample + " - " + pseudoAllelesInFile[0] + "/" + pseudoAllelesInFile[1]
									+ "\t" + vcfSamples.get(kidId) + "\t" + vcfSamples.get(momId) + "\t" + vcfSamples.get(dadId)
							);
						}
					}
				}

			}
		}
	}

	// makes pseudocontrols from VCF and FAM path.
	public void make(String vcfIn, String vcfOut, String famfile, String famfileout) throws IOException {


		System.out.println("VCF in: " + vcfIn);
		System.out.println("VCF out: " + vcfOut);
		System.out.println("FAM in: " + famfile);
		System.out.println("FAM out: " + famfileout);

		VCFGenotypeData data = new VCFGenotypeData(vcfIn);
		ArrayList<String> vcfSamples = data.getSamples();

		System.out.println(vcfSamples.size());
		// load fam path
		// get the trios (only trios where kid is a case)
		ArrayList<Pair<String, Triple<String, String, String>>> famData = getTrios(famfile); // format: familyname<kid, mom, dad>

		HashMap<String, Integer> sampleMap = new HashMap<String, Integer>();
		for (int i = 0; i < vcfSamples.size(); i++) {
			sampleMap.put(vcfSamples.get(i), i);
		}

		// boolean[] excludeSample = new boolean[vcfSamples.size()];
		Triple<Pair<HashMap<String, Integer>, Integer[][]>, ArrayList<String>, Integer> trioDataInVCF = buildTrios(vcfOut, famData, sampleMap, vcfSamples, false);

		int nrTrios = trioDataInVCF.getRight();
		System.out.println("Nr trios in genotypes: " + nrTrios);
		if (nrTrios == 0) {
			System.out.println("Trying to find trios by appending famid");
			// retry by appending family IDs (plink output yay)
			trioDataInVCF = buildTrios(vcfOut, famData, sampleMap, vcfSamples, true);
			nrTrios = trioDataInVCF.getRight();
			if (nrTrios == 0) {
				System.out.println("No trios found in dataset. Exiting..");
				System.exit(-1);
			}
		}

		Pair<HashMap<String, Integer>, Integer[][]> meh = trioDataInVCF.getLeft();
		HashMap<String, Integer> kidToPseudoControl = meh.getLeft();
		Integer[][] momAndDad = meh.getRight();
		ArrayList<String> pseudoCCNames = trioDataInVCF.getMiddle();

		writeFamFile(vcfSamples, sampleMap, pseudoCCNames, famfile, famfileout);
		System.out.println("Nr trios in genotypes: " + nrTrios);

		TextFile tfOut = new TextFile(vcfOut, TextFile.W);
		// write the header.. extend with pseudocontrol names
		TextFile tfin = new TextFile(vcfIn, TextFile.R);
		String ln = tfin.readLine();
		while (ln != null) {
			if (ln.startsWith("##")) {
				tfOut.writeln(ln);
			} else if (ln.startsWith("#")) {
				String[] elems = ln.split("\t");
				String out = Strings.concat(elems, Strings.tab, 0, 9);
				for (int i = 9; i < elems.length; i++) {
					out += "\t" + elems[i];
				}

				for (int i = 0; i < pseudoCCNames.size(); i++) {
					out += "\t" + pseudoCCNames.get(i);
				}
				tfOut.writeln(out);
			} else {
				break;
			}
			ln = tfin.readLine();
		}
		data.close();

		// index pseudocontrols
		data = new VCFGenotypeData(vcfIn);
		Integer[] kidToPseudoIdArr = new Integer[vcfSamples.size()];
		for (int i = 0; i < vcfSamples.size(); i++) {
			String kidName = vcfSamples.get(i);
			Integer pseudoId = kidToPseudoControl.get(kidName);
			kidToPseudoIdArr[i] = pseudoId;
		}

		int finalNrSamples = vcfSamples.size() + nrTrios;
		System.out.println("Final Nr of samples: " + finalNrSamples);

		while (data.hasNext()) {

			VCFVariant current = data.next();
			DoubleMatrix2D alleles = current.getGenotypeAllelesAsMatrix2D();
			byte[][] finalAlleles = new byte[2][finalNrSamples];

			for (int kidId = 0; kidId < alleles.rows(); kidId++) {
				Integer momId = momAndDad[0][kidId];
				Integer dadId = momAndDad[1][kidId];
				fixMendelianErrors(kidId, momId, dadId, alleles);
			}


			for (int i = 0; i < finalAlleles[0].length; i++) {
				finalAlleles[0][i] = -1;
				finalAlleles[1][i] = -1;
			}

			for (int kidId = 0; kidId < alleles.rows(); kidId++) {
				Integer pseudoId = kidToPseudoIdArr[kidId];
				Integer momId = momAndDad[0][kidId];
				Integer dadId = momAndDad[1][kidId];

				// copy original alleles
				finalAlleles[0][kidId] = (byte) alleles.getQuick(kidId, 0);
				finalAlleles[1][kidId] = (byte) alleles.getQuick(kidId, 1);

				// make the pseudo control
				if (pseudoId != null) {
					// first round.. find mendelian errors: some parents can be children as well, mendelian errors in these
					// trios will propagate to the pseudocontrols


					// second round // fix mendelian errors
					makePseudoControl(kidId, momId, dadId, pseudoId, alleles, finalAlleles, current);
				}

//				if (pseudoId != null) {
//					byte[] allelesMom = getAlleles(alleles, momId);
//					byte[] allelesDad = getAlleles(alleles, dadId);
//					byte[] allelesKid = getAlleles(alleles, kidId);
//					System.out.println(current.getId() + "\tPSEUDO!: dad - " + allelesDad[0] + "/" + allelesDad[1]
//							+ "\tmom - " + allelesMom[0] + "/" + allelesMom[1]
//							+ "\tkid - " + kidId + " - " + allelesKid[0] + "/" + allelesKid[1]
//							+ "\tPseudo - " + pseudoId + " - " + finalAlleles[0][pseudoId] + "/" + finalAlleles[1][pseudoId]);
//				}
			}

			// write the results to the disque..
			StringBuilder builder = new StringBuilder(10000);
			// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
			builder.append(current.getChr()).append("\t")
					.append(current.getPos()).append("\t")
					.append(current.getId()).append("\t")
					.append(current.getAlleles()[0]).append("\t")
					.append(Strings.concat(current.getAlleles(), Strings.comma, 1, current.getAlleles().length)).append("\t")
					.append(current.getQual()).append("\t")
					.append(current.getFilter()).append("\t")
					.append(current.getInfoString()).append("\t")
					.append("GT");

			for (int j = 0; j < finalAlleles[0].length; j++) {
				if (finalAlleles[0][j] == -1) {
					builder.append("\t").append("./.");
				} else {
					builder.append("\t").append(finalAlleles[0][j]).append("/").append(finalAlleles[1][j]);
				}


			}

			tfOut.writeln(builder.toString());
//			}

		}
		System.out.println("Done parsing. Output is here: " + vcfOut);
		tfOut.close();
	}

	private Triple<Pair<HashMap<String, Integer>, Integer[][]>, ArrayList<String>, Integer> buildTrios(String vcfOut,
																									   ArrayList<Pair<String, Triple<String, String, String>>> famData,
																									   HashMap<String, Integer> sampleMap,
																									   ArrayList<String> vcfSamples,
																									   boolean appendFamId) throws IOException {
		int nrTrios = 0;
		Integer[][] momAndDad = new Integer[2][vcfSamples.size()];
		HashMap<String, Integer> kidToPseudoControl = new HashMap<>();
		ArrayList<String> pseudoCCNames = new ArrayList<>();
		TextFile individualsToExclude = new TextFile(vcfOut + "-parents.txt", TextFile.W);

		for (Pair<String, Triple<String, String, String>> family : famData) {
			String famId = family.getLeft();
			Triple<String, String, String> trio = family.getRight();
			String kid = trio.getLeft();
			String mom = trio.getMiddle();
			String dad = trio.getRight();

			Integer kidId = sampleMap.get(kid);
			if (appendFamId) {
				kidId = sampleMap.get(famId + "_" + kid);
			}
			if (kidId != null) {

				Integer momId = sampleMap.get(mom);
				Integer dadId = sampleMap.get(dad);

				if (momId != null) {
					individualsToExclude.writeln(mom);
				}
				if (dadId != null) {
					individualsToExclude.writeln(dad);
				}
				if (momId != null && dadId != null) {
					if (kidToPseudoControl.containsKey(kid)) {
						System.err.println("ERROR:  kid id already used: " + kid);
					}
					String pseudoname = kid + "-pseudo";
					kidToPseudoControl.put(kid, nrTrios + vcfSamples.size());

					pseudoCCNames.add(pseudoname);

					momAndDad[0][kidId] = momId;
					momAndDad[1][kidId] = dadId;
					nrTrios++;
				}
			}
		}
		individualsToExclude.close();
		return new Triple<>(new Pair<>(kidToPseudoControl, momAndDad), pseudoCCNames, nrTrios);
	}


	private void writeFamFile(ArrayList<String> vcfSamples, HashMap<String, Integer> sampleMap, ArrayList<String> pseudoControlNames, String famfile, String famfileout) throws IOException {
		TextFile famIn = new TextFile(famfile, TextFile.R);


		HashMap<String, String> sampleToFam = new HashMap<String, String>();
		String line = famIn.readLine();
		while (line != null) {
			String[] elems = Strings.whitespace.split(line);
			String sample = elems[1];
			sampleToFam.put(sample, line);
			line = famIn.readLine();
		}


		TextFile famout = new TextFile(famfileout, TextFile.W);
		int nr = 0;
		for (int i = 0; i < vcfSamples.size(); i++) {
			String sample = vcfSamples.get(i);
			String fam = sampleToFam.get(sample);
			if (fam != null) {
				famout.writeln(fam);
				nr++;
			}
		}
		System.out.println(nr + " fam lines written, " + vcfSamples.size() + " samples in vcf");


		/*
		Family ID
     Individual ID
     Paternal ID
     Maternal ID
     Sex (1=male; 2=female; other=unknown)
     Phenotype
		 */
		for (int i = 0; i < pseudoControlNames.size(); i++) {
			String pseudoControl = pseudoControlNames.get(i);
			String origSample = pseudoControl.replaceAll("-pseudo", "");
			Integer id = sampleMap.get(origSample);
			if (id == null) {
				System.err.println("ERROR while writing FAM path");
			} else {
				String fam = sampleToFam.get(origSample);
				String[] elems = Strings.whitespace.split(fam);

				elems[1] = pseudoControl;

				int sex = Integer.parseInt(elems[4]);
				if (sex == 1) {
					elems[4] = "" + 2;
				} else if (sex == 2) {
					elems[4] = "" + 1;
				} else {
					elems[4] = "" + (-9);
				}

				elems[5] = "" + 1;

				String out = Strings.concat(elems, Pattern.compile(" "));

				famout.writeln(out);


			}

		}

		famIn.close();
		famout.close();

	}

	private void fixMendelianErrors(int kidId, Integer momId, Integer dadId, DoubleMatrix2D alleles) {
		if (momId == null || dadId == null) {
			// incomplete trio..
			// set pseudo to missing

		} else {

			byte[] allelesMom = getAlleles(alleles, momId);
			byte[] allelesDad = getAlleles(alleles, dadId);
			byte[] allelesKid = getAlleles(alleles, kidId);

			if (allelesMom[0] == -1 || allelesDad[0] == -1 || allelesKid[0] == -1) {
				// missing genotypes
				// set pseudo to missing..

			} else {
				HashSet<Pair<Byte, Byte>> allowedGenotypes = new HashSet<Pair<Byte, Byte>>(4);
				HashMap<Byte, Byte> uniqueAlleles = new HashMap<Byte, Byte>(3);

				// DETECT MENDELIAN ERRORS:
				// make a list of available alleles
				for (int a = 0; a < allelesDad.length; a++) {
					Byte nr = uniqueAlleles.get(allelesDad[a]);
					if (nr == null) {
						nr = 1;
					} else {
						nr++;
					}
					uniqueAlleles.put(allelesDad[a], nr);
					nr = uniqueAlleles.get(allelesMom[a]);
					if (nr == null) {
						nr = 1;
					} else {
						nr++;
					}
					uniqueAlleles.put(allelesMom[a], nr);
				}

				// cross with mom: generate all possible genotypes
				for (int a = 0; a < allelesDad.length; a++) {
					for (int b = 0; b < allelesMom.length; b++) {
						allowedGenotypes.add(new Pair<Byte, Byte>(allelesDad[a], allelesMom[b]));
						allowedGenotypes.add(new Pair<Byte, Byte>(allelesMom[b], allelesDad[a]));
					}
				}

				Pair<Byte, Byte> allelesPairKid = new Pair<Byte, Byte>(allelesKid[0], allelesKid[1]);
				if (!allowedGenotypes.contains(allelesPairKid)) {
					// mendelian error detected
//					System.err.println(current.getId() + "\tMendelian error: dad - " + allelesDad[0] + "/" + allelesDad[1]
//							+ "\tmom - " + allelesMom[0] + "/" + allelesMom[1]
//							+ "\tkid - " + allelesKid[0] + "/" + allelesKid[1]);

					// should we set them all to missing in this case?

					// set pseudo and kid to missing

					alleles.set(kidId, 0, -1);//[0][kidId] = -1;
					alleles.set(kidId, 1, -1);// alleles[1][kidId] = -1;
				}
			}
		}
	}

	private void makePseudoControl(int kidId, Integer momId, Integer dadId, int pseudoId, DoubleMatrix2D alleles, byte[][] finalAlleles, VCFVariant current) {
		if (momId == null || dadId == null) {
			// incomplete trio..
			// set pseudo to missing
			finalAlleles[0][pseudoId] = -1;
			finalAlleles[1][pseudoId] = -1;
		} else {

			byte[] allelesMom = getAlleles(alleles, momId);
			byte[] allelesDad = getAlleles(alleles, dadId);
			byte[] allelesKid = getAlleles(alleles, kidId);

			if (allelesMom[0] == -1 || allelesDad[0] == -1 || allelesKid[0] == -1) {
				// missing genotypes
				// set pseudo to missing..
				finalAlleles[0][pseudoId] = -1;
				finalAlleles[1][pseudoId] = -1;
			} else {
				HashSet<Pair<Byte, Byte>> allowedGenotypes = new HashSet<Pair<Byte, Byte>>(4);
				HashMap<Byte, Byte> uniqueAlleles = new HashMap<Byte, Byte>(3);

				// DETECT MENDELIAN ERRORS:
				// make a list of available alleles
				for (int a = 0; a < allelesDad.length; a++) {
					Byte nr = uniqueAlleles.get(allelesDad[a]);
					if (nr == null) {
						nr = 1;
					} else {
						nr++;
					}
					uniqueAlleles.put(allelesDad[a], nr);
					nr = uniqueAlleles.get(allelesMom[a]);
					if (nr == null) {
						nr = 1;
					} else {
						nr++;
					}
					uniqueAlleles.put(allelesMom[a], nr);
				}

				// cross with mom: generate all possible genotypes
				for (int a = 0; a < allelesDad.length; a++) {
					for (int b = 0; b < allelesMom.length; b++) {
						allowedGenotypes.add(new Pair<Byte, Byte>(allelesDad[a], allelesMom[b]));
						allowedGenotypes.add(new Pair<Byte, Byte>(allelesMom[b], allelesDad[a]));
					}
				}

				Pair<Byte, Byte> allelesPairKid = new Pair<Byte, Byte>(allelesKid[0], allelesKid[1]);
				if (!allowedGenotypes.contains(allelesPairKid)) {
					// mendelian error detected
//					System.err.println(current.getId() + "\tMendelian error: dad - " + allelesDad[0] + "/" + allelesDad[1]
//							+ "\tmom - " + allelesMom[0] + "/" + allelesMom[1]
//							+ "\tkid - " + allelesKid[0] + "/" + allelesKid[1]);

					// should we set them all to missing in this case?

					// set pseudo and kid to missing
					finalAlleles[0][pseudoId] = -1;
					finalAlleles[1][pseudoId] = -1;
					finalAlleles[0][kidId] = -1;
					finalAlleles[1][kidId] = -1;
				} else {
					Set<Byte> keys = uniqueAlleles.keySet();
					for (int a = 0; a < allelesKid.length; a++) {
						int nr = uniqueAlleles.get(allelesKid[a]); // get count for each allele
						nr--; // claim one
						uniqueAlleles.put(allelesKid[a], (byte) nr); // put it back
					}


					int allelesSet = 0;
					for (Byte k : keys) { // key to an allele
						int nr = uniqueAlleles.get(k); // the number of this type of allele that is left
						if (nr == 2) { // 2 alleles of this type left: pseudoCC must be homozygote for this allele
							finalAlleles[allelesSet][pseudoId] = k;
							allelesSet++;
							nr--; // claim one
							uniqueAlleles.put(k, (byte) nr); // put it back
							finalAlleles[allelesSet][pseudoId] = k;
							allelesSet++;
							nr--; // claim one
							uniqueAlleles.put(k, (byte) nr); // put it back
							break;
						} else if (nr == 1) {
							finalAlleles[allelesSet][pseudoId] = k;
							nr--; // claim one
							uniqueAlleles.put(k, (byte) nr); // put it back
							allelesSet++;

						}
					}


				}
			}


		}
	}

	private byte[] getAlleles(DoubleMatrix2D alleles, int id) {
		byte[] al = new byte[2];
		al[0] = (byte) alleles.getQuick(id, 0);
		al[1] = (byte) alleles.getQuick(id, 1);
		return al;
	}

	private ArrayList<Pair<String, Triple<String, String, String>>> getTrios(String famfile) throws IOException {
		System.out.println("Loading trios from FAM path: " + famfile);
		ArrayList<Pair<String, Triple<String, String, String>>> output = new ArrayList<Pair<String, Triple<String, String, String>>>();
		TextFile tf = new TextFile(famfile, TextFile.R);
		String[] elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {

			String family = elems[0];

			String kid = elems[1];
			String dad = elems[2];
			String mom = elems[3];

			if (!dad.equals("0") && !mom.equals("0")) {
				Integer control = Integer.parseInt(elems[5]);
//				System.out.println(kid + "\t" + control);
				if (control.equals(2)) {
					Triple<String, String, String> t = new Triple<String, String, String>(kid, mom, dad);
					Pair<String, Triple<String, String, String>> p = new Pair<String, Triple<String, String, String>>(family, t);
					output.add(p);
				}
			}
			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();
		System.out.println(output.size() + " trios found in FAM path");
		return output;
	}

}
