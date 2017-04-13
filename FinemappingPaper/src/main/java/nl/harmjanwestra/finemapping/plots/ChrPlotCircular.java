package nl.harmjanwestra.finemapping.plots;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.GFFFile;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.CircularHeatmapPanel;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Triple;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by hwestra on 9/13/16.
 */
public class ChrPlotCircular {

	public static void main(String[] args) {
		String[] gfffiles = new String[]{
				"/Data/ImmunoChip/ImmunoBase/RA.csv",
				"/Data/ImmunoChip/ImmunoBase/T1D.csv"
		};

		String[] datasetname = new String[]{"RA", "T1D"};
		String output = "/Data/tmp/circchrplot.pdf";


		ChrPlotCircular chrPlotCircular = new ChrPlotCircular();
		try {
			chrPlotCircular.run(gfffiles, true, datasetname, output);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
	}

	public void run(String[] gffLocusFiles, boolean onlySuggestedLoci, String[] datasetnames, String output) throws IOException, DocumentException {
// plot the chromosomes


		Chromosome[] allChr = Chromosome.values();

		int ctr = 0;

		ArrayList<ArrayList<Feature>> associatedLociPerDisease = new ArrayList<ArrayList<Feature>>();
		GFFFile gff = new GFFFile();
		boolean sortstuff = true;
		for (int f = 0; f < gffLocusFiles.length; f++) {

			String gfffile = gffLocusFiles[f];
			if (gfffile.endsWith("GFF")) {
				associatedLociPerDisease.add(gff.readGFF(gfffile, onlySuggestedLoci));
			} else {
				sortstuff = false;
				ArrayList<Feature> regions = new ArrayList<>();
				TextFile tf = new TextFile(gfffile, TextFile.R);
				tf.readLine();
				String[] elems = tf.readLineElems(TextFile.comma);
				while (elems != null) {
					if (elems.length > 2) {

						String pos = elems[1].replaceAll("\"", "");
						// chr1:37731957-38002404
						String[] posElems = pos.split(":");

						String[] posElems2 = posElems[1].split("-");

						Chromosome chr = Chromosome.parseChr(posElems[0]);
						Integer start = Integer.parseInt(posElems2[0]);
						Integer stop = Integer.parseInt(posElems2[1]);
						String genes = elems[2].replaceAll("\"", "");
						String[] geneElems = Strings.whitespace.split(genes);

						if (elems[0].contains("MHC")) {
							Feature feat = new Feature();
							feat.setChromosome(chr);
							feat.setName("MHC / HLA");
							feat.setStart(0);
							feat.setStop(1000);
							regions.add(feat);
							System.out.println(gfffile + "\t MHC regions");
						} else {
							for (String gene : geneElems) {
								if (gene.trim().length() > 0) {
									Feature feat = new Feature();
									feat.setChromosome(chr);
									feat.setName(gene);
									feat.setStart(0);
									feat.setStop(1000);
									regions.add(feat);
								}
							}
						}

					}
					elems = tf.readLineElems(TextFile.comma);
				}
				tf.close();
				associatedLociPerDisease.add(regions);
				System.out.println(gfffile + "\t" + regions.size() + " regions");
			}
		}


		double[][][] data = new double[associatedLociPerDisease.size()][allChr.length][];
//		String[] catnames = new String[allChr.length];
//		String[] subcatnames = new String[];

		String[][] subcatnames = new String[allChr.length][];
		for (int chrn = 0; chrn < allChr.length; chrn++) {
			Chromosome chr = allChr[chrn];

			if (chr.getNumber() > 0 && chr.getNumber() < 23) {

				// make list of associated genes in set
				HashSet<String> geneNames = new HashSet<>();
				ArrayList<Feature> geneObj = new ArrayList<>();


				for (int dataset = 0; dataset < associatedLociPerDisease.size(); dataset++) {
					ArrayList<Feature> lociForDs = associatedLociPerDisease.get(dataset);
					for (Feature feature : lociForDs) {
						if (feature.getChromosome().equals(chr)) {
							String name = feature.getName();
							if (name == null || name.length() == 0) {
								name = feature.toString();
							}

							if (!geneNames.contains(name)) {
								geneNames.add(name);
								geneObj.add(feature);
							}
						}
					}
				}


				if (!geneObj.isEmpty()) {
					if (sortstuff) {
						Collections.sort(geneObj, new FeatureComparator());
					}
					HashMap<Feature, Integer> featMap = new HashMap<Feature, Integer>();
					int gctr = 0;
					subcatnames[chrn] = new String[geneObj.size()];
					for (Feature f : geneObj) {
						featMap.put(f, gctr);
						subcatnames[chrn][gctr] = f.getName();
						gctr++;
					}


					for (int dataset = 0; dataset < associatedLociPerDisease.size(); dataset++) {
						data[dataset][chrn] = new double[featMap.size()];
						ArrayList<Feature> lociForDs = associatedLociPerDisease.get(dataset);
						for (Feature feature : lociForDs) {
							if (feature.getChromosome().equals(chr)) {
								Integer loc = featMap.get(feature);
								if (loc != null) {
									data[dataset][chrn][loc] = 1;
								}
							}
						}
					}
				}
			}
		}

		// trim the dataset
		int missing = 0;
		boolean[] missingb = new boolean[allChr.length];
		for (int chrn = 0; chrn < allChr.length; chrn++) {
			for (int dataset = 0; dataset < associatedLociPerDisease.size(); dataset++) {
				if (data[dataset][chrn] == null || data[dataset][chrn].length == 0) {

					missingb[chrn] = true;
				}
			}
			if (missingb[chrn]) {
				missing++;
			}
		}

		System.out.println("missing: " + missing);

		int locusctr = 0;
		for (int chrn = 0; chrn < allChr.length; chrn++) {
			if (!missingb[chrn]) {

				locusctr += data[0][chrn].length;

			}
		}

		double[][] data2 = new double[associatedLociPerDisease.size()][locusctr];
//		String[] catNames = new String[allChr.length - missing];
		ArrayList<String> subcatnames2 = new ArrayList<String>();
		int chrctr = 0;
		int lctr = 0;
		ArrayList<Triple<Integer, Integer, String>> groups = new ArrayList<Triple<Integer, Integer, String>>();
		int prevlctr = 0;
		for (int chrn = 0; chrn < allChr.length; chrn++) {
			if (!missingb[chrn]) {

				for (int i = 0; i < data[0][chrn].length; i++) {
					for (int dataset = 0; dataset < associatedLociPerDisease.size(); dataset++) {
						data2[dataset][lctr] = data[dataset][chrn][i];
					}
					subcatnames2.add(subcatnames[chrn][i]);
					lctr++;
				}
				Triple<Integer, Integer, String> group = new Triple<Integer, Integer, String>(prevlctr, lctr, allChr[chrn].getName());
				groups.add(group);
				prevlctr = lctr;

				//				catNames[chrctr] = allChr[chrn].getName();
				System.out.println(allChr[chrn].toString());
				chrctr++;
			}
		}

		Grid grid = new Grid(750, 750, 1, 1, 100, 0);
		CircularHeatmapPanel panel = new CircularHeatmapPanel(1, 1);


		panel.setData(datasetnames, subcatnames2.toArray(new String[0]), data2);
		panel.setGroups(groups);
		grid.addPanel(panel);

		grid.draw(output);


	}
}
