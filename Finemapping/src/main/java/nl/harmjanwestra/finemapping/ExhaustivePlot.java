package nl.harmjanwestra.finemapping;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.panels.LDPanel;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 8/30/16.
 */
public class ExhaustivePlot {


	public static void main(String[] args) {
		// make the actual plot
//		Grid grid = new Grid(400, 400, 1, 1, 100, 100);
//
//		LDPanel panel = new LDPanel(1, 1);
//		Feature region = new Feature();
//		region.setStart(0);
//		region.setStop(10);
//		ArrayList<Pair<Integer, Integer>> positions = new ArrayList<>();
//		ArrayList<Double> pvals = new ArrayList<>();
//
//		positions.add(new Pair<Integer, Integer>(1, 1));
//		positions.add(new Pair<Integer, Integer>(2, 1));
//		positions.add(new Pair<Integer, Integer>(3, 1));
//		positions.add(new Pair<Integer, Integer>(2, 2));
//		positions.add(new Pair<Integer, Integer>(2, 3));
//		positions.add(new Pair<Integer, Integer>(3, 3));
//
//		pvals.add(1d);
//		pvals.add(0.2);
//		pvals.add(0.2);
//		pvals.add(1d);
//		pvals.add(0.5);
//		pvals.add(1d);
//
//
//		String output = "/Data/tmp/exhaustive/test.png";
//
//		panel.setData(region, positions, pvals);
//
//		try {
//			grid.addPanel(panel);
//			grid.draw(output);
//		} catch (IOException e) {
//			e.printStackTrace();
//		} catch (DocumentException e) {
//			e.printStackTrace();
//		}

		ExhaustivePlot p = new ExhaustivePlot();
		String output = "/Data/tmp/exhaustive/test.png";
		Feature region = new Feature();
		// Chr2_204446380-204816382
		region.setChromosome(Chromosome.TWO);
		region.setStart(204446380);
		region.setStop(204816382);
		String assocfile = "/Data/tmp/exhaustive/RA-assoc0.3-COSMO-chr2-pairwise.txt.gz";
		try {
			p.run(region, assocfile, output);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}

	}


	public static void run(Feature region, String assocfile, String output) throws IOException, DocumentException {

		String[] variants = new String[]{"", ""};

		TextFile tf2 = new TextFile(assocfile, TextFile.R);
		String headerln = tf2.readLine();
		String[] elems = tf2.readLineElems(TextFile.tab);
		double maxP = 0;
		ArrayList<Pair<Integer, Integer>> positions = new ArrayList<Pair<Integer, Integer>>();
		ArrayList<Double> pvals = new ArrayList<>();
		while (elems != null) {
//			lineIsWhatWereLookingFor = false;
			Feature s1 = new Feature();
			Chromosome chr = Chromosome.parseChr(elems[0]);
			Integer pos1 = Integer.parseInt(elems[1]);
			String snp1Id = elems[2];
			Integer pos2 = Integer.parseInt(elems[5]);
			String snp2Id = elems[6];

			SNPFeature snp1 = new SNPFeature();
			snp1.setChromosome(chr);
			snp1.setStart(pos1);
			snp1.setStop(pos1);
			snp1.setName(snp1Id);

			SNPFeature snp2 = new SNPFeature();
			snp2.setChromosome(chr);
			snp2.setStart(pos2);
			snp2.setStop(pos2);
			snp2.setName(snp2Id);


			if (region.overlaps(snp1) && region.overlaps(snp2)) {


				String snp1str = snp1.getChromosome().getNumber() + "_" + snp1.getStart() + "_" + snp1.getName();
				String snp2str = snp2.getChromosome().getNumber() + "_" + snp2.getStart() + "_" + snp2.getName();
//					System.out.println(snp1str);
//					System.exit(0);
				Double p = Double.parseDouble(elems[elems.length - 1]);
				if ((snp1str.equals(variants[0]) && snp2str.equals(variants[1]))
						|| (snp1str.equals(variants[1]) && snp2str.equals(variants[0]))) {
					System.out.println("found it.");

				}

				if (p > maxP) {
					maxP = p;
				}

				positions.add(new Pair<>(snp1.getStart(), snp2.getStart()));
				pvals.add(p);


			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();

		System.out.println(pvals.size()+" pvals");


		ArrayList<Double> pvals2 = new ArrayList<>();
		for (int i = 0; i < pvals.size(); i++) {
			pvals2.add(pvals.get(i) / maxP);
		}

		// make the actual plot
		Grid grid = new Grid(400, 400, 1, 1, 100, 100);

		LDPanel panel = new LDPanel(1, 1);
		panel.setData(region, positions, pvals2);
		grid.addPanel(panel);
		panel.scaleQuadratic(true);
		grid.draw(output + ".png");
	}


}
