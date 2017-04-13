package nl.harmjanwestra.finemappingtools.gwas;

import cern.colt.matrix.tbit.BitVector;
import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.graphics.DefaultGraphics;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.themes.DefaultTheme;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;

import java.awt.*;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

/**
 * Created by hwestra on 11/6/16.
 */
public class HaplotypeMultivariatePlot {

	Color colorA1 = new Color(10, 150, 181);
	Color colorA2 = new Color(100, 100, 100);
	Color colorAX = new Color(224, 224, 224);
	Color colorL = colorA2;


	class PairSorter implements Comparator<Pair<Integer, Double>> {

		@Override
		public int compare(Pair<Integer, Double> o1, Pair<Integer, Double> o2) {
			return o1.getRight().compareTo(o2.getRight());
		}
	}

	public void run(
			HaplotypeDataset dataset,
			ArrayList<BitVector> testedHaps,
			AssociationResult associationResult,
			Pair<Double, Double> minmax,
			String outputfile) throws IOException, DocumentException {
		int margin = 100;

		int blockSize = 10;
		int yMarginBetweenHaps = 2;
		int marginBetweenHaps = 2;
		int marginBetweenConditions = 20;

		int witdhOfForestplot = 150;


		int maxNrHaplotypes = 0;
		int maxGroup0 = 0;
		double minOR = Double.MAX_VALUE;
		double maxOr = 0;


		double[] orhi = associationResult.getConfHi();
		double[] orlo = associationResult.getConfLo();

		for (double d : orhi) {
			if (d < minOR) {
				minOR = d;
			}
			if (d > maxOr) {
				maxOr = d;
			}
		}

		for (double d : orlo) {
			if (d < minOR) {
				minOR = d;
			}
			if (d > maxOr) {
				maxOr = d;
			}
		}

		if (maxOr < 1) {
			maxOr = 1;
		}
		if (minOR > 1) {
			minOR = 1;
		}

		if (minmax != null) {
			minOR = minmax.getLeft();
			maxOr = minmax.getRight();
		}

		maxNrHaplotypes = testedHaps.size() + 1;
		ArrayList<VCFVariant> variants = dataset.getVariants();


		Range range = new Range(minOR, 0, maxOr, 1);

		// calculate size of plot
		int hapWidth = blockSize * variants.size();

		int widthLeft = maxGroup0 * hapWidth + ((maxGroup0 - 1) * marginBetweenHaps);


		int nrHaps = testedHaps.size() + 1;

		int totalHapWidth = widthLeft;
		int totalHapHeight = (nrHaps * blockSize + ((nrHaps - 1) * yMarginBetweenHaps));
		int totalHeight = (margin * 2) + totalHapHeight;
		int totalWidth = (margin * 2) + marginBetweenConditions + witdhOfForestplot;

		DefaultGraphics graphics = new DefaultGraphics(outputfile, totalWidth, totalHeight);
		Graphics2D g2d = graphics.getG2d();


		int midpoint = margin + (totalHapWidth / 2);
		int halfmarginbetweenconditions = marginBetweenConditions / 2;

		DefaultTheme theme = new DefaultTheme();
		g2d.setFont(theme.getSmallFont());


		// draw forestplot line
		int forestPlotStartX = margin + totalHapWidth + marginBetweenConditions;
		g2d.setColor(colorL);
		g2d.drawLine(forestPlotStartX, margin - 10, forestPlotStartX + witdhOfForestplot, margin - 10);
		double perc = range.getRelativePositionX(1d);
		int forestmidpoint = forestPlotStartX + (int) Math.floor(perc * witdhOfForestplot);
		Stroke prevStroke = g2d.getStroke();
		System.out.println(forestmidpoint + " forest midpoint");
		g2d.setStroke(theme.getStrokeDashed());
		g2d.drawLine(forestmidpoint, margin - 15, forestmidpoint, margin + totalHapHeight);
		g2d.setStroke(prevStroke);

		DecimalFormat format = new DecimalFormat("#.##");
		g2d.drawString("" + format.format(range.getMinX()), forestPlotStartX, margin - 20);
		g2d.drawString("" + 1, forestmidpoint, margin - 20);
		g2d.drawString("" + format.format(range.getMaxX()), forestPlotStartX + witdhOfForestplot, margin - 20);


		BitVector reference = dataset.getReferenceHaplotype();
		int referenceId = dataset.getHaplotypeId(reference);

		int hapctr = 0;
		double[] freqsCases = dataset.getHaplotypeFrequenciesCases();
		double[] freqsControls = dataset.getHaplotypeFrequenciesControls();

		drawHap(variants,
				reference,
				hapctr,
				1d,
				1d,
				1d,
				g2d,
				blockSize,
				margin,
				yMarginBetweenHaps,
				midpoint,
				halfmarginbetweenconditions,
				hapWidth,
				range,
				forestPlotStartX,
				witdhOfForestplot,
				freqsCases[referenceId],
				freqsControls[referenceId]
		);


		ArrayList<Pair<Integer, Double>> d = new ArrayList<Pair<Integer, Double>>();
		for (int q = 0; q < testedHaps.size(); q++) {
			d.add(new Pair<Integer, Double>(q, associationResult.getORs()[q]));
		}

		Collections.sort(d, new PairSorter());

		for (int z = 0; z < d.size(); z++) {
			int q = d.get(z).getLeft();
			int hapid = dataset.getHaplotypeId(testedHaps.get(q));
			drawHap(variants,
					testedHaps.get(q),
					z + 1,
					associationResult.getORs()[q],
					associationResult.getConfHi()[q],
					associationResult.getConfLo()[q],
					g2d,
					blockSize,
					margin,
					yMarginBetweenHaps,
					midpoint,
					halfmarginbetweenconditions,
					hapWidth,
					range,
					forestPlotStartX,
					witdhOfForestplot,
					freqsCases[hapid],
					freqsControls[hapid]);


		}


		// draw alleles
		int x0 = 10;
		int y0 = 30;
		int y1 = 30 + marginBetweenHaps + blockSize;

		for (int v = 0; v < variants.size(); v++) {

			int x = 10 + (v * blockSize);

			VCFVariant variant = variants.get(v);
			String allele1 = variant.getAlleles()[0];
			String allele2 = variant.getAlleles()[1];

			// draw boxes for Allele1
			g2d.setColor(colorA1);
			g2d.fillRect(x, y0, blockSize, blockSize);
			g2d.setColor(colorAX);
			g2d.drawRect(x, y0, blockSize, blockSize);

			g2d.setColor(colorA2);
			g2d.fillRect(x, y1, blockSize, blockSize);
			g2d.setColor(colorAX);
			g2d.drawRect(x, y1, blockSize, blockSize);

			g2d.setColor(colorL);
			g2d.drawString(allele1, x, y0 - 5);
			g2d.drawString(allele2, x, y1 + blockSize + 15);
			drawRotate(g2d, x, y1 + 30, 90, variant.getId());
		}
		g2d.setColor(colorL);
		g2d.drawRect(x0, y0, hapWidth, blockSize);
		g2d.drawRect(x0, y1, hapWidth, blockSize);


		graphics.close();

	}

	private void drawHap(
			ArrayList<VCFVariant> variants,
			BitVector hap,
			int hapctr,
			double or,
			double orhi,
			double orlo,
			Graphics2D g2d,
			int blockSize,
			int margin,
			int yMarginBetweenHaps,
			int midpoint,
			int halfmarginbetweenconditions,
			int hapWidth,
			Range range,
			int forestPlotStartX,
			int witdhOfForestplot,
			double freqsCase,
			double freqsControl) {

		ArrayList<VCFVariant> hapvars = variants;
		ArrayList<Integer> varIds = new ArrayList<>();

		HashMap<VCFVariant, Integer> variantToInt = new HashMap<VCFVariant, Integer>();
		for (int v = 0; v < variants.size(); v++) {
			variantToInt.put(variants.get(v), v);
		}

		HashSet<Integer> variantsPresent = new HashSet<Integer>();
		for (int v = 0; v < hapvars.size(); v++) {
			Integer id = variantToInt.get(hapvars.get(v));
			varIds.add(id);
			variantsPresent.add(id);
		}
		ArrayList<Integer> varIdsNotPresent = new ArrayList<>();
		for (int v = 0; v < variants.size(); v++) {
			if (!variantsPresent.contains(v)) {
				varIdsNotPresent.add(v);
			}
		}

		int y = margin + (hapctr * blockSize + (hapctr * yMarginBetweenHaps));

		BitVector v1 = hap;
		int x1 = midpoint - halfmarginbetweenconditions;
		x1 -= hapWidth;


		// draw blocks not present
		for (int j = 0; j < varIdsNotPresent.size(); j++) {
			Integer z = varIdsNotPresent.get(j);
			int x2 = x1 + (z * blockSize);
			g2d.setColor(colorAX);
			g2d.fillRect(x2, y, blockSize, blockSize);
			g2d.drawRect(x2, y, blockSize, blockSize);
		}

		for (int j = 0; j < v1.size(); j++) {
			Integer z = varIds.get(j);
			int x2 = x1 + (z * blockSize);
			if (v1.get(j)) {
				g2d.setColor(colorA2);
				g2d.fillRect(x2, y, blockSize, blockSize);
			} else {
				g2d.setColor(colorA1);
				g2d.fillRect(x2, y, blockSize, blockSize);
			}
			g2d.setColor(colorAX);
			g2d.drawRect(x2, y, blockSize, blockSize);

		}

		g2d.setColor(colorL);
		g2d.drawRect(x1, y, hapWidth, blockSize);

		double perc1 = range.getRelativePositionX(or);
		double perc2 = range.getRelativePositionX(orhi);
		double perc3 = range.getRelativePositionX(orlo);

		int y2 = y + (blockSize / 2);
		int xa = forestPlotStartX + (int) Math.floor(perc2 * witdhOfForestplot);
		int xb = forestPlotStartX + (int) Math.floor(perc3 * witdhOfForestplot);
		int xc = forestPlotStartX + (int) Math.floor(perc1 * witdhOfForestplot);

		g2d.setColor(colorL);
		g2d.drawLine(xa, y2, xb, y2);
		int halfblock = blockSize / 2;
		g2d.fillOval(xc - halfblock, y2 - halfblock, blockSize, blockSize);


		DecimalFormat format = new DecimalFormat("#.###");

		String t = format.format(freqsCase);
		int width = g2d.getFontMetrics().stringWidth(t);
		int xt = forestPlotStartX + witdhOfForestplot + 10;
		g2d.drawString(t, xt, y);
		String t2 = format.format(freqsControl);
		xt = forestPlotStartX + witdhOfForestplot + 10 + width + 10;
		g2d.drawString(t2, xt, y);


	}

	public void drawRotate(Graphics2D g2d, double x, double y, int angle, String text) {
		g2d.translate((float) x, (float) y);
		g2d.rotate(Math.toRadians(angle));
		g2d.drawString(text, 0, 0);
		g2d.rotate(-Math.toRadians(angle));
		g2d.translate(-(float) x, -(float) y);
	}
}
