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
 * Created by hwestra on 11/3/16.
 */
public class HaplotypeGroupPlot {


	public void run(ArrayList<Pair<HaplotypeGroup,
			AssociationResult>> allResults,
					ArrayList<VCFVariant> variants,
					String outputfile) throws IOException, DocumentException {
		int margin = 100;

		int blockSize = 10;
		int yMarginBetweenHaps = 2;
		int marginBetweenHaps = 2;
		int marginBetweenConditions = 20;

		int witdhOfForestplot = 150;

		int maxNrVariants = variants.size();
		int maxNrHaplotypes = 0;
		int maxGroup0 = 0;
		int maxGroup1 = 0;
		double minOR = Double.MAX_VALUE;
		double maxOr = 0;


		for (Pair<HaplotypeGroup, AssociationResult> p : allResults) {
			AssociationResult r = p.getRight();
			if (r.getConfHi()[0] > maxOr) {
				maxOr = r.getConfHi()[0];
			}
			if (r.getConfLo()[0] < minOR) {
				minOR = r.getConfLo()[0];
			}
			HaplotypeGroup g = p.getLeft();


			int nrHaps = g.getGroup0().size() + g.getGroup1().size();
			if (g.getGroup0().size() > maxGroup0) {
				maxGroup0 = g.getGroup0().size();
			}
			if (g.getGroup0().size() > maxGroup1) {
				maxGroup1 = g.getGroup1().size();
			}
			if (nrHaps > maxNrHaplotypes) {
				maxNrHaplotypes = nrHaps;
			}
		}

		Range range = new Range(minOR, 0, maxOr, 1);

		Collections.sort(allResults, new HaplogroupComparator());

		// calculate size of plot
		int hapWidth = blockSize * variants.size();

		int widthLeft = maxGroup0 * hapWidth + ((maxGroup0 - 1) * marginBetweenHaps);
		int widthRight = maxGroup1 * hapWidth + ((maxGroup1 - 1) * marginBetweenHaps);


		HashMap<VCFVariant, Integer> variantToInt = new HashMap<VCFVariant, Integer>();
		for (int v = 0; v < variants.size(); v++) {
			variantToInt.put(variants.get(v), v);
		}

		int totalHapWidth = (widthLeft + marginBetweenConditions + widthRight);
		int totalHapHeight = (allResults.size() * blockSize + ((allResults.size() - 1) * yMarginBetweenHaps));
		int totalHeight = (margin * 2) + totalHapHeight;
		int totalWidth = (margin * 2) + marginBetweenConditions + witdhOfForestplot;

		DefaultGraphics graphics = new DefaultGraphics(outputfile, totalWidth, totalHeight);
		Graphics2D g2d = graphics.getG2d();


		int midpoint = margin + (totalHapWidth / 2);
		int halfmarginbetweenconditions = marginBetweenConditions / 2;

		DefaultTheme theme = new DefaultTheme();
		g2d.setFont(theme.getSmallFont());


		Color colorA1 = new Color(10, 150, 181);
		Color colorA2 = new Color(100, 100, 100);
		Color colorAX = new Color(224, 224, 224);
		Color colorL = colorA2;

		// draw forestplot line
		int forestPlotStartX = margin + totalHapWidth + marginBetweenConditions;
		g2d.setColor(colorL);
		g2d.drawLine(forestPlotStartX, margin - 10, forestPlotStartX + witdhOfForestplot, margin - 10);
		double perc = range.getRelativePositionX(1d);
		int forestmidpoint = forestPlotStartX + (int) Math.floor(perc * witdhOfForestplot);
		Stroke prevStroke = g2d.getStroke();
		g2d.setStroke(theme.getStrokeDashed());
		g2d.drawLine(forestmidpoint, margin - 15, forestmidpoint, margin + totalHapHeight);
		g2d.setStroke(prevStroke);

		DecimalFormat format = new DecimalFormat("#.##");
		g2d.drawString("" + format.format(range.getMinX()), forestPlotStartX, margin - 20);
		g2d.drawString("" + 1, forestmidpoint, margin - 20);
		g2d.drawString("" + format.format(range.getMaxX()), forestPlotStartX + witdhOfForestplot, margin - 20);


		for (int i = 0; i < allResults.size(); i++) {
			Pair<HaplotypeGroup, AssociationResult> p = allResults.get(i);
			HaplotypeGroup g = p.getLeft();
			AssociationResult a = p.getRight();

			ArrayList<BitVector> b1 = g.getGroup0();
			ArrayList<BitVector> b2 = g.getGroup1();

			ArrayList<VCFVariant> hapvars = g.getVariants();
			ArrayList<Integer> varIds = new ArrayList<>();

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

			int y = margin + (i * blockSize + (i * yMarginBetweenHaps));
			for (int k = 0; k < b1.size(); k++) {
				BitVector v1 = b1.get(k);
				int x1 = midpoint - halfmarginbetweenconditions - (k * hapWidth) - (k * marginBetweenHaps);
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
			}

			for (int k = 0; k < b2.size(); k++) {
				BitVector v1 = b2.get(k);
				int x1 = midpoint + halfmarginbetweenconditions + (k * hapWidth) + (k * marginBetweenHaps);

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
			}

			double or = a.getORs()[0];
			double orhi = a.getConfHi()[0];
			double orlo = a.getConfLo()[0];

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
			g2d.fillOval(xc, y2 - halfblock, blockSize, blockSize);

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

	public void drawRotate(Graphics2D g2d, double x, double y, int angle, String text) {
		g2d.translate((float) x, (float) y);
		g2d.rotate(Math.toRadians(angle));
		g2d.drawString(text, 0, 0);
		g2d.rotate(-Math.toRadians(angle));
		g2d.translate(-(float) x, -(float) y);
	}

}

class HaplogroupComparator implements Comparator<Pair<HaplotypeGroup, AssociationResult>> {

	@Override
	public int compare(Pair<HaplotypeGroup, AssociationResult> o1, Pair<HaplotypeGroup, AssociationResult> o2) {
		AssociationResult r1 = o1.getRight();
		AssociationResult r2 = o2.getRight();

		double d1 = r1.getORs()[0];
		double d2 = r2.getORs()[0];

		if (d1 > d2) {
			return -1;
		} else if (d1 < d2) {
			return 1;
		} else {
			return 0;
		}


	}
}
