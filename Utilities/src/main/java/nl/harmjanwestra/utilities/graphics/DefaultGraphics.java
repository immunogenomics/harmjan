package nl.harmjanwestra.utilities.graphics;

import com.lowagie.text.Document;
import com.lowagie.text.DocumentException;
import com.lowagie.text.pdf.PdfWriter;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Locale;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * @author hwestra
 */
public class DefaultGraphics {

	protected Output output = Output.PDF;
	protected BufferedImage bi;
	protected Graphics2D g2d;
	protected Document document;
	protected PdfWriter writer;
	protected int figureWidth;
	protected int figureHeight;
	protected String outputFileName;
	private Locale defaultLocale;



	public enum Output {

		PDF, PNG
	}

	;

	protected static final Font LARGE_FONT = new Font("Verdana", Font.PLAIN, 14);
	protected static final Font LARGE_FONT_BOLD = new Font("Verdana", Font.BOLD, 14);
	protected static final Font SMALL_FONT = new Font("Verdana", Font.PLAIN, 10);
	protected static final Font SMALL_FONT_BOLD = new Font("Verdana", Font.BOLD, 10);

	protected static final Stroke dashed = new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 0, new float[]{4}, 0);
	protected static final Stroke line2pt = new BasicStroke(2, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
	protected static final Stroke line = new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
	protected com.lowagie.text.pdf.PdfContentByte cb = null;

	protected DefaultGraphics(){

	}

	protected DefaultGraphics(String outputFileName, int width, int height) throws FileNotFoundException, DocumentException {
		this.initializePlot(outputFileName,width,height);
	}

	protected void initializePlot(String outputFileName, int width, int height) throws DocumentException, FileNotFoundException {
		defaultLocale = Locale.getDefault();
		Locale.setDefault(Locale.US);
		// set up Graphics2D depending on required format using iText in case PDF
		g2d = null;
		document = null;
		writer = null;
		bi = null;
		this.outputFileName = outputFileName;
		bi = new BufferedImage(1, 1, BufferedImage.TYPE_INT_RGB);
		g2d = bi.createGraphics();

		// initialize plot
		if (output == Output.PDF) {
			com.lowagie.text.Rectangle rectangle = new com.lowagie.text.Rectangle(width, height);
			document = new com.lowagie.text.Document(rectangle);
			writer = com.lowagie.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(outputFileName));

			document.open();
			cb = writer.getDirectContent();
			cb.saveState();

			g2d = cb.createGraphics(width, height);
		} else {
			bi = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
			g2d = bi.createGraphics();
		}

		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setColor(Color.white);
		g2d.fillRect(0, 0, width, height);
		g2d.setStroke(line);
	}

	public void close() throws IOException {
		// dispose
		g2d.dispose();
		if (output == Output.PDF) {
			g2d.dispose();
			cb.restoreState();
			document.close();
			writer.close();
		} else {
			bi.flush();
			ImageIO.write(bi, output.toString().toLowerCase(), new File(outputFileName));
		}

		Locale.setDefault(defaultLocale);
	}
}
