/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.visualizationtool;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import umcg.genetica.containers.Exon;
import umcg.genetica.containers.Gene;
import umcg.genetica.containers.Transcript;
import umcg.genetica.ensembl.Features;
import umcg.genetica.gwas.Dependifier;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.ChrAnnotation;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.geom.QuadCurve2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author Harm-Jan
 */
public class GosiaViz {

	private Viz.Output output = Viz.Output.PNG;

	public enum Output {

		PDF, PNG
	}

	;
	private static final Font LARGE_FONT = new Font("Verdana", Font.PLAIN, 14);
	private static final Font LARGE_FONT_BOLD = new Font("Verdana", Font.BOLD, 14);
	private static final Font SMALL_FONT = new Font("Verdana", Font.PLAIN, 10);
	private static final int upMaxR = 108;
	private static final int upMaxG = 189;
	private static final int upMaxB = 69;

	private static final int downMaxR = 0;
	private static final int downMaxG = 174;
	private static final int downMaxB = 239;

	private static final float downMaxH = 196f / 360f;
	private static final float downMaxS = 1f;
	private static final float downMaxBr = 0.93f;

	private static final float upMaxH = 100f / 360f;
	private static final float upMaxS = 0.63f;
	private static final float upMaxBr = 0.74f;

	private static final Color red1Color = new Color(242, 101, 94, 50);
	private static final Color red2Color = new Color(242, 101, 94, 25);
	private static final Color gray1Color = new Color(237, 237, 237);
	private static final Color gray2Color = new Color(247, 247, 247);
	private static final Logger LOGGER = Logger.getLogger(Viz.class.getName());
	private static final Stroke dashed = new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 0, new float[]{4}, 0);
	private static final Stroke line2pt = new BasicStroke(2, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
	private static final Stroke line = new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);

	public static void main(String[] args) {
		try {

//            if (args.length < 5) {
//                System.out.println("Usage: GosiaViz.jar gwasCatalog trait features bedfiledir outdir");
//            } else {
			GosiaViz vis = new GosiaViz();
			double minP = 5E-8;
			String gwasCatalog = "C:\\Work\\DNASe\\2014-05-19-GWASCatalog.txt";
			String loci = "";
			String trait = "Rheumatoid arthritis";
			String features = "C:\\Work\\Ensembl\\structures.txt";
			String bedFileDir = "C:\\Work\\DNASe\\";
			String outdir = "C:\\Work\\test\\";
			String hapmap = "C:\\Work\\GGD\\HapMap2r24-CEU\\";

			vis.run(gwasCatalog, trait, minP, features, bedFileDir, hapmap, outdir);
//            }
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException ex) {
			Logger.getLogger(GosiaViz.class.getName()).log(Level.SEVERE, null, ex);
		}
	}

	public void run(String gwasCatalogLoc, String trait, double gwasMinimalPval, String features, String bedFileDir, String hapmap, String outdir) throws IOException, DocumentException {
		GWASCatalog catalog = new GWASCatalog(gwasCatalogLoc);
		GWASTrait[] traits = catalog.getTraitsForCertainKey(trait);

		String[] lsof = Gpio.getListOfFiles(bedFileDir);
		ArrayList<String> bedFiles = new ArrayList<String>();
		for (String file : lsof) {
			if (file.endsWith(".bed.gz")) {
				bedFiles.add(file);
				System.out.println("Found bedfile: " + file);
			}
		}
		System.out.println(bedFiles.size() + " bedFiles in total.");

		// load gwas catalog
		HashSet<GWASSNP> allSNPs = new HashSet<GWASSNP>();
		for (GWASTrait t : traits) {
			allSNPs.addAll(Arrays.asList(t.getSNPs(gwasMinimalPval)));
		}
		outdir = Gpio.formatAsDirectory(outdir);
		Gpio.createDir(outdir);
		// load ensembl features
		Features feat = new Features();
		feat.loadAnnotation(features);
		HashMap<String, Gene> genHash = feat.getGeneHash();

		// iterate SNPs
		int maxWindowSize = 1000000;
		int bpMargin = 50;
		int figureMargin = 350;
		int exonPlotSize = 20;
		int distPlotSize = 20;

		int distplotMargin = 5;
		int maxPlotWidth = 1000;
		int ticksEveryBp = 1000;
		Set<Entry<String, Gene>> geneSets = genHash.entrySet();

		TriTyperGenotypeData ds = new TriTyperGenotypeData(hapmap);
		SNPLoader loader = ds.createSNPLoader();
		Dependifier proxyfinder = new Dependifier(ds, loader);

		for (GWASSNP leadSNP : allSNPs) {
			// get the nearest genes
			byte chr = leadSNP.getChr();
			int pos = leadSNP.getPosition();
			boolean plot = true;
			String snp = leadSNP.getName();

			HashSet<String> proxies = proxyfinder.findProxiesForSNP(snp, 0.8, 1000000);
			int maxSNPPos = pos;
			int minSNPPos = pos;

			for (String proxy : proxies) {
				Integer id = ds.getSnpToSNPId().get(proxy);

				int snpPos = ds.getChrPos(id);
				if (snpPos > maxSNPPos) {
					maxSNPPos = snpPos;
				}
				if (snpPos < minSNPPos) {
					minSNPPos = snpPos;
				}
			}

			if (plot) {

				System.out.println(leadSNP.getName() + "\t" + chr + "\t" + pos);

				if (chr > 1 && pos > 1) {
					int minimumDist = Integer.MAX_VALUE;
					Gene closestGene = null;
					for (Entry<String, Gene> geneSet : geneSets) {
						Gene gene = geneSet.getValue();

						if (gene.getParentChromosome().getName().equals("" + chr)) {

							int distance = Math.abs(gene.getStart() - pos);
							if (distance < minimumDist) {
								minimumDist = distance;
								closestGene = gene;
							}
						}
					}

					if (closestGene == null) {
						System.out.println("No gene found");
					} else {
						System.out.println("Closest gene: " + minimumDist + " - " + closestGene.getName());

						// plot the thing
						// determine whether the snp is within the gene....
						boolean isInGene = false;
						if (pos >= closestGene.getStart() && pos <= closestGene.getEnd()) {
							isInGene = true;
						}
						boolean upstream = true;
						if (pos > closestGene.getEnd()) {
							upstream = false;
						}

						// set window size
						int windowStart = minSNPPos - bpMargin;
						int windowStop = maxSNPPos + bpMargin;

						// make the windowsize a fraction of the maxPlotWidth
						int windowSizeInBp = windowStop - windowStart;
						if (windowSizeInBp < maxPlotWidth) {

							int remainder = maxPlotWidth - windowSizeInBp;
							windowStart -= (remainder / 2);
							windowSizeInBp = windowStop - windowStart;
							remainder = maxPlotWidth - windowSizeInBp;
							windowStop += remainder;
						} else if (windowSizeInBp > maxPlotWidth) {
							// adjust start and stop for maxWindowSize
							int remainder = maxPlotWidth - (windowSizeInBp % maxPlotWidth);
							windowStart -= (remainder / 2);
							windowSizeInBp = windowStop - windowStart;
							remainder = maxPlotWidth - (windowSizeInBp % maxPlotWidth);
							windowStop += remainder;
						}

						windowSizeInBp = windowStop - windowStart; // windowsize should now be %maxPlotWidth.
						System.out.println(windowSizeInBp);
						// determine whether there are other genes in the window
						HashSet<Transcript> otherGenes = new HashSet<Transcript>();
						for (Entry<String, Gene> geneSet : geneSets) {
							Gene gene = geneSet.getValue();
							if (gene.getParentChromosome().getName().equals(closestGene.getParentChromosome().getName())) {
//                                

								if (gene.getStart() >= windowStart && gene.getEnd() <= windowStop) {
									otherGenes.addAll(gene.getTranscripts().values());
								} else if (gene.getStart() <= windowStart && gene.getEnd() >= windowStart) {
									otherGenes.addAll(gene.getTranscripts().values());
								} else if (gene.getStart() <= windowStop && gene.getEnd() >= windowStop) {
									otherGenes.addAll(gene.getTranscripts().values());
								}
							}

						}

						// initialize figure
						int height = figureMargin * 2
								+ otherGenes.size() * exonPlotSize
								+ distPlotSize
								+ bedFiles.size() * (distPlotSize + distplotMargin);
						int width = figureMargin * 2 + maxPlotWidth;

						Locale defaultLocale = Locale.getDefault();
						Locale.setDefault(Locale.US);
						// set up Graphics2D depending on required format using iText in case PDF
						Graphics2D g2d = null;
						com.itextpdf.text.Document document = null;
						com.itextpdf.text.pdf.PdfWriter writer = null;
						BufferedImage bi = null;

						bi = new BufferedImage(1, 1, BufferedImage.TYPE_INT_RGB);
						g2d = bi.createGraphics();

						g2d.setFont(LARGE_FONT);
						FontMetrics fontmetrics = g2d.getFontMetrics();

						int fontheight = fontmetrics.getHeight();

						String outputFileName = outdir + leadSNP.getName() + "-" + closestGene.getName() + ".png";

						// initialize plot
						com.itextpdf.text.pdf.PdfContentByte cb = null;
						if (output == Viz.Output.PDF) {
							com.itextpdf.text.Rectangle rectangle = new com.itextpdf.text.Rectangle(width, height);
							document = new com.itextpdf.text.Document(rectangle);
							writer = com.itextpdf.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(outputFileName));

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

						// title, subtitle
						int bpPerPixel = windowSizeInBp / maxPlotWidth;
						g2d.setColor(Color.black);
						g2d.drawString("SNP: " + leadSNP.getName() + " Chr. " + ChrAnnotation.parseByte(chr) + " Pos: " + pos, figureMargin, figureMargin / 2);

						g2d.drawString("Showing: " + windowStart + " - " + windowStop + " (" + windowSizeInBp + "bp) " + bpPerPixel + " bp/pixel", figureMargin, (figureMargin / 2) + 15);
//                    

						g2d.drawLine(figureMargin, figureMargin, width - figureMargin, figureMargin);

						// scale ticks..
						int nrTicks = windowSizeInBp / ticksEveryBp;
						while (nrTicks > 10) {
							ticksEveryBp *= 10;
							nrTicks = windowSizeInBp / ticksEveryBp;
						}

						Stroke s = g2d.getStroke();
						g2d.setStroke(dashed);
						int tickStart = windowStart;
						int remainderTickStart = ticksEveryBp - windowStart % ticksEveryBp;
						tickStart = windowStart + remainderTickStart;
						while (tickStart < windowStop) {
							int posX = figureMargin + (tickStart - windowStart) / bpPerPixel;

							g2d.setColor(Color.black);
							g2d.drawString("" + tickStart, posX, figureMargin - 10);

							g2d.setColor(new Color(128, 128, 128, 64));
							g2d.drawLine(posX, figureMargin - 5, posX, height - figureMargin + 5);
							tickStart += ticksEveryBp;
						}
						g2d.setStroke(s);
						// plot the exons themselves
						// -- plot a line
						// -- plot arcs for the gaps between exons, unless exons are close together

						int y = 0;
						int exonPlotY = figureMargin;
						g2d.setColor(new Color(0, 0, 0, 128));
						for (Transcript t : otherGenes) {
							int strand = t.getStrand();
							Exon[] exons = t.getExonsRanked();
							Exon previousExon = null;
							int posY = exonPlotY + (y * exonPlotSize) + exonPlotSize / 2;
							int transcriptEnd = t.getEnd();

							if (transcriptEnd > windowStop) {
								transcriptEnd = windowStop;
							}
							int maxXPos = figureMargin + ((transcriptEnd - windowStart) / bpPerPixel) + 5;
							int exonPlotHeight = 5;
							g2d.drawString(t.getParentGene().getAnnotation() + " - " + t.getName(), figureMargin + maxPlotWidth + 5, posY + exonPlotHeight);
							for (Exon e : exons) {
								int exonStart = e.getStart();
								int exonStop = e.getEnd();

								if ((exonStart < windowStart && exonStop < windowStart)
										|| (exonStart > windowStop && exonStop > windowStop)) {
									// don't draw box
								} else {
									boolean overlapsEnd = false;
									boolean overlapsStart = false;
									if (exonStart < windowStart) {
										exonStart = windowStart;
										overlapsStart = true;
									}
									if (exonStop > windowStop) {
										exonStop = windowStop;
										overlapsEnd = true;
									}

									// draw box
									int posX = figureMargin + (exonStart - windowStart) / bpPerPixel;
									int len = (exonStop - exonStart) / bpPerPixel;
									if (len == 0) {
										len = 1;
									}

									g2d.fillRect(posX, posY, len, exonPlotHeight);
									if (previousExon != null) {

										if (e.getEnd() < previousExon.getStart()) {
											/* draw arc
											 curr    prev
                                             ||-------||
                                             */
											int previousExonStart = previousExon.getStart();
											int curveXStart = posX + len;
											int curveXStop = figureMargin + (previousExonStart - windowStart) / bpPerPixel;
											int curveYStart = posY;
											int xControlpoint = (curveXStart + curveXStop) / 2;
											int yControlpoint = exonPlotY + (y * exonPlotSize);

											QuadCurve2D curve = new QuadCurve2D.Double();

											curve.setCurve(curveXStart,
													curveYStart,
													xControlpoint,
													yControlpoint,
													curveXStop,
													curveYStart);
											g2d.draw(curve);

										} else {
                                            /* draw arc
                                             prev    curr
                                             ||-------||
                                             */

											int previousExonStop = previousExon.getEnd();
											int curveXStart = posX;
											int curveXStop = figureMargin + (previousExonStop - windowStart) / bpPerPixel;
											int curveYStart = posY;
											int xControlpoint = (curveXStart + curveXStop) / 2;
											int yControlpoint = exonPlotY + (y * exonPlotSize);

											QuadCurve2D curve = new QuadCurve2D.Double();
											curve.setCurve(curveXStart,
													curveYStart,
													xControlpoint,
													yControlpoint,
													curveXStop,
													curveYStart);
											g2d.draw(curve);
										}

									}
									previousExon = e;
								}

							}

							y++;
						}
						// -- plot boxes for the exons
						// -- plot the gene name

						int distPlotStartY = exonPlotY + (y * exonPlotSize) + distplotMargin + distPlotSize;

						// load the BED files
						BedFileReader reader = new BedFileReader();
						Color lineColor = new Color(0, 0, 0, 128);
						y = 0;
//                        for (String path : bedFiles) {
//                            // load the track
//                            String chrStr = "chr" + ChrAnnotation.parseByte(chr);
//                            Track t = reader.readAsTrack(bedFileDir + path, path, Chromosome.parseChr(chrStr), windowStart, windowStop, true);
//                            // bin the peaks
//
//                            Set<Feature> featureSet = t.getFeatureSet(Chromosome.parseChr(chrStr), windowStart, windowStop);
//
//                            int yPos = distPlotStartY + y * (distPlotSize + distplotMargin);
//                            System.out.println(featureSet.size());
//                            Stroke defaultStroke = g2d.getStroke();
//                            for (Feature f : featureSet) {
//                                BedFileFeature f2 = (BedFileFeature) f;
//                                g2d.setColor(new Color(0, 125, 255, 64));
//                                int start =  f2.getStart();
//                                int stop =  f2.getStop();
//
//                                if (start < windowStart) {
//                                    start = windowStart;
//                                }
//                                if (stop > windowStop) {
//                                    stop = windowStop;
//                                }
//                                int x0 = figureMargin + (start - windowStart) / bpPerPixel;
//                                int stopX = figureMargin + (stop - windowStart) / bpPerPixel;
//
////                                double peakQ =  f2.getPeakQ();
////                                if (peakQ > 100) {
////                                    peakQ = 100;
////                                }
////                                int peakHeight = (int) Math.ceil((peakQ / 100) * distPlotSize);
////
////                                g2d.fillRect(x0, yPos + distPlotSize - peakHeight, stopX - x0, peakHeight);
////
////                                int peakPos =  f2.getPeakPos();
////                                int peakPosX = figureMargin + (peakPos - windowStart) / bpPerPixel;
//                                g2d.setColor(new Color(0, 125, 255, 255));
//                                g2d.setStroke(line2pt);
////                                g2d.drawLine(peakPosX, yPos + distPlotSize - peakHeight, peakPosX, yPos + distPlotSize);
//                                System.out.println(stopX - x0);
//                                g2d.setStroke(defaultStroke);
//                            }
//                            g2d.setStroke(defaultStroke);
//                            g2d.setColor(lineColor);
//                            g2d.drawLine(figureMargin, yPos + distPlotSize, (width - figureMargin), yPos + distPlotSize);
//                            // plot the bins
//                            g2d.drawString(path, (width - figureMargin) + 5, yPos + distPlotSize);
//                            y++;
//                        }

						// plot a red line for the position of the variant
						g2d.setColor(new Color(255, 125, 0, 255));
						g2d.setStroke(line2pt);
						int snpXPos = figureMargin + (pos - windowStart) / bpPerPixel;
						g2d.drawLine(snpXPos, figureMargin - 45, snpXPos, (height - figureMargin) + 5);
						g2d.drawString(leadSNP.getName(), snpXPos, figureMargin - 50);

						// plot the proxyfinder
						g2d.setColor(new Color(255, 125, 0, 64));
						g2d.setStroke(dashed);
						for (String proxy : proxies) {
							Integer id = ds.getSnpToSNPId().get(proxy);
							int snpPos = ds.getChrPos(id);
							snpXPos = figureMargin + (snpPos - windowStart) / bpPerPixel;
							g2d.drawLine(snpXPos, figureMargin - 45, snpXPos, (height - figureMargin) + 5);
							// g2d.drawString(proxy, snpXPos, figureMargin - 25);
						}

						// dispose
						g2d.dispose();
						if (output == Viz.Output.PDF) {
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
			}
		}
		loader.close();
	}

}
