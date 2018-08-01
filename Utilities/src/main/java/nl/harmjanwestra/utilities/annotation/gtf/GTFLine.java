/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.utilities.annotation.gtf;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.enums.Strand;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.UnsupportedEncodingException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author Harm-Jan
 */
public class GTFLine {
	
	public String getExonId() {
		return exonId;
	}
	
	public String getTypeStr() {
		return typeStr;
	}
	
	private String exonId;
	private String attribStr;
	private Chromosome chr;
	private Strand str;
	private int start;
	private int stop;
	private Double score;
	private Double frame;
	private String geneId;
	private String geneName;
	private String transcriptId;
	private String tssId;
	private String typeStr;
	
	
	public GTFLine(String ln) {
		try {
			String[] elems = Strings.tab.split(ln);
			String sequenceStr = elems[0];
			String sourceStr = elems[1];
			typeStr = new String(elems[2].getBytes("UTF-8")).intern();
			String startStr = elems[3];
			String stopStr = elems[4];
			String scoreStr = elems[5];
			String strandStr = elems[6];
			String frameStr = elems[7];
			attribStr = elems[8];
			
			chr = Chromosome.parseChr(sequenceStr);
			str = Strand.parseStr(strandStr);
			start = Integer.parseInt(startStr);
			stop = Integer.parseInt(stopStr);
			score = null;
			frame = null;
			try {
				score = Double.parseDouble(scoreStr);
			} catch (NumberFormatException e) {
			
			}
			try {
				frame = Double.parseDouble(frameStr);
			} catch (NumberFormatException e) {
			
			}
			
			String[] attribElems = attribStr.split("; ");
			
			for (String attribElem : attribElems) {
				while (attribElem.startsWith(" ")) {
					attribElem = attribElem.substring(1);
				}
				String[] attribSubElems = attribElem.split(" ");
				String property = attribSubElems[0].toLowerCase().replaceAll(" ", "");
				String value = attribSubElems[1].replaceAll("\"", "");
				value = value.replaceAll(";", "");
				if (property.equals("gene_id")) {
					geneId = new String(value.getBytes("UTF-8")).intern();
				} else if (property.equals("gene_name")) {
					geneName = new String(value.getBytes("UTF-8")).intern();
				} else if (property.equals("transcript_id")) {
					transcriptId = new String(value.getBytes("UTF-8")).intern();
				} else if (property.equals("tss_id")) {
					tssId = new String(value.getBytes("UTF-8")).intern();
				} else if (property.equals("exon_id")) {
					exonId = new String(value.getBytes("UTF-8")).intern();
				}
			}
			if (geneName != null && tssId != null) {
				geneName += "_" + tssId;
				
			}
		} catch (UnsupportedEncodingException ex) {
			Logger.getLogger(GTFLine.class.getName()).log(Level.SEVERE, null, ex);
		}
		
	}
	
	public String getTranscriptId() {
		return transcriptId;
	}
	
	public String getTssId() {
		return tssId;
	}
	
	public String getAttribStr() {
		return attribStr;
	}
	
	public Chromosome getChr() {
		return chr;
	}
	
	public Strand getStr() {
		return str;
	}
	
	public int getStart() {
		return start;
	}
	
	public int getStop() {
		return stop;
	}
	
	public Double getScore() {
		return score;
	}
	
	public Double getFrame() {
		return frame;
	}
	
	public String getGeneId() {
		return geneId;
	}
	
	public String getGeneName() {
		return geneName;
	}
	
	public String getType() {
		return typeStr;
	}
	
	@Override
	public String toString() {
		return "GTFLine{" + "attribStr=" + attribStr + ", chr=" + chr + ", str=" + str + ", start=" + start + ", stop=" + stop + ", score=" + score + ", frame=" + frame + ", geneId=" + geneId + ", geneName=" + geneName + ", transcriptId=" + transcriptId + ", tssId=" + tssId + ", typeStr=" + typeStr + '}';
	}
	
}
