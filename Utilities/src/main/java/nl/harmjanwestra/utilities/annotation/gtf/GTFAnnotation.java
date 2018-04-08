/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.utilities.annotation.gtf;

import nl.harmjanwestra.utilities.annotation.Annotation;
import nl.harmjanwestra.utilities.features.Exon;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.features.Transcript;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashMap;

/**
 * @author Harm-Jan
 */
public class GTFAnnotation extends Annotation {
	
	
	public static void main(String[] args) {
		String ensemblannotation = "D:\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz";
		try {
			GTFAnnotation g = new GTFAnnotation(ensemblannotation);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public GTFAnnotation(String annotation) throws IOException {
		this.annotationLocation = annotation;
		readAnnotation();
	}
	
	private void readAnnotation() throws IOException {
		
		System.out.println("Reading annotation: " + annotationLocation);
		
		TextFile tf = new TextFile(annotationLocation, TextFile.R);
		String ln = tf.readLine();
		
		strToGene = new HashMap<String, Gene>();
		
		HashMap<String, Transcript> strToTranscript = new HashMap<String, Transcript>();
		HashMap<String, Exon> strToExon = new HashMap<String, Exon>();

// this all assumes the path is sorted on genomic coordinates..
		while (ln != null) {
			if (ln.startsWith("#")) {
				// skip
			} else {
				GTFLine lineObj = new GTFLine(ln);
				
				
				if (lineObj.getType().equals("exon") || lineObj.getType().equals("start_codon") || lineObj.getType().equals("stop_codon") || lineObj.getType().equals("transcript")
						|| lineObj.getType().equals("gene")) {
					Gene currentGene = null;
					Transcript currentTranscript = null;
					
					String geneName = lineObj.getGeneId();
					String transcriptName = lineObj.getTranscriptId();
					Gene tmpGene = new Gene(geneName, lineObj.getChr(), lineObj.getStr());
					tmpGene.setGeneSymbol(lineObj.getGeneName());
					if (strToGene.containsKey(tmpGene.getName())) {
						// continue annotating transcripts and exons with this gene
//                System.out.println("getting gene: " + geneName);
						
						currentGene = strToGene.get(tmpGene.getName());
						
					} else {
//                System.out.println("new Gene obj");
						currentGene = tmpGene;
						strToGene.put(tmpGene.getName(), currentGene);
					}
					
					
					
					
					Exon e = new Exon(lineObj.getType(), lineObj.getChr(), lineObj.getStr(), currentGene, lineObj.getStart(), lineObj.getStop());
					e.setName(lineObj.getExonId());
					
					if (strToTranscript.containsKey(transcriptName)) {
						currentTranscript = strToTranscript.get(transcriptName);
						if (strToExon.containsKey(e.getName())) {
							// no work to be done
							e = strToExon.get(e.getName());
							e.addTranscript(currentTranscript);
							//System.out.println(e.toString());
							//System.out.println("Duplicate entry in path: " + lineObj.toString());
						} else {
							currentTranscript.addExon(e);
							e.addTranscript(currentTranscript);
							strToExon.put(e.getName(), e);
						}
					} else {
						currentTranscript = new Transcript(lineObj.getTranscriptId(), lineObj.getChr(), lineObj.getStr(), currentGene);
						currentGene.addTranscript(currentTranscript);
						strToTranscript.put(lineObj.getTranscriptId(), currentTranscript);
//                System.out.println(e.toString());
						if (strToExon.containsKey(e.getName())) {
							e = strToExon.get(e.getName());
							e.addTranscript(currentTranscript);
							currentTranscript.addExon(e);
						} else {
							currentTranscript.addExon(e);
							e.addTranscript(currentTranscript);
							strToExon.put(e.getName(), e);
						}
					}
					
				}
			}
			
			
			ln = tf.readLine();
		}
		tf.close();
		System.out.println(annotationLocation + ", Genes: " + strToGene.size() + "\tTranscripts: " + strToTranscript.size() + "\tExons: " + strToExon.size());
		
		// set the relative start and end positions of each gene and transcript
		genes = strToGene.values();
		genes.forEach(Gene::getBounds);
		
	}
	
	
}
