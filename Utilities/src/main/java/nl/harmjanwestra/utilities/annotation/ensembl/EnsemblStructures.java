package nl.harmjanwestra.utilities.annotation.ensembl;

import nl.harmjanwestra.utilities.annotation.Annotation;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.enums.Strand;
import nl.harmjanwestra.utilities.features.Exon;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.features.Transcript;
import nl.harmjanwestra.utilities.features.UTR;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashMap;

/**
 * Created by hwestra on 11/14/16.
 */
public class EnsemblStructures extends Annotation {


	public EnsemblStructures(String location) throws IOException {
		this.annotationLocation = location;
		loadAnnotation();
	}

	/*
		* 0	Ensembl Gene ID 1	Ensembl Transcript ID 2	Ensembl Protein ID
		*
		* 3	Chromosome Name 4	Gene Start (bp) 5	Gene End (bp) 6	Transcript Start (bp) 7	Transcript End (bp) 8	Strand 9	Associated Gene Name 10	Associated Gene DB
		* 11	5' UTR Start 12	5' UTR End 13	3' UTR Start 14	3' UTR End 15	CDS Length 16	Transcript count 17	Description 18	Gene Biotype 19	Ensembl Exon ID 20	Exon
		* Chr Start (bp) 21	Exon Chr End (bp) 22	Constitutive Exon 23	Exon Rank in Transcript 24	phase 25	cDNA coding start 26	Genomic coding start 27	cDNA coding
		* end 28	Genomic coding end 29	CDS End 30	CDS Start
		*/
	public void loadAnnotation() throws IOException {
		System.out.println("Loading sequence feature annotation from: " + annotationLocation);

		TextFile in = new TextFile(annotationLocation, TextFile.R);

		String[] linesplit = in.readLineElems(TextFile.tab);
		HashMap<String, Integer> elementToColId = new HashMap<String, Integer>();
		for (int i = 0; i < linesplit.length; i++) {
//			System.out.println(linesplit[i] + "\t" + i);
			elementToColId.put(linesplit[i], i);
		}

		HashMap<String, Gene> geneHash = new HashMap<String, Gene>();
		HashMap<String, Transcript> transcriptHash = new HashMap<String, Transcript>();
		HashMap<String, Exon> exonHash = new HashMap<String, Exon>();

		int lnCounter = 0;
		String[] elems = in.readLineElems(TextFile.tab);

		Integer colgeneId = elementToColId.get("Ensembl Gene ID");
		Integer coltranscriptId = elementToColId.get("Ensembl Transcript ID");
		Integer colprotein = elementToColId.get("Ensembl Protein ID");
		Integer colchromosome = elementToColId.get("Chromosome Name");
		Integer colgenestart = elementToColId.get("Gene Start (bp)");
		Integer colgenestop = elementToColId.get("Gene End (bp)");
		Integer coltranscriptstart = elementToColId.get("Transcript Start (bp)");
		Integer coltranscriptstop = elementToColId.get("Transcript End (bp)");
		Integer colstrand = elementToColId.get("Strand");
		Integer colgeneSymbol = elementToColId.get("Associated Gene Name");
		Integer colexon = elementToColId.get("Ensembl Exon ID");
		Integer colexonstart = elementToColId.get("Exon Chr Start (bp)");
		Integer colexonend = elementToColId.get("Exon Chr End (bp)");
		Integer colexonrankintranscript = elementToColId.get("Exon Rank in Transcript");

		Integer colUTR3start = elementToColId.get("3' UTR Start");
		Integer colUTR3stop = elementToColId.get("3' UTR End");
		Integer colUTR5start = elementToColId.get("5' UTR Start");
		Integer colUTR5stop = elementToColId.get("5' UTR End");


		while (elems != null) {

            /*
			 * Ensembl Gene ID Ensembl Transcript ID Ensembl Protein ID Chromosome Name Gene Start (bp) Gene End (bp) Transcript Start (bp) Transcript End (bp)
             * Strand Associated Gene Name Gene Biotype Description Transcript count CDS Length 3' UTR End 3' UTR Start 5' UTR End 5' UTR Start Ensembl Exon ID
             * Exon Chr Start (bp) Exon Chr End (bp) Exon Rank in Transcript phase Associated Gene DB
             */
			if (elems.length > 1) {
				String geneId = new String(elems[colgeneId]).intern();
				String transcriptId = new String(elems[coltranscriptId]).intern();
//				String protein = new String(elems[elementToColId.get("Ensembl Protein ID")]).intern();
				Chromosome chromosome = Chromosome.parseChr(elems[colchromosome]);
				Integer genestart = Integer.parseInt(elems[colgenestart]);
				Integer genestop = Integer.parseInt(elems[colgenestop]);
				Integer transcriptstart = Integer.parseInt(elems[coltranscriptstart]);
				Integer transcriptstop = Integer.parseInt(elems[coltranscriptstop]);
				Strand strand = Strand.parseStr(elems[colstrand]);
				if (strand.equals(Strand.NA)) {
					System.out.println("Could not parse strand str: " + elems[colstrand]);
				}
				String geneSymbol = new String(elems[colgeneSymbol]).intern();
				String exon = elems[colexon];
				Integer exonstart = Integer.parseInt(elems[colexonstart]);
				Integer exonend = Integer.parseInt(elems[colexonend]);
				Integer exonrankintranscript = Integer.parseInt(elems[colexonrankintranscript]);
				// String phase                   = elems[24];
//                    Integer cDNACodingStart         = Integer.parseInt(elems[25]);
//                    Integer genomicCodingStart      = Integer.parseInt(elems[26]);
//                    Integer cDNACodingEnd           = Integer.parseInt(elems[27]);
//                    Integer genomicCodingEnd        = Integer.parseInt(elems[28]);
//                    Integer CDSStart                = Integer.parseInt(elems[29]);
//                    Integer CDSEnd                  = Integer.parseInt(elems[30]);


				Gene currGen = null;
				Transcript currTra = null;
				Exon currExo = null;

				if (geneId.trim().length() > 0) {
					if (geneHash.get(geneId) == null) {
						Gene tmpGen = new Gene(geneId, chromosome, strand);
						tmpGen.setName(geneId);
						tmpGen.setGeneSymbol(geneSymbol);
						tmpGen.setStart(genestart);
						tmpGen.setStop(genestop);
						tmpGen.setStrand(strand);
						geneHash.put(geneId, tmpGen);
					}
					currGen = geneHash.get(geneId);
				}

				if (transcriptId.trim().length() > 0) {
					if (transcriptHash.get(transcriptId) == null) {
						Transcript tmpTra = new Transcript(transcriptId, chromosome, strand, currGen);
						tmpTra.setStart(transcriptstart);
						tmpTra.setStop(transcriptstop);
						transcriptHash.put(transcriptId, tmpTra);
					}

					currTra = transcriptHash.get(transcriptId);

					currGen.addTranscript(currTra);
				}

				if (exon.trim().length() > 0) {
					if (exonHash.get(exon) == null) {
						Exon tmpExo = new Exon(exon, chromosome, strand, currGen, exonstart, exonend);
						exonHash.put(exon, tmpExo);
					}
					currExo = exonHash.get(exon);
					currExo.addTranscript(currTra);
				}

				if (currTra != null && currExo != null) {
					currTra.addExon(currExo);
					currTra.setExonRank(currExo, exonrankintranscript);
				}

				if (colUTR3start != null && colUTR3stop != null) {
					String utrStaStr = elems[colUTR3start];
					String utrStoStr = elems[colUTR3stop];

					if (utrStaStr.length() > 0 && utrStoStr.length() > 0) {
						Integer sta = Integer.parseInt(utrStaStr);
						Integer sto = Integer.parseInt(utrStoStr);
						UTR utr = new UTR(chromosome, sta, sto);
						utr.setType(UTR.TYPE.PRIME3);
						utr.setTranscript(currTra);
						currTra.addUTR(utr);
					}

				}

				if (colUTR5start != null && colUTR5stop != null) {
					String utrStaStr = elems[colUTR5start];
					String utrStoStr = elems[colUTR5stop];
					if (utrStaStr.length() > 0 && utrStoStr.length() > 0) {
						Integer sta = Integer.parseInt(elems[colUTR5start]);
						Integer sto = Integer.parseInt(elems[colUTR5stop]);
						UTR utr = new UTR(chromosome, sta, sto);
						utr.setType(UTR.TYPE.PRIME5);
						utr.setTranscript(currTra);
						currTra.addUTR(utr);
					}
				}

			}

			if (lnCounter % 100000 == 0) {
				System.out.print(".");
			}
			elems = in.readLineElems(TextFile.tab);
			lnCounter++;
		}
		System.out.println("\tDone.");

		in.close();
		System.out.println("Loaded " + geneHash.size() + " genes, " + transcriptHash.size() + " transcripts, " + exonHash.size() + " exons.");
		strToGene = geneHash;
		genes = strToGene.values();
		genes.forEach(Gene::getBounds);


	}

}
