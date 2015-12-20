package nl.harmjanwestra.ngs;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import nl.harmjanwestra.utilities.bamfile.BamFileReader;
import nl.harmjanwestra.utilities.features.Chromosome;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashMap;
import java.util.Set;

/**
 * Created by hwestra on 9/10/15.
 */
public class CountNumberOfReadsPerPair {

	public static void main(String[] args) {
		try {
			if (args.length < 2) {
				System.out.println("Usage: bamfile outfile");
			} else {
				CountNumberOfReadsPerPair c = new CountNumberOfReadsPerPair();
				c.run(args[0], args[1]);
			}
		} catch (IOException e) {

		}
	}

	public void run(String bamin, String fileout) throws IOException {
		HashMap<String, Integer> alignmentsPerRead = new HashMap<String, Integer>();

		BamFileReader reader = new BamFileReader(bamin);
		SAMRecordIterator it = reader.iterator();
		int ctr = 0;


		while (it.hasNext() && ctr < 10000000) {
			SAMRecord record = it.next();

			String reg = record.getReferenceName();
			Chromosome chr = Chromosome.parseChr(reg);

			String name = record.getReadName();
			if (!chr.equals(Chromosome.MT) && !chr.equals(Chromosome.NA)) {
				Integer i = alignmentsPerRead.get(name);
				if (i == null) {
					i = 1;
				} else {
					i++;
				}

				ctr++;
				if (ctr % 10000000 == 0) {
					System.out.println(ctr + " alignments processed.");
				}

				alignmentsPerRead.put(name, i);

			}

		}

		Set<String> keys = alignmentsPerRead.keySet();
		TextFile out = new TextFile(fileout, TextFile.W);
		int nrWMoreThan2Alignments = 0;
		int nrWithSingleAlignment = 0;
		int normal = 0;
		for (String key : keys) {
			int i = alignmentsPerRead.get(key);
			if (i == 1) {
				nrWithSingleAlignment++;

			} else if (i == 2) {
				normal++;
			} else {

				nrWMoreThan2Alignments++;
			}

			out.writeln(key + "\t" + i);
		}
		System.out.println(nrWithSingleAlignment + " ==1 \t" + normal + " == 2\t" + nrWMoreThan2Alignments + " ==>2");
		out.close();
	}

}
