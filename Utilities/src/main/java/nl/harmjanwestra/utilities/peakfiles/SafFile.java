package nl.harmjanwestra.utilities.peakfiles;

import nl.harmjanwestra.utilities.features.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.Strand;
import umcg.genetica.io.text.TextFile;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by hwestra on 10/21/15.
 */
public class SafFile extends TextFile {

	public SafFile(String file, boolean mode) throws IOException {
		super(file, mode);
	}

	public SafFile(File file, boolean mode) throws IOException {
		super(file, mode);
	}

	public SafFile(File file, boolean mode, int buffersize) throws IOException {
		super(file, mode, buffersize);
	}

	public SafFile(String file, boolean mode, int buffersize) throws IOException {
		super(file, mode, buffersize);
	}


	public ArrayList<Feature> read() throws IOException {
		close();
		open();

		String[] elems = readLineElems(TextFile.tab);
		ArrayList<Feature> output = new ArrayList<Feature>();
		while (elems != null) {

			if (elems.length >= 5) {
				String name = elems[0];
				String chr = elems[1];
				String start = elems[2];
				String stop = elems[3];
				String strand = elems[4];

				Strand s = Strand.parseStr(strand);
				Integer sta = Integer.parseInt(start);
				Integer sto = Integer.parseInt(stop);

				Chromosome c = Chromosome.parseChr(chr);
				Feature f = new Feature();
				f.setStart(sta);
				f.setStop(sto);
				f.setChromosome(c);
				f.setStrand(s);
				f.setName(name);
				output.add(f);
			}
			elems = readLineElems(TextFile.tab);


		}
		return output;
	}
}
