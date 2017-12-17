package nl.harmjanwestra.finemappingtools.gwas;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import org.apache.commons.cli.*;

import java.io.IOException;

public class PLINKPCAConvert {
	
	public PLINKPCAConvert(String[] args) throws IOException {
		
		CommandLineParser parser = new DefaultParser();
		
		Options options = new Options();
		
		Option plopt = Option.builder("p").desc("Meh.").longOpt("plink").build();
		Option inopt = Option.builder("i").argName("path").hasArg().desc("Location and name of the input file.").longOpt("in").build();
		Option outopt = Option.builder("o").argName("path").hasArg().desc("Location and name of the output file.").longOpt("out").build();
		Option nrpcasopt = Option.builder("n").argName("path").hasArg().desc("Nr of PCAs to convert").longOpt("ncpcas").build();
		
		options.addOption(plopt).addOption(inopt).addOption(outopt).addOption(nrpcasopt);
		
		
		String in = "";
		String out = "";
		int nrpcas = 0;
		
		
		try {
			CommandLine cmd = parser.parse(options, args, false);
			if (cmd.hasOption("i")) {
				in = cmd.getOptionValue("i");
			}
			if (cmd.hasOption("o")) {
				out = cmd.getOptionValue("o");
			}
			if (cmd.hasOption("n")) {
				nrpcas = Integer.parseInt(cmd.getOptionValue("n"));
			}
		} catch (ParseException e) {
			e.printStackTrace();
		}
		
		
		run(in, out, nrpcas);
	}
	
	public void run(String in, String out, int nrpcaz) throws IOException {
		TextFile tf = new TextFile(in, TextFile.R);
		
		String header = "-";
		String[] elems = tf.readLineElems(Strings.whitespace);
		
		int stop = elems.length;
		if (nrpcaz > 0) {
			stop = 2 + nrpcaz;
		}
		
		for (int i = 2; i < stop; i++) {
			header += "\tPC" + (i - 1);
		}
		
		TextFile tfw = new TextFile(out, TextFile.W);
		tfw.writeln(header);
		while (elems != null) {
			String ln = elems[1];
			for (int i = 2; i < stop; i++) {
				ln += "\t" + elems[i];
			}
			tfw.writeln(ln);
			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();
		tfw.close();
	}
}
