package nl.harmjanwestra.txtr;

import org.apache.commons.cli.*;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;

/**
 * Created by Harm-Jan on 02/01/16.
 */
public class TXTr {

	private static Options OPTIONS;

	static {
		OPTIONS = new Options();

		Option option;

		option = Option.builder()
				.desc("Split textfile by lines")
				.longOpt("split")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Multiple line hashtag delimited header")
				.longOpt("multilineheader")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Merge while skipping first line of all files except 1st")
				.longOpt("merge")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("i")
				.desc("Input")
				.hasArg().required()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("n")
				.desc("Nr lines for splitting")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("o")
				.desc("Input")
				.hasArg().required()
				.build();
		OPTIONS.addOption(option);
	}

	public static void main(String[] args) {
		TXTr t = new TXTr();
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, true);

			String input = "";
			String output = "";
			if (cmd.hasOption("i")) {
				input = cmd.getOptionValue("i");
			}
			if (cmd.hasOption("o")) {
				output = cmd.getOptionValue("o");
			}

			if (cmd.hasOption("split")) {

				if (cmd.hasOption("n")) {
					try {
						int nrLines = Integer.parseInt(cmd.getOptionValue("n"));
						try {
							t.split(input, output, nrLines);
						} catch (IOException e) {
							e.printStackTrace();
						}
					} catch (NumberFormatException e) {
						System.out.println(cmd.getOptionValue("n") + " is not an integer");
					}
				} else {
					System.out.println("Use -n for --split");

					printHelp();
				}
			} else if (cmd.hasOption("merge")) {
				try {
					boolean multilineheader = false;
					if (cmd.hasOption("multilineheader")) {
						multilineheader = true;
					}
					t.mergeSkipHeader(input, output, multilineheader);
				} catch (IOException e) {
					e.printStackTrace();
				}

			}


		} catch (ParseException e) {
			printHelp();
			e.printStackTrace();
		}

	}

	public static void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
		System.exit(-1);
	}

	public void split(String file, String fileout, int lns) throws IOException {

		TextFile in = new TextFile(file, TextFile.R);
		int ctr = 0;
		int ctr2 = 1;
		String ln = in.readLine();
		TextFile out = new TextFile(fileout + "-" + ctr2, TextFile.W);
		while (ln != null) {
			ln = in.readLine();
			out.writeln(ln);
			ctr++;
			if (ctr % lns == 0) {
				out.close();
				out = new TextFile(fileout + "-" + ctr2, TextFile.W);
				ctr2++;
			}
		}
		out.close();

		in.close();
	}

	public void mergeSkipHeader(String fileList, String output, boolean multilinehashtagheader) throws IOException {

		TextFile tf1 = new TextFile(fileList, TextFile.R);
		String[] files = tf1.readAsArray();
		tf1.close();

		boolean headerwritten = false;

		TextFile out = new TextFile(output, TextFile.W);
		int fctr = 0;
		for (String file : files) {
			TextFile in = new TextFile(file, TextFile.R);
			String ln = in.readLine();
			int lnctr = 0;
			while (ln != null) {
				if (multilinehashtagheader) {
					if (ln.startsWith("#")) {
						if (fctr == 0) {
							out.writeln(ln);
						}
					} else {
						out.writeln(ln);
					}
				} else {
					if (!headerwritten && lnctr == 0) {
						out.writeln(ln);
						headerwritten = true;
					} else if (lnctr > 0) {
						out.writeln(ln);
					}
				}

				ln = in.readLine();
				lnctr++;
			}

			in.close();
			fctr++;
		}
		out.close();

	}

}
