package nl.harmjanwestra.finemappingtools.gwas.CLI;

import org.apache.commons.cli.*;

/**
 * Created by Harm-Jan on 02/28/16.
 */
public class CaviarOptions {
	private static final Options OPTIONS;

	static {
		OPTIONS = new Options();
		Option option = Option.builder().longOpt("caviar").build();
		OPTIONS.addOption(option);

		option = Option.builder()

				.longOpt("convert")
				.desc("Convert association path")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()

				.longOpt("filter")
				.desc("Convert association path")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("i")
				.hasArg()
				.desc("Association path")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("m")
				.hasArg()
				.desc("Matrix")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("o")
				.hasArg()
				.desc("Output path")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("l")
				.hasArg()
				.desc("Variant list")
				.build();
		OPTIONS.addOption(option);
	}

	public String assoc;
	public String variantlist;
	public String matrix;
	public String out;

	public boolean filter;
	public boolean convert;

	public CaviarOptions(String[] args) {

		boolean run = true;
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, false);

			if (cmd.hasOption("i")) {
				assoc = cmd.getOptionValue("i");
			}
			if (cmd.hasOption("m")) {
				matrix = cmd.getOptionValue("m");
			}
			if (cmd.hasOption("o")) {
				out = cmd.getOptionValue("o");
			}

			if (cmd.hasOption("l")) {
				variantlist = cmd.getOptionValue("l");
			}

			if (cmd.hasOption("filter")) {
				filter = true;
			} else if (cmd.hasOption("convert")) {
				convert = true;
			} else {
				run = false;
			}
		} catch (ParseException e) {
			e.printStackTrace();

		}

		if (!run) {
			printHelp();
			System.exit(-1);
		}
	}

	public void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
		System.exit(-1);
	}
}
