package nl.harmjanwestra.gwas.CLI;

import org.apache.commons.cli.*;

/**
 * Created by Harm-Jan on 04/20/16.
 */
public class CountVariantsOptions {

	private static final Options OPTIONS;

	static {
		OPTIONS = new Options();
		Option option = Option.builder().longOpt("gwas").build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("countvariants")
				.desc("Count variants")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("i")
				.hasArg()
				.desc("A list of log files to parse")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("r")
				.hasArg()
				.desc("A list of regions of interest")
				.build();
		OPTIONS.addOption(option);
	}

	public String bedfile;
	public String input;

	public CountVariantsOptions(String[] args) {

		boolean run = true;
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, false);

			if (cmd.hasOption("r")) {
				bedfile = cmd.getOptionValue("r");
			} else {
				System.out.println("Please provide input: -r");
				run = false;
			}

			if (cmd.hasOption("i")) {
				input = cmd.getOptionValue("i");
			} else {
				System.out.println("Please provide input: -i");
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
