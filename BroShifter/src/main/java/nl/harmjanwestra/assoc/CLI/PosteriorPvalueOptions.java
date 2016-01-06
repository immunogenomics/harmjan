package nl.harmjanwestra.assoc.CLI;

import org.apache.commons.cli.*;

/**
 * Created by hwestra on 11/23/15.
 */
public class PosteriorPvalueOptions {
	private String assocFile;
	private String regionFile;
	private String outputPrefix;
	private double bayesThreshold;

	public String getAssocFile() {
		return assocFile;
	}

	public String getRegionFile() {
		return regionFile;
	}

	public String getOutputPrefix() {
		return outputPrefix;
	}

	public double getBayesThreshold() {
		return bayesThreshold;
	}

	private static final Options OPTIONS;

	static {
		OPTIONS = new Options();
		Option option = Option.builder().longOpt("posterior").build();
		OPTIONS.addOption(option);


		option = Option.builder("i")
				.hasArg()
				.desc("Input VCF")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("o")
				.hasArg()
				.desc("Output prefix")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("b")
				.hasArg()
				.desc("Bayesian Threshold")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("r")
				.hasArg()
				.desc("Region file")
				.build();
		OPTIONS.addOption(option);


	}

	public PosteriorPvalueOptions(String[] args) {


		boolean run = true;
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, false);

			if (cmd.hasOption("i")) {
				assocFile = cmd.getOptionValue("i");
			} else {
				System.out.println("Please provide input: -i");
				run = false;
			}

			if (cmd.hasOption("o")) {
				outputPrefix = cmd.getOptionValue("o");
			} else {
				System.out.println("Please provide output: -o");
				run = false;
			}

			if (cmd.hasOption("o")) {
				regionFile = cmd.getOptionValue("r");
			} else {
				System.out.println("Please provide regionfile: -r");
				run = false;
			}

			if (cmd.hasOption("b")) {
				String nstr = cmd.getOptionValue("b");
				try {
					bayesThreshold = Double.parseDouble(nstr);
				} catch (NumberFormatException e) {
					System.out.println(nstr + " is not a decimal/double for -b");
					run = false;
				}
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
