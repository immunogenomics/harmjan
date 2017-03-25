package nl.harmjanwestra.finemappingtools.gwas.CLI;

import org.apache.commons.cli.*;

/**
 * Created by hwestra on 1/5/16.
 */
public class AssociationResultMergerOptions {

	private static final Options OPTIONS;

	static {
		OPTIONS = new Options();
		Option option = Option.builder().longOpt("merge").build();
		OPTIONS.addOption(option);

		option = Option.builder("i")
				.hasArg()
				.desc("Input VCF (or comma separated list when not using --concat)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("o")
				.hasArg()
				.desc("Output prefix")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("concat")
				.desc("Output prefix")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("c")
				.desc("Recalculate posteriors")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("n")
				.hasArg()
				.desc("Names of datasets (comma separated)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("b")
				.hasArg()
				.desc("Bayes threshold for credible sets (default = 0.95)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("r")
				.hasArg()
				.desc("Region bed path")
				.build();
		OPTIONS.addOption(option);

	}

	private boolean concat = false;
	private double bayesthreshold = 0.95;
	private String input;
	private String outputprefix;
	private String refStr;
	private String regions;
	private boolean recalculatePosteriors = false;

	public AssociationResultMergerOptions(String[] args) {
		boolean run = true;
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, false);

			if (cmd.hasOption("concat")) {
				concat = true;
			}
			if (cmd.hasOption("c")) {
				recalculatePosteriors = true;
			}

			if (cmd.hasOption("i")) {
				input = cmd.getOptionValue("i");
			} else {
				System.out.println("Please provide input with -i");
				run = false;
			}

			if (cmd.hasOption("o")) {
				outputprefix = cmd.getOptionValue("o");
			} else {
				System.out.println("Please provide output prefix with -o");
				run = false;
			}

			if (cmd.hasOption("n")) {
				refStr = cmd.getOptionValue("n");
			} else if (!concat) {
				System.out.println("Please provide datasetnames with -n");
				run = false;
			}

			if (cmd.hasOption("b")) {
				String nstr = cmd.getOptionValue("b");
				try {
					bayesthreshold = Double.parseDouble(nstr);
				} catch (NumberFormatException e) {
					System.out.println(nstr + " is not a decimal/double for -b");
				}
			}

			if (cmd.hasOption("r")) {
				regions = cmd.getOptionValue("r");
				System.out.println("Please provide regions with -r");
			} else if (!concat) {
				run = false;
			}

			if (!run) {
				printHelp();
				System.exit(-1);
			}

		} catch (ParseException e) {
			e.printStackTrace();
		}

	}

	public boolean isConcat() {
		return concat;
	}

	public double getBayesthreshold() {
		return bayesthreshold;
	}

	public String getInput() {
		return input;
	}

	public String getOutputprefix() {
		return outputprefix;
	}

	public String getRefStr() {
		return refStr;
	}

	public String getRegions() {
		return regions;
	}

	public boolean isRecalculatePosteriors() {
		return recalculatePosteriors;
	}

	public void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
		System.exit(-1);
	}
}