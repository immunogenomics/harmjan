package nl.harmjanwestra.broshifter.CLI;

import org.apache.commons.cli.*;

/**
 * Created by hwestra on 11/23/15.
 */
public class MainOptions {


	public MODE mode = MODE.NA;

	public enum MODE {
		BROSHIFTER,
		POSTERIORPVAL,
		PLOT,
		MERGE,
		ASSOC, NA
	}

	public static final Options OPTIONS;

	static {
		OPTIONS = new Options();

		Option option;
		option = Option.builder()
				.desc("Run broshifter - enrichment analysis with posteriors")
				.longOpt("broshifter")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Merge association files")
				.longOpt("merge")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Plot association results")
				.longOpt("plot")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Calculate posteriors")
				.longOpt("posterior")
				.build();
		OPTIONS.addOption(option);

	}


	public MainOptions(String[] args) {
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, true);

			if (cmd.hasOption("broshifter")) {
				System.out.println("About to run broshifter");
				mode = MODE.BROSHIFTER;
			} else if (cmd.hasOption("merge")) {
				mode = MODE.MERGE;
			} else if (cmd.hasOption("plot")) {
				mode = MODE.PLOT;
			} else if (cmd.hasOption("posterior")) {
				mode = MODE.POSTERIORPVAL;
			} else {
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp(" ", OPTIONS);
			}
		} catch (ParseException ex) {
			System.err.println();
			mode = MODE.NA;
			printHelp();
		}


	}

	public void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
		System.exit(-1);
	}

}
