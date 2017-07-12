package nl.harmjanwestra.finemappingtools.broshifter.CLI;

import org.apache.commons.cli.*;

/**
 * Created by hwestra on 11/23/15.
 */
public class MainOptions {


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
				.desc("Run goshifter - enrichment analysis")
				.longOpt("goshifter")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Run guess converter")
				.longOpt("guess")
				.build();
		OPTIONS.addOption(option);


		option = Option.builder()
				.desc("Merge association files")
				.longOpt("merge")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("QTL analysis")
				.longOpt("qtl")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("CaviarBF helper functions")
				.longOpt("caviar")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Plot association results")
				.longOpt("plotposteriors")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Plot annotation overlap results")
				.longOpt("plotoverlap")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Filter bed region path for significant results")
				.longOpt("bedfilter")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Calculate posteriors")
				.longOpt("posterior")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Perform association analysis")
				.longOpt("gwas")
				.build();
		OPTIONS.addOption(option);


		option = Option.builder()
				.desc("Find LD partners of variants given a VCF path")
				.longOpt("proxy")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Parse GWAS logs to determine number of tested variants")
				.longOpt("countvariants")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Update RS ids for association files")
				.longOpt("updaters")
				.build();
		OPTIONS.addOption(option);

	}

	public MODE mode = MODE.NA;

	public MainOptions(String[] args) {
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, true);

			if (cmd.hasOption("broshifter")) {
				mode = MODE.BROSHIFTER;
			} else if (cmd.hasOption("merge")) {
				mode = MODE.MERGE;
			} else if (cmd.hasOption("plotoverlap")) {
				mode = MODE.ANNOTATIONOVERLAPPLOT;
			} else if (cmd.hasOption("bedfilter")) {
				mode = MODE.BEDFILTER;
			} else if (cmd.hasOption("plotposteriors")) {
				mode = MODE.PLOTPOSTERIORS;
			} else if (cmd.hasOption("posterior")) {
				mode = MODE.POSTERIORPVAL;
			} else if (cmd.hasOption("goshifter")) {
				mode = MODE.GOSHIFTER;
			} else if (cmd.hasOption("gwas")) {
				mode = MODE.ASSOC;
			} else if (cmd.hasOption("caviar")) {
				mode = MODE.CAVIAR;
			} else if (cmd.hasOption("proxy")) {
				mode = MODE.PROXYFINDER;
			} else if (cmd.hasOption("qtl")) {
				mode = MODE.QTL;
			} else if (cmd.hasOption("guess")) {
				mode = MODE.GUESS;
			} else if (cmd.hasOption("countvariants")) {
				mode = MODE.COUNTVARIANTS;
			} else if (cmd.hasOption("updaters")) {
				mode = MODE.UPDATERS;
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

	public enum MODE {
		BROSHIFTER,
		POSTERIORPVAL,
		PLOTPOSTERIORS,
		MERGE,
		ASSOC, ANNOTATIONOVERLAPPLOT, BEDFILTER, CAVIAR, QTL, PROXYFINDER, GOSHIFTER, COUNTVARIANTS, UPDATERS, STATS, GUESS, NA
	}

}
