package nl.harmjanwestra.coverage.CLI;

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
				.desc("Convert bam to bed file")
				.longOpt("bamtobed")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Correlate Coverage Within Regions")
				.longOpt("correlatebamregions")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Coverage per chromosome")
				.longOpt("coverageperchr")
				.build();
		OPTIONS.addOption(option);




		/*
		c.bamToBed(options.bamfile, options.outdir);
			c.bamToBedWithinRegionsForList(options.list, options.outdir, options.targetregions, options.outputcoverageperregion, options.threads);
			c.correlateCoverageWithinRegionFiles(options.targetregions, options.inputdir, options.outdir);
			c.coveragePerChromosome(options.bamfile, options.outdir, options.posshift, options.negshift);
		 */

		option = Option.builder()
				.desc("Bam file")
				.longOpt("bam")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Outdir")
				.longOpt("out")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("List file")
				.longOpt("list")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Target regions")
				.longOpt("regions")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Output a file per region")
				.longOpt("outputregions")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Threads")
				.longOpt("threads")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Positive shift")
				.longOpt("shiftpos")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.desc("Negative shift")
				.longOpt("shiftneg")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

	}

	public String bamfile;
	public String outdir;
	public String list;
	public String targetregions;
	public boolean outputcoverageperregion;
	public int threads = 1;
	public String inputdir;
	public int posshift = 0;
	public int negshift = 0;

	public MODE mode = MODE.NA;

	public MainOptions(String[] args) {
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, true);

			//				.longOpt("bam")
			if (cmd.hasOption("bam")) {
				bamfile = cmd.getOptionValue("bam");
			}
//				.longOpt("out")
			if (cmd.hasOption("out")) {
				outdir = cmd.getOptionValue("out");
			}
//				.longOpt("list")
			if (cmd.hasOption("list")) {
				list = cmd.getOptionValue("list");
			}
//				.longOpt("regions")
			if (cmd.hasOption("regions")) {
				targetregions = cmd.getOptionValue("regions");
			}
//				.longOpt("outputregions")
			if (cmd.hasOption("outputregions")) {
				outputcoverageperregion = true;
			}
//				.longOpt("threads")
			if (cmd.hasOption("threads")) {
				threads = Integer.parseInt(cmd.getOptionValue("threads"));
			}
//				.longOpt("shiftpos")
			if (cmd.hasOption("shiftpos")) {
				posshift = Integer.parseInt(cmd.getOptionValue("shiftpos"));

			}
//				.longOpt("shiftneg")
			if (cmd.hasOption("shiftneg")) {
				negshift = Integer.parseInt(cmd.getOptionValue("shiftneg"));
			}



			boolean run = true;
			if (cmd.hasOption("bamtobed")) {
				if (targetregions != null) {
					mode = MODE.BAMTOBEDREGIONS;
					if (list == null || outdir == null) {
						System.out.println("Please provide bam file list and outdir");
						run = false;
					}
				} else {
					mode = MODE.BAMTOBED;
					if (bamfile == null || outdir == null) {
						System.out.println("Please provide bam file and outdir");
						run = false;
					}
				}
			} else if (cmd.hasOption("correlatebamregions")) {
				mode = MODE.CORRELATEBAMS;
				if (targetregions == null || outdir == null || inputdir == null) {
					System.out.println("Please provide target regions, outdir and inputdir");
					run = false;
				}
			} else if (cmd.hasOption("coverageperchr")) {
				mode = MODE.COVERAGEPERCHR;
				if (bamfile == null || outdir == null ) {
					System.out.println("Please provide outdir and inputfile");
					run = false;
				}
			} else {
				run = false;
			}

			if (!run) {
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp(" ", OPTIONS);
			}
		} catch (ParseException ex) {
			System.err.println();
			printHelp();
		}


	}

	public void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
		System.exit(-1);
	}

	public enum MODE {
		NA, BAMTOBEDREGIONS, BAMTOBED, CORRELATEBAMS, COVERAGEPERCHR
	}


}
