package nl.harmjanwestra.proxyfinder;

import org.apache.commons.cli.*;

/**
 * Created by hwestra on 3/22/16.
 */
public class ProxyFinderOptions {
	private static final Options OPTIONS;

	static {
		OPTIONS = new Options();
		Option option = Option.builder().longOpt("proxy").build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("tabix")
				.hasArg()
				.desc("Prefix for tabix path [format /path/to/chrCHR.vcf.gz]. Replace the chromosome number with CHR (will be replaced by chr number depending on input SNP or region).")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("samplefilter")
				.hasArg()
				.desc("Limit samples to individuals in this list (one sample per line)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("w")
				.longOpt("windowsize")
				.hasArg()
				.desc("Window size [default 1000000]")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("m")
				.longOpt("maf")
				.hasArg()
				.desc("MAF threshold [default: 0.005]")
				.build();
		OPTIONS.addOption(option);


		option = Option.builder("h")
				.longOpt("hwep")
				.hasArg()
				.desc("HWE-P threshold [default: 0.0001]")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("t")
				.longOpt("threshold")
				.hasArg()
				.desc("R-squared threshold [default: 0.8]")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("i")
				.longOpt("snps")
				.hasArg()
				.desc("SNP path (format: 3 or 6 columns, tab separated, one or two snps per line: chr pos rsid) ")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("r")
				.longOpt("regions")
				.hasArg()
				.desc("Region bed file path")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("pairwise")
				.desc("Perform Pairwise LD calculation")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("vcf")
				.hasArg()
				.desc("Use non-indexed VCF as input")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("locusld")
				.desc("Perform Pairwise LD calculation within a region (provide region with --regions)")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("o")
				.longOpt("out")
				.hasArg()
				.desc("Output path")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("threads")
				.hasArg()
				.desc("Nr of threads [default: 1]")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("matchrsid")
				.desc("Match variants on RS id")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("dontaddself")
				.desc("Don't add input as proxy to itself.")
				.build();
		OPTIONS.addOption(option);

	}

	enum MODE {
		PAIRWISE, LOCUSLD, CMWINDOW, PROXY
	}

	public MODE mode;
	public boolean matchrsid = false;
	public String tabixrefprefix;
	public String samplefilter;
	public int windowsize = 1000000;
	public double threshold = 0.8;
	public String snpfile;
	public String output;
	public int nrthreads = 1;
	public String regionfile;
	public double mafthreshold = 0.005;
	public String vcf;
	public double hwepthreshold = 0.0001;
	public boolean addSNPasProxyToItself = true;

	public ProxyFinderOptions(String[] args) {

		boolean run = true;
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, false);

			if (cmd.hasOption("tabix")) {
				tabixrefprefix = cmd.getOptionValue("tabix");
			} else if (cmd.hasOption("vcf")) {
				vcf = cmd.getOptionValue("vcf");
			} else {

				System.out.println("Provide reference");
				run = false;
			}

			if (cmd.hasOption("matchrsid")) {
				matchrsid = true;
			}

			if (cmd.hasOption("w")) {
				windowsize = Integer.parseInt(cmd.getOptionValue("w"));
			}

			if (cmd.hasOption("samplefilter")) {
				samplefilter = cmd.getOptionValue("samplefilter");
			}

			if (cmd.hasOption("proxy")) {
				mode = MODE.PROXY;
			} else if (cmd.hasOption("pairwise")) {
				mode = MODE.PAIRWISE;
			} else if (cmd.hasOption("locusld")) {
				mode = MODE.LOCUSLD;
			} else {
				System.out.println("Please specify mode.");
				run = false;
			}

			if (cmd.hasOption("t")) {
				threshold = Double.parseDouble(cmd.getOptionValue("t"));
			}

			if (cmd.hasOption("i")) {
				snpfile = cmd.getOptionValue("i");
			}

			if (cmd.hasOption("r")) {
				regionfile = cmd.getOptionValue("r");
			}

			if (cmd.hasOption("o")) {
				output = cmd.getOptionValue("o");
			} else {
				System.out.println("Provide output");
				run = false;
			}

			if (cmd.hasOption("maf")) {
				mafthreshold = Double.parseDouble(cmd.getOptionValue("maf"));
			}

			if (cmd.hasOption("hwep")) {
				hwepthreshold = Double.parseDouble(cmd.getOptionValue("hwep"));
			}

			if (cmd.hasOption("threads")) {
				nrthreads = Integer.parseInt(cmd.getOptionValue("threads"));
			}


		} catch (ParseException e) {
			e.printStackTrace();


		}

		if (!run) {
			printHelp();
			System.exit(-1);
		} else {
			printOptions();
		}
	}

	private void printOptions() {
		System.out.println("Parameters:");
		System.out.println("MODE:\t" + mode);
		System.out.println("matchrsid:\t" + matchrsid);
		System.out.println("tabixrefprefix:\t" + tabixrefprefix);
		System.out.println("samplefilter:\t" + samplefilter);
		System.out.println("windowsize:\t" + windowsize);
		System.out.println("threshold:\t" + threshold);
		System.out.println("snpfile:\t" + snpfile);
		System.out.println("output:\t\t" + output);
		System.out.println("nrthreads:\t" + nrthreads);
		System.out.println("regionfile:\t" + regionfile);
		System.out.println("mafthreshold:\t" + mafthreshold);
		System.out.println("hwepthreshold:\t" + hwepthreshold);
		System.out.println("vcf:\t\t" + vcf);
		System.out.println();

	}

	public void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
		System.exit(-1);
	}

}
