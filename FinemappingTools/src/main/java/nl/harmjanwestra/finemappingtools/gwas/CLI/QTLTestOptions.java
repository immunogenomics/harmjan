package nl.harmjanwestra.finemappingtools.gwas.CLI;

import org.apache.commons.cli.*;

/**
 * Created by hwestra on 3/15/16.
 */
public class QTLTestOptions {
	private static final Options OPTIONS;

	static {
		OPTIONS = new Options();
		Option option = Option.builder().longOpt("qtl").build();
		OPTIONS.addOption(option);


		option = Option.builder("o")
				.longOpt("out")
				.desc("Output path name")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("gtf")
				.desc("GTF Annotation path")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("e")
				.longOpt("exp")
				.desc("Expression matrix")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("c")
				.longOpt("cov")
				.desc("Covariate matrix")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("g")
				.longOpt("gen")
				.desc("Genotype VCF")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("gte")
				.desc("Coupling between genotype and expression samples")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("glimit")
				.desc("Limit to these genes")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("vlimit")
				.desc("Limit to these variants")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("climit")
				.desc("Limit to these covariates")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("cis")
				.desc("Cis distance window [default 1mb]")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("maf")
				.desc("Minor allele frequency threshold [default 0.01]")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("callrate")
				.desc("Callrate threshold [default 0.95]")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("impqual")
				.desc("Imputation quality (AR2) threshold [default: 0.3]")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("perm")
				.desc("Number of permutations per gene [default: 1000]")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("threads")
				.desc("Number of threads to use [default: 1]")
				.hasArg()
				.build();
		OPTIONS.addOption(option);
	}

	public String out;
	public String annotation;
	public String expression;
	public String covariates;
	public String genotype;
	public String genotypeToExpression;
	public String genelimit;
	public String variantlimit;
	public String covariatelimit;
	public int ciswindow = 1000000;
	public double mafthreshold = 0.01;
	public double callratethreshold = 0.95;
	public int nrpermutationspergene = 1000;
	public double impqualthreshold = 0.3;
	public int nrThreads = 1;
	public int distsize = 1000000;


	public QTLTestOptions(String[] args) {

		boolean run = true;
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, false);

			if (cmd.hasOption("o")) {
				this.out = cmd.getOptionValue("o");
			} else {
				System.out.println("Specify output");
				run = false;
			}

			if (cmd.hasOption("gtf")) {
				this.annotation = cmd.getOptionValue("gtf");
			} else {
				System.out.println("Specify annotation");
				run = false;
			}

			if (cmd.hasOption("exp")) {
				this.expression = cmd.getOptionValue("exp");
			} else {
				System.out.println("Specify expression");
				run = false;
			}

			if (cmd.hasOption("cov")) {
				this.covariates = cmd.getOptionValue("cov");
			}


			if (cmd.hasOption("gen")) {
				this.genotype = cmd.getOptionValue("gen");
			} else {
				System.out.println("Specify genotyping vcf");
				run = false;
			}

			if (cmd.hasOption("gte")) {
				this.genotypeToExpression = cmd.getOptionValue("gte");
			}

			if (cmd.hasOption("glimit")) {
				this.genelimit = cmd.getOptionValue("glimit");
			}

			if (cmd.hasOption("vlimit")) {
				this.variantlimit = cmd.getOptionValue("vlimit");
			}

			if (cmd.hasOption("climit")) {
				this.covariatelimit = cmd.getOptionValue("climit");
			}

			if (cmd.hasOption("cis")) {
				this.ciswindow = Integer.parseInt(cmd.getOptionValue("cis"));
			}

			if (cmd.hasOption("callrate")) {
				this.callratethreshold = Double.parseDouble(cmd.getOptionValue("callrate"));
			}


			if (cmd.hasOption("impqual")) {
				this.impqualthreshold = Double.parseDouble(cmd.getOptionValue("impqual"));
			}

			if (cmd.hasOption("perm")) {
				this.nrpermutationspergene = Integer.parseInt(cmd.getOptionValue("perm"));
			}

			if (cmd.hasOption("threads")) {
				this.nrThreads = Integer.parseInt(cmd.getOptionValue("threads"));
			}


			if (!run) {
				printHelp();
				System.exit(-1);
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
