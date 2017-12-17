package nl.harmjanwestra.finemappingtools.gwas.CLI;

import org.apache.commons.cli.*;

public class AssociationFilterOptions {
	private static final Options OPTIONS;
	
	static {
		OPTIONS = new Options();
		Option option = Option.builder().longOpt("filterassoc").build();
		OPTIONS.addOption(option);
		
		option = Option.builder("i")
				.hasArg()
				.desc("Input assoc file")
				.build();
		OPTIONS.addOption(option);
		
		option = Option.builder("o")
				.hasArg()
				.desc("Output assoc file")
				.build();
		OPTIONS.addOption(option);
		
		option = Option.builder("s")
				.hasArg()
				.desc("List of SNP combinedIDs to filter out.")
				.build();
		OPTIONS.addOption(option);
		
		option = Option.builder("p")
				.desc("Association file is a pairwise file")
				.build();
		OPTIONS.addOption(option);
	}
	
	private boolean isPairwise;
	
	String input;
	String outputprefix;
	String snploc;
	
	public AssociationFilterOptions(String[] args) {
		boolean run = true;
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, false);
			
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
			
			
			if (cmd.hasOption("s")) {
				snploc = cmd.getOptionValue("s");
			} else {
				System.out.println("Please provide snp list with -s");
				run = false;
			}
			
			if (cmd.hasOption("p")) {
				isPairwise = true;
			}
			
			
			if (!run) {
				printHelp();
				System.exit(-1);
			}
			
		} catch (ParseException e) {
			e.printStackTrace();
		}
		
	}
	
	
	public void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
		System.exit(-1);
	}
	
	
	public boolean isPairwise() {
		return isPairwise;
	}
	
	public String getInput() {
		return input;
	}
	
	public String getOutputprefix() {
		return outputprefix;
	}
	
	public String getSnploc() {
		return snploc;
	}
}
