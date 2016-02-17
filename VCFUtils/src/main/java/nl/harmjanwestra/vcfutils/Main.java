package nl.harmjanwestra.vcfutils;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.oxford.HapSample;
import nl.harmjanwestra.utilities.vcf.VCFFunctions;
import org.apache.commons.cli.*;

import java.io.IOException;
import java.util.Random;

/**
 * Created by hwestra on 11/24/15.
 */
public class Main {

    private static Options OPTIONS;

    static {
        OPTIONS = new Options();

        Option option;

        option = Option.builder()
                .desc("Batch splitter")
                .longOpt("split")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Replace sample names")
                .longOpt("samplereplace")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Sample filter")
                .longOpt("filtersample")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Detect genotyping sample mix-ups")
                .longOpt("mixups")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Plot r-squared")
                .longOpt("plotrsq")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Plot correlations")
                .longOpt("plotcor")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Update rs names in a VCF")
                .longOpt("updaters")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("DB SNP VCF file")
                .longOpt("dbsnp")
                .hasArg()
                .build();
        OPTIONS.addOption(option);


        option = Option.builder()
                .desc("Filter variants on maf,callrate,genotypingquals and allelic depth")
                .longOpt("filter")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Filter for bed regions")
                .longOpt("bedfilter")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Reintroduce unimputed variants (can also be used to merge two VCF files)")
                .longOpt("reintroduce")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Merge two VCF files (with non-overlapping samples) as a big chunk of swiss cheese; with holes that is.")
                .longOpt("mergecheese")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Merge overlapping variants only. Do strand checks etc.")
                .longOpt("merge")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Correlate imputation output")
                .longOpt("correlate")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Generate Random seed")
                .longOpt("rng")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Dedup (very naive)")
                .longOpt("dedup")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Replace missing codes")
                .longOpt("rmc")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Convert Hap/sample")
                .longOpt("hapsample")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Variant sampler")
                .longOpt("sample")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Chr splitter")
                .longOpt("splitchr")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Merge imputation batches")
                .longOpt("mergeimputation")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Make pseudo controls")
                .longOpt("pseudo")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Combine correlation results")
                .longOpt("combine")
                .build();
        OPTIONS.addOption(option);


        option = Option.builder()
                .desc("RSquared threshold")
                .argName("double")
                .longOpt("rsquared")
                .hasArg()
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Callrate threshold")
                .argName("double")
                .longOpt("callrate")
                .hasArg()
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Genotype qual threshold")
                .argName("double")
                .longOpt("gqual")
                .hasArg()
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Read depth threshold")
                .argName("double")
                .longOpt("readdepth")
                .hasArg()
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Allelic balance")
                .argName("double")
                .longOpt("allelicbalance")
                .hasArg()
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Maf threshold")
                .argName("double")
                .longOpt("maf")
                .hasArg()
                .build();
        OPTIONS.addOption(option);

        option = Option.builder("i")
                .desc("Input VCF")
                .argName("file")
                .hasArg()
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Reference name")
                .argName("string")
                .longOpt("refname")
                .hasArg()
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Number of iterations")
                .argName("int")
                .longOpt("nriter")
                .hasArg()
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .desc("Number of batches")
                .argName("int")
                .longOpt("nrbatches")
                .hasArg()
                .build();
        OPTIONS.addOption(option);

        option = Option.builder("i2")
                .desc("Input VCF #2")
                .argName("file")
                .hasArg()
                .build();
        OPTIONS.addOption(option);

        option = Option.builder("f")
                .argName("file")
                .hasArg()
                .desc("Input family information (PED or FAM file)")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder("b")
                .hasArg()
                .desc("Region BED file")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder("l")
                .hasArg()
                .desc("List file")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder("r")
                .hasArg()
                .desc("Reference VCF")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder("o")
                .hasArg()
                .desc("Output prefix")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder("s")
                .hasArg()
                .desc("Random number generator seed")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder("p")
                .argName("double")
                .hasArg()
                .desc("Sampling percentage for sampling")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder("n")
                .hasArg()
                .argName("int")
                .desc("Batch size for splitting")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .longOpt("vcfsort")
                .hasArg()
                .argName("path")
                .desc("Location of vcfsort")
                .build();
        OPTIONS.addOption(option);

        option = Option.builder()
                .longOpt("linux")
                .desc("Is machine a linux machine (default: no)")
                .build();
        OPTIONS.addOption(option);


    }

    public static void main(String[] args) {

        try {
            CommandLineParser parser = new DefaultParser();
            final CommandLine cmd = parser.parse(OPTIONS, args, true);

            boolean run = true;
            String ref = null;
            String input = null;
            String out = null;
            double perc = 1;


            if (cmd.hasOption("rng")) {
                Random random = new Random();

                System.out.println("Here's your random number :)");
                long seed = random.nextLong();
                System.out.println(seed);
                System.exit(-1);
            }


            if (!cmd.hasOption("i")) {
                run = false;
                System.out.println("Please provide input vcf");
            } else {
                input = cmd.getOptionValue("i");
            }


            if (!cmd.hasOption("o")) {
                System.out.println("Please provide output prefix");
                run = false;
            } else {
                out = cmd.getOptionValue("o");

            }

            if (cmd.hasOption("mixups")) {

                MixupTest test = new MixupTest();
                try {
                    if (cmd.hasOption("i2")) {
                        test.test(input, cmd.getOptionValue("i2"), out);
                    } else {
                        System.out.println("Use -i2 for --mixups");
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }

            } else if (cmd.hasOption("samplereplace")) {
                try {
                    if (cmd.hasOption("l")) {
                        VCFSampleNameReplace vpl = new VCFSampleNameReplace();
                        vpl.replace(cmd.getOptionValue("l"), input, out);
                    } else {
                        System.out.println("Use -l with --samplereplace");
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
            } else if (cmd.hasOption("filtersample")) {
                try {
                    VCFSampleFilter filter = new VCFSampleFilter();
                    String samplefile = null;
                    if (cmd.hasOption("l")) {
                        samplefile = cmd.getOptionValue("l");
                        filter.filter(input, out, samplefile);
                    } else {
                        System.out.println("Use params -i -o and -l");
                    }

                } catch (IOException e) {
                    e.printStackTrace();
                }
            } else if (cmd.hasOption("updaters")) {
                try {
                    VCFVariantRSNameUpdater updatr = new VCFVariantRSNameUpdater();
                    if (cmd.hasOption("dbsnp")) {
                        updatr.updateRSNames(cmd.getOptionValue("dbsnp"), input, out);
                    } else {
                        System.out.println("Use --dbsnp, -i and -o with --updaters");
                    }

                } catch (IOException e) {
                    e.printStackTrace();
                }
            } else if (cmd.hasOption("filter")) {


                try {

                    VCFFilter filter = new VCFFilter();
                    Integer gqual = null;
                    Integer readdepth = null;
                    Double allelicBalance = null;
                    double maf = 0;
                    double callrate = 0;
                    boolean onlyautosomes = false;
                    if (cmd.hasOption("readdepth")) {
                        readdepth = Integer.parseInt(cmd.getOptionValue("readdepth"));
                    }
                    if (cmd.hasOption("gqual")) {
                        gqual = Integer.parseInt(cmd.getOptionValue("gqual"));
                    }

                    if (cmd.hasOption("allelicbalance")) {
                        allelicBalance = Double.parseDouble("allelicbalance");
                    }

                    if (cmd.hasOption("maf")) {
                        maf = Double.parseDouble(cmd.getOptionValue("maf"));
                    }

                    if (cmd.hasOption("callrate")) {
                        callrate = Double.parseDouble(cmd.getOptionValue("callrate"));
                    }

                    if (cmd.hasOption("autosomes")) {
                        onlyautosomes = true;
                    }


                    filter.filter(input, out, maf, callrate, readdepth, gqual, allelicBalance, onlyautosomes);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            } else if (cmd.hasOption("plotrsq")) {

                RSquaredPlot plot = new RSquaredPlot();
                try {

                    double threshold = 0.8;
                    if (cmd.hasOption("rsquared")) {
                        threshold = Double.parseDouble(cmd.getOptionValue("rsquared"));
                    }

                    double mafthreshold = 0.0;
                    if (cmd.hasOption("maf")) {
                        mafthreshold = Double.parseDouble(cmd.getOptionValue("maf"));
                    }

                    String bedfile = null;
                    if (cmd.hasOption("b")) {
                        bedfile = cmd.getOptionValue("b");
                    }
                    plot.plot(input, out, threshold, mafthreshold, bedfile);
                } catch (IOException e) {
                    e.printStackTrace();
                } catch (DocumentException e) {
                    e.printStackTrace();
                }

            } else if (cmd.hasOption("plotcor")) {
                CorrelationResultPlotter plotter = new CorrelationResultPlotter();

                try {

                    if (cmd.hasOption("refname")) {
                        plotter.run(input, cmd.getOptionValue("refname"), out);
                    } else {
                        plotter.run(input, out);
                    }


                } catch (IOException e) {
                    e.printStackTrace();
                } catch (DocumentException e) {
                    e.printStackTrace();
                }

            } else if (cmd.hasOption("correlate")) {

                String v2 = null;
                String li = null;

                if (cmd.hasOption("i2")) {
                    v2 = cmd.getOptionValue("i2");
                } else {
                    System.err.println("Provide -i2 with --correlate");
                    run = false;
                }

                if (cmd.hasOption("l")) {
                    li = cmd.getOptionValue("l");
                } else {
                    System.err.println("Provide -l with --correlate");
                    run = false;
                }
                if (run) {
                    VCFCorrelator c = new VCFCorrelator();
                    try {
                        c.run(input, v2, li, out);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                } else {
                    printHelp();
                }
            } else if (cmd.hasOption("rmc")) {
                VCFFunctions f = new VCFFunctions();
                try {
                    f.replaceMissingCodes(input, out);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            } else if (cmd.hasOption("reintroduce")) {
                VCFMerger merger = new VCFMerger();
                boolean linux = false;
                String v2 = "";
                String vcfsort = "";

                if (cmd.hasOption("linux")) {
                    linux = true;
                }
                if (cmd.hasOption("i2")) {
                    v2 = cmd.getOptionValue("i2");
                } else {
                    System.err.println("Provide -i2 with --correlate");
                    run = false;
                }
                if (cmd.hasOption("vcfsort")) {
                    vcfsort = cmd.getOptionValue("vcfsort");
                } else {
                    System.out.println("Please supply --vcfsort");
                }
                if (run) {
                    try {
                        merger.reintroducteNonImputedVariants(input, v2, out, linux, vcfsort);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                } else {
                    printHelp();
                }
            } else if (cmd.hasOption("dedup")) {
                VCFFunctions f = new VCFFunctions();
                try {
                    f.removedups(input, out);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            } else if (cmd.hasOption("combine")) {
                CorrelationResultCombiner combine = new CorrelationResultCombiner();

                String refname = null;
                int nriter = 0;
                int nrbatches = 0;
                run = true;
                if (cmd.hasOption("refname")) {
                    refname = cmd.getOptionValue("refname");
                } else {
                    run = false;
                    System.out.println("Provide --refname");
                }

                if (cmd.hasOption("nriter")) {
                    nriter = Integer.parseInt(cmd.getOptionValue("nriter"));
                } else {
                    run = false;
                    System.out.println("Provide --nriter");
                }

                if (cmd.hasOption("nrbatches")) {
                    nrbatches = Integer.parseInt(cmd.getOptionValue("nrbatches"));
                } else {
                    run = false;
                    System.out.println("Provide --nrbatches");
                }

                if (run) {
                    try {
                        combine.run(input, out, refname, nriter, nrbatches);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            } else if (cmd.hasOption("merge")) {
//				VCFMerger merger = new VCFMerger();
//				merger.mergeAndIntersect(linux,chrint,vcfsort,vcf1,vcf2,out,sep);
            } else if (cmd.hasOption("mergeimputation")) {
                VCFMerger merger = new VCFMerger();
                int nrbatches = 0;
                if (cmd.hasOption("nrbatches")) {
                    nrbatches = Integer.parseInt(cmd.getOptionValue("nrbatches"));
                } else {
                    run = false;
                    System.out.println("Provide --nrbatches");
                }

                String variantList = null;
                if (cmd.hasOption("l")) {
                    variantList = cmd.getOptionValue("l");
                }

                if (run) {
                    try {
                        merger.mergeImputationBatches(input, out, variantList, nrbatches);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            } else if (cmd.hasOption("mergecheese")) {
                VCFMerger merger = new VCFMerger();
                String vcf1 = input;
                if (cmd.hasOption("i2")) {
                    String vcf2 = cmd.getOptionValue("i2");
                    try {
                        merger.merge(vcf1, vcf2, out);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                } else {
                    System.out.println("Please provide -i2 for --mergecheese");
                    printHelp();
                }
            } else if (cmd.hasOption("hapsample")) {
                HapSample s = new HapSample();
                boolean linux = false;

                if (cmd.hasOption("linux")) {
                    linux = true;
                }
                if (cmd.hasOption("vcfsort")) {
                    try {
                        s.covertHapsSampleToVCF(input, out, linux, cmd.getOptionValue("vcfsort"));
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                } else {
                    System.out.println("Please supply --vcfsort");
                    printHelp();
                }
            } else if (cmd.hasOption("split")) {

                int n = 0;
                String ped = null;
                if (cmd.hasOption("f")) {
                    ped = cmd.getOptionValue("f");
                } else {
                    System.out.println("Please provide pedigree data (ped or fam file)");
                    run = false;
                }

                if (cmd.hasOption("n")) {
                    String dstr = cmd.getOptionValue("n");
                    try {
                        n = Integer.parseInt(dstr);
                        if (n <= 0) {
                            System.out.println("Value " + dstr + " provided should be larger than 0");
                            run = false;
                        }
                    } catch (NumberFormatException e) {
                        System.out.println("Value " + dstr + " provided for -n is not an int");
                        run = false;
                    }
                } else {
                    run = false;
                    System.out.println("Please provide batch size");
                }
                Random random = new Random();

                long seed = random.nextLong();
                if (cmd.hasOption("s")) {
                    seed = Long.parseLong(cmd.getOptionValue("s"));
                    System.out.println("Using preset seed: " + seed);
                }

                System.out.println(run);
                if (run) {
                    try {

                        BatchSplitter f = new BatchSplitter();
                        System.out.println(input);
                        System.out.println(ped);
                        System.out.println(out);
                        System.out.println(seed);
                        f.splitVCFOverRandomBatches(input, ped, out, n, seed);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                } else {
                    printHelp();
                }


            } else if (cmd.hasOption("sample")) {
                if (cmd.hasOption("p")) {
                    String dstr = cmd.getOptionValue("p");
                    try {
                        perc = Double.parseDouble(dstr);
                        if (perc >= 1 || perc <= 0) {
                            System.out.println("Value " + dstr + " provided should be between 0 and 1");
                            run = false;
                        }
                    } catch (NumberFormatException e) {
                        System.out.println("Value " + dstr + " provided for -p is not a double.");
                        run = false;
                    }
                } else {
                    run = false;
                    System.out.println("Please provide sampling percentage");
                }

                if (cmd.hasOption("r")) {
                    ref = cmd.getOptionValue("r");
                }


                if (run) {
                    try {
                        VariantSampler sampler = new VariantSampler();
                        sampler.sample(input, ref, perc, out);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                } else {
                    printHelp();
                }
            } else if (cmd.hasOption("splitchr")) {
                VCFFunctions f = new VCFFunctions();
                try {
                    f.splitPerChromosome(input, out);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            } else if (cmd.hasOption("bedfilter")) {
                VCFFunctions f = new VCFFunctions();
                if (cmd.hasOption("b")) {
                    try {
                        f.filterVCFForBedRegions(input, out, cmd.getOptionValue("b"));
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                } else {
                    System.out.println("Please provide -b for --filter");
                    printHelp();
                }

            } else if (cmd.hasOption("pseudo")) {

                String fam = null;

                if (cmd.hasOption("f")) {
                    fam = cmd.getOptionValue("f");
                } else {
                    System.out.println("Please provide -f for --pseudo");
                    printHelp();
                    System.exit(0);
                }

                try {
                    PseudoControls c = new PseudoControls();
                    c.make(input, out, fam, out + ".fam");
                } catch (IOException e) {
                    e.printStackTrace();
                }


            } else {
                printHelp();
            }


        } catch (ParseException ex) {
            printHelp();
        }
    }


    public static void printHelp() {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp(" ", OPTIONS);
        System.exit(-1);
    }

}
