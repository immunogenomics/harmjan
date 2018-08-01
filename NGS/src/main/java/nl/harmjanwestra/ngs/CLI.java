/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.harmjanwestra.ngs;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 *
 * @author hwestra
 */
public class CLI {

    private String coverageMatrix;
    private double coverageThreshold = 0.9;
    private int coverageThresholdColumn = 3;
    private String bamFileIndir;
    private String bamFileOutdir;
    private String bamFileSuffix = "bam";

    public CLI(String[] args) {

        Options options = new Options();
        options.addOption("c", true, "Coverage matrix");
        options.addOption("ct", true, "Coverage Threshold");
        options.addOption("cc", true, "Coverage Threshold column");
        options.addOption("i", true, "BAM path indir");
        options.addOption("o", true, "BAM path outdir");
        options.addOption("s", true, "BAM path suffix");

        CommandLineParser parser = new BasicParser();

        try {
            CommandLine cmd = parser.parse(options, args);
            if (cmd.hasOption("c")) {
                coverageMatrix = cmd.getOptionValue("c");
            }
            if (cmd.hasOption("i")) {
                bamFileIndir = cmd.getOptionValue("i");
            }
            if (cmd.hasOption("o")) {
                bamFileOutdir = cmd.getOptionValue("o");
            }
            if (cmd.hasOption("o")) {
                bamFileSuffix = cmd.getOptionValue("s");
            }

            if (cmd.hasOption("ct")) {
                String val = cmd.getOptionValue("ct");
                try {
                    coverageThreshold = Double.parseDouble(val);
                } catch (NumberFormatException e) {
                    System.out.println("Value for -ct is not a double");
                }
            }
            if (cmd.hasOption("cc")) {
                String val = cmd.getOptionValue("cc");
                try {
                    coverageThresholdColumn = Integer.parseInt(val);
                } catch (NumberFormatException e) {
                    System.out.println("Value for -cc is not a double");
                }
            }
        } catch (ParseException ex) {
            System.out.println("CLI parse exception.");
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("CLI", options);
        }

        if (bamFileIndir == null || bamFileOutdir == null || coverageMatrix == null) {
            if (bamFileIndir == null) {
                System.out.println("provide -i");
            }
            if (bamFileOutdir == null) {
                System.out.println("provide -o");
            }

            if (coverageMatrix == null) {
                System.out.println("provide -c");
            }
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("CLI", options);
            System.exit(0);
        }

    }

    public String getCoverageMatrix() {
        return coverageMatrix;
    }

    public double getCoverageThreshold() {
        return coverageThreshold;
    }

    public int getCoverageThresholdColumn() {
        return coverageThresholdColumn;
    }

    public String getBamFileIndir() {
        return bamFileIndir;
    }

    public String getBamFileOutdir() {
        return bamFileOutdir;
    }

    public String getBamSuffix() {
        return bamFileSuffix;
    }
}
