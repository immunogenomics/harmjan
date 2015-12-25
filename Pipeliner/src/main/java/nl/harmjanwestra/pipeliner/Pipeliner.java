package nl.harmjanwestra.pipeliner;


import org.apache.commons.cli.*;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;



/**
 * Created by hwestra on 11/25/15.
 */
public class Pipeliner {

	private static Options OPTIONS;

	static {
		OPTIONS = new Options();


		Option option;

		option = Option.builder("k")
				.desc("File with keys and values")
				.longOpt("keyfile")
				.hasArg()
				.required()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("t")
				.desc("Template")
				.longOpt("template")
				.required()
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("f")
				.desc("Require all keys to be present in filename template")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder("c")
				.desc("Use a counter in stead of filename template")
				.build();
		OPTIONS.addOption(option);


		option = Option.builder("o")
				.desc("Outdir")
				.longOpt("out")
				.required()
				.hasArg()
				.build();
		OPTIONS.addOption(option);
	}


	public static void main(String[] args) {
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, true);

			if (cmd.hasOption("o") && cmd.hasOption("t") && cmd.hasOption("k")) {
				String out = cmd.getOptionValue("o");
				String tpl = cmd.getOptionValue("t");
				String key = cmd.getOptionValue("k");


				boolean requireAllVarsInFileNametemplate = false;
				boolean useCounterAsFileName = false;
				if (cmd.hasOption("f")) {
					requireAllVarsInFileNametemplate = true;
				}
				if (cmd.hasOption("c")) {
					useCounterAsFileName = true;
				}

				Pipeliner p = new Pipeliner();
				p.run(key, tpl, out, requireAllVarsInFileNametemplate, useCounterAsFileName);
			} else {
				printHelp();
			}


		} catch (ParseException ex) {
			System.out.println(ex.getMessage());

			printHelp();
		}

//		nl.harmjanwestra.pipeliner.Pipeliner p = new nl.harmjanwestra.pipeliner.Pipeliner();
//		String key = "/Data/Pipelines/liftoverIC.txt";
//		String tpl = "/Data/Pipelines/liftoverIC.sh";
//		String out = "/Data/Pipelines/liftoverIC/";
//		boolean requireAllVarsInFileNametemplate = false;
//		p.run(key, tpl, out, requireAllVarsInFileNametemplate);


		System.exit(0);
	}


	public static void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);
		System.exit(-1);
	}

	public void run(String keyfile, String templateFile, String dirout, boolean requireAllVarsInFileNametemplate, boolean useCounterAsFilename) {

		System.out.println("Key: " + keyfile);
		System.out.println("Template: " + templateFile);
		System.out.println("Dirout: " + dirout);

		try {
			boolean continueRun = true;
			if (!Gpio.exists(keyfile)) {
				System.out.println("Key file: " + keyfile + " not found.");
				continueRun = false;
			}
			if (!Gpio.exists(templateFile)) {
				System.out.println("Template file: " + templateFile + " not found.");
				continueRun = false;
			}

			if (continueRun) {
				Gpio.createDir(dirout);

				Triple<String, ArrayList<String>, HashMap<String, ArrayList<String>>> p = readKeysAndValues(keyfile);
				String combotemplate = p.getLeft();
				ArrayList<String> keys = p.getMiddle();
				HashMap<String, ArrayList<String>> replacements = p.getRight();
				ArrayList<String> template = readTemplate(templateFile);

				keys = checkKeys(template, keys);

				if (requireAllVarsInFileNametemplate) {
					ArrayList<String> tmp = new ArrayList<String>();
					tmp.add(combotemplate);
					keys = checkKeys(tmp, keys);
				}

				if (keys.isEmpty()) {
					System.err.println("Error: no keys remain");
				} else {
					System.out.println(keys.size() + " keys will be used.");

					// check whether all keys are in template
					Pair<ArrayList<String>, ArrayList<ArrayList<String>>> combos = combine(keys, replacements, template, combotemplate);

					Gpio.createDir(dirout);
					ArrayList<String> combonames = combos.getLeft();
					ArrayList<ArrayList<String>> pipelines = combos.getRight();


					for (int i = 0; i < combonames.size(); i++) {
						String name = combonames.get(i);

						if (useCounterAsFilename) {
							name = "" + (i + 1);
						}

						ArrayList<String> pipeline = pipelines.get(i);

						String fullnameOut = dirout + name + ".sh";
						File f = new File(fullnameOut);
						String parent = f.getParent();

						if (!Gpio.exists(parent)) {
							Gpio.createDir(parent);
							System.out.println(parent);
						}
						System.out.println("Writing: " + dirout + name + ".sh");
						writePipeline(dirout + name + ".sh", pipeline);
						f.setExecutable(true, false);

					}

				}
			}


		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private void writePipeline(String outfile, ArrayList<String> pipeline) throws IOException {
		TextFile out = new TextFile(outfile, TextFile.W);
		for (String s : pipeline) {
			out.writeln(s);
		}
		out.close();
	}

	private Pair<ArrayList<String>, ArrayList<ArrayList<String>>> combine(ArrayList<String> keys,
																		  HashMap<String, ArrayList<String>> replacements,
																		  ArrayList<String> template,
																		  String combotemplate) {
		ArrayList<ArrayList<String>> allcombos = new ArrayList<ArrayList<String>>();
		ArrayList<String> combonames = new ArrayList<String>();

		for (int k = 0; k < keys.size(); k++) {
			String key = keys.get(k);
			ArrayList<String> values = replacements.get(key);
			if (!allcombos.isEmpty()) {
				ArrayList<ArrayList<String>> allcombosTmp = new ArrayList<ArrayList<String>>();
				ArrayList<String> comboNamesTmp = new ArrayList<String>();
				for (int c = 0; c < allcombos.size(); c++) {

					for (int v = 0; v < values.size(); v++) {
						ArrayList<String> tpl = allcombos.get(c);
						String comboName = combonames.get(c);

						String val = values.get(v);
						tpl = update(tpl, key, val);
						comboName = update(comboName, key, val);
						allcombosTmp.add(tpl);
						comboNamesTmp.add(comboName);
					}
				}
				allcombos = allcombosTmp;
				combonames = comboNamesTmp;
			} else {
				if (values == null) {
					System.out.println("No replacements for key: " + key);
					System.exit(-1);
				}
				for (int v = 0; v < values.size(); v++) {
					String val = values.get(v);
					ArrayList<String> tpl = update(template, key, val);
					String comboTemplateCopy = update(combotemplate, key, val);
					combonames.add(comboTemplateCopy);
					allcombos.add(tpl);
				}
			}
		}
		return new Pair<ArrayList<String>, ArrayList<ArrayList<String>>>(combonames, allcombos);
	}

	private ArrayList<String> update(ArrayList<String> template, String key, String val) {
		ArrayList<String> updated = new ArrayList<String>();
		for (int i = 0; i < template.size(); i++) {
			updated.add(update(template.get(i), key, val));
		}
		return updated;
	}

	private String update(String template, String key, String val) {
		String replacedLn = template.replaceAll("\\{" + key + "\\}", val);
		return replacedLn;
	}

	private ArrayList<String> checkKeys(ArrayList<String> template, ArrayList<String> keys) {
		ArrayList<String> output = new ArrayList<String>();

		for (String key : keys) {
			boolean present = false;
			Pattern p = Pattern.compile(".*\\{" + key + "\\}.*");
			for (String ln : template) {
				if (p.matcher(ln).matches()) {
					present = true;
				} else {
					// System.out.println(ln + " does not contain " + key);
				}
			}
			if (present) {
				output.add(key);
				System.out.println(key + " found in template");
			} else {
				System.out.println(key + " not found in template");
			}
		}
		return output;
	}

	private Triple<String, ArrayList<String>, HashMap<String, ArrayList<String>>> readKeysAndValues(String f) throws IOException {
		TextFile tf = new TextFile(f, TextFile.R);
		String combotemplate = tf.readLine();
		ArrayList<String> keys = new ArrayList<String>();
		HashMap<String, ArrayList<String>> values = new HashMap<String, ArrayList<String>>();
		String[] elems = tf.readLineElems(Strings.whitespace);
		while (elems != null) {
			if (elems.length > 1) {
				String key = elems[0];
				keys.add(key);
				ArrayList<String> vals = new ArrayList<String>();
				for (int i = 1; i < elems.length; i++) {
					vals.add(elems[i]);
				}
				System.out.println("Key " + key + " has " + vals.size() + " values");
				values.put(key, vals);
			}
			elems = tf.readLineElems(Strings.whitespace);
		}

		tf.close();
		return new Triple<String, ArrayList<String>, HashMap<String, ArrayList<String>>>(combotemplate, keys, values);
	}

	private ArrayList<String> readTemplate(String filename) throws IOException {
		ArrayList<String> output = new ArrayList<String>();
		TextFile tf = new TextFile(filename, TextFile.R);
		String ln = tf.readLine();
		while (ln != null) {
			output.add(ln);
			ln = tf.readLine();
		}
		tf.close();
		return output;
	}
}
