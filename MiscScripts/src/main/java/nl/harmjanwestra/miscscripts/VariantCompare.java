package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.vcf.VCFVariant;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by hwestra on 1/13/17.
 */
public class VariantCompare {
	public static void main(String[] args) {

		String[] files = new String[]{
				"/Data/Projects/RachelKnevel/stats/1kg-COSMO-Chr1-stats.vcf.gz",
				"/Data/Projects/RachelKnevel/stats/1kg-EUR-Chr1-stats.vcf.gz",
				"/Data/Projects/RachelKnevel/stats/leiden-Chr1-stats.vcf.gz"
		};
		String[] names = new String[]{
				"COSMO",
				"EUR",
				"Leiden"
		};
		String out = "/Data/Projects/RachelKnevel/stats/out-test.txt";

		VariantCompare c = new VariantCompare();
		try {
			c.run(files, names, out);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String[] files, String[] names, String out) throws IOException {

		ArrayList<ArrayList<VCFVariant>> allVars = new ArrayList<>();

		ArrayList<HashMap<Integer, ArrayList<VCFVariant>>> pointers = new ArrayList<HashMap<Integer, ArrayList<VCFVariant>>>();
		HashSet<Integer> positions = new HashSet<Integer>();
		for (int i = 0; i < files.length; i++) {
			allVars.add(load(files[i]));
			pointers.add(hash(allVars.get(i)));
			Set<Integer> pos = pointers.get(i).keySet();
			positions.addAll(pos);
		}

		String header = "Id";
		for (int i = 0; i < names.length; i++) {
			header += "\tRefAf-" + names[i] + "\tImpQ-" + names[i];
		}

		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln(header);
		for (Integer position : positions) {


			HashSet<String> varsOnPos = new HashSet<>();
			for (int i = 0; i < pointers.size(); i++) {
				HashMap<Integer, ArrayList<VCFVariant>> set = pointers.get(i);
				ArrayList<VCFVariant> vars = set.get(position);
				for (VCFVariant v : vars) {
					varsOnPos.add(v.getPos() + "-" + v.getId() + "-" + Strings.concat(v.getAlleles(), Strings.semicolon));
				}
			}
			int ctr = 0;
			HashMap<String, Integer> map = new HashMap<>();

			for (String s : varsOnPos) {
				map.put(s, ctr);
				ctr++;
			}

			String[][] varstr = new String[map.size()][names.length];
			for (int i = 0; i < pointers.size(); i++) {
				HashMap<Integer, ArrayList<VCFVariant>> set = pointers.get(i);
				ArrayList<VCFVariant> vars = set.get(position);
				for (VCFVariant v : vars) {
					Integer id = map.get(v.getPos() + "-" + v.getId() + "-" + Strings.concat(v.getAlleles(), Strings.semicolon));
					String info = v.getInfo().get("AF");
					String[] infoelems = info.split(",");
					double remain = 1;
					for (String s : infoelems) {
						double d = Double.parseDouble(s);
						remain -= d;
					}

					varstr[id][i] = remain + "\t" + v.getImputationQualityScore();
				}
			}
			for (String s : varsOnPos) {
				Integer id = map.get(s);
				String outstr = s;
				for (int i = 0; i < pointers.size(); i++) {
					String str = varstr[id][i];
					if (str == null) {
						outstr += "\t-\t-";
					} else {
						outstr += "\t" + str;
					}
				}
				outf.writeln(outstr);
			}


		}
		outf.close();


	}

	private HashMap<Integer, ArrayList<VCFVariant>> hash(ArrayList<VCFVariant> vcfVariants) {
		HashMap<Integer, ArrayList<VCFVariant>> output = new HashMap<>();

		for (VCFVariant v : vcfVariants) {
			ArrayList<VCFVariant> q = output.get(v.getPos());
			if (q == null) {
				q = new ArrayList<>();
			}
			q.add(v);
			output.put(v.getPos(), q);
		}

		return output;
	}

	public ArrayList<VCFVariant> load(String file) throws IOException {

		System.out.println("Loading " + file);
		ArrayList<VCFVariant> v = new ArrayList<>();
		TextFile tf = new TextFile(file, TextFile.R);
		String ln = tf.readLine();
		int ctr = 0;
		while (ln != null) {
			if (!ln.startsWith("#")) {
				VCFVariant var = new VCFVariant(ln, VCFVariant.PARSE.HEADER);
				v.add(var);
			}
			ctr++;
			if (ctr % 100000 == 0) {
				System.out.println(ctr + " lines parsed.");
			}
			ln = tf.readLine();
		}
		tf.close();
		return v;
	}
}
