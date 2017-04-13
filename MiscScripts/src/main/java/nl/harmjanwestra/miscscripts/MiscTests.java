/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package nl.harmjanwestra.miscscripts;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.Set;

/**
 * @author hwestra
 */
public class MiscTests {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
//
//		File f = new File("/Data/Pipelines/tmp/blaat/test/koediekoedie.txt");
//		System.out.println(f.exists());
//		System.out.println(f.getParent());

		String[] meeh = Strings.tab.split(null);

//        int lastI = 0;
////        for(int i=0;i<123456789; i+=100000){
////            System.out.println(i);
////            lastI = i;
////        }
////        System.out.println(lastI);
////        int a = 0x1;
////        for(int i=0;i<10; i++){
////
//////            System.out.println((a << 8)+"\t"+Integer.toBinaryString(a << 8)+"\t"+Integer.bitCount(a<<8)+"\t"+Integer.);
////        }
////        int v = 1<<5;
//        DecimalFormat f = new DecimalFormat("#.##");
//
//        double a = 0.999;
//        double b = 0.899;
//        double c = 99.999;
//        double d = 999.99;
//        System.out.println(f.format(a));
//        System.out.println(f.format(b));
//        System.out.println(f.format(c));

//		String test = "split|this|string";
//		Pattern p = Pattern.compile("\\|");
//		String[] elems = p.split(test);
//		System.out.println(elems.length);
//		p = Pattern.compile("|");
//		elems = p.split(test);
//		System.out.println(elems.length);

//		HashMap<String, String> q = new HashMap<String, String>();
//		q.put("kaas", "cracker");
//
//		q.remove("kaas");
//		System.out.println(q.containsKey("kaas"));
//
//		try {
//
//			String f1 = "/Data/tmp/snps.txt";
//			String f2 = "/Data/tmp/eqtls.txt";
//			String out = "/Data/tmp/eqtls-repl.txt";
//
//			MiscTests m = new MiscTests();
//			m.filterTransEQTL(f1, f2, out);
//
//
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

		int num = 5;
		int div = 2;
		System.out.println((num/div));
	}

	public void filterTransEQTL(String f1, String f2, String out) throws IOException {
		TextFile tf1 = new TextFile(f1, TextFile.R);
		Set<String> set = tf1.readAsSet(0, TextFile.tab);
		tf1.close();


		TextFile outf = new TextFile(out, TextFile.W);
		TextFile tf2 = new TextFile(f2, TextFile.R);
		outf.writeln(tf2.readLine());
		String[] elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {

			String snp = elems[0];
			if (set.contains(snp)) {
				outf.writeln(Strings.concat(elems, Strings.tab));
			}
			elems = tf2.readLineElems(TextFile.tab);
		}

		outf.close();
		tf2.close();

	}


}
