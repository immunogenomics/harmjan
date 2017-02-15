package nl.harmjanwestra.finemapping;

import umcg.genetica.io.text.TextFile;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

/**
 * Created by hwestra on 11/26/16.
 */
public class AnnotationAnnotation {

	public static void main(String[] args) {

		String annotationfile = "/Data/Enhancers/Roadmap/dnase.txt";
		String out = "/Data/Enhancers/Roadmap/dnase-groups.txt";
		String groupAnnotFile = "/Data/Enhancers/ChromHMM/celltypenames.txt";

		AnnotationAnnotation meh = new AnnotationAnnotation();
		try {
			meh.run(annotationfile, groupAnnotFile, out);
			annotationfile = "/Data/Enhancers/Roadmap/h3k4me3.txt";
			out = "/Data/Enhancers/Roadmap/h3k4me3-groups.txt";

			meh.run(annotationfile, groupAnnotFile, out);

			annotationfile = "/Data/Enhancers/ChromHMM/ChromHMMEnhancers.txt";
			out = "/Data/Enhancers/ChromHMM/ChromHMMEnhancers-groups.txt";
			meh.run(annotationfile, groupAnnotFile, out);

			annotationfile = "/Data/Enhancers/ChromHMM/ChromHMMPromotors.txt";
			out = "/Data/Enhancers/ChromHMM/ChromHMMPromotors-groups.txt";
			meh.run(annotationfile, groupAnnotFile, out);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void run(String annotationfile, String groupAnnotFile, String out) throws IOException {


		TextFile tf = new TextFile(groupAnnotFile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, String> eidToName = new HashMap<String, String>();
		HashMap<String, String> eidToGroup = new HashMap<String, String>();

		while (elems != null) {

			String eid = elems[0];
			String name = elems[1];
			String group = elems[2];

			eidToName.put(eid, name);
			eidToGroup.put(eid, group);

			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		TextFile tf2 = new TextFile(annotationfile, TextFile.R);
		TextFile tfout = new TextFile(out, TextFile.W);

		String ln = tf2.readLine();
		while (ln != null) {


			String outln = ln;

			File file = new File(ln);
			String filename = file.getName();


			String eid = filename.substring(0,4);
			String name = eidToName.get(eid);
			String group = eidToGroup.get(eid);

			outln += "\t" + name + "\t" + group;
			tfout.writeln(outln);
			ln = tf2.readLine();
		}
		tfout.close();
		tf2.close();

	}


}
