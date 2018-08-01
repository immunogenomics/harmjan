package nl.harmjanwestra.utilities.shell;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.SequenceInputStream;
import java.util.Map;

/**
 * Created by hwestra on 11/30/15.
 */
public class ProcessRunner {

	public static void run(ProcessBuilder builder) throws IOException {
		Map<String, String> environ = builder.environment();

		final Process process = builder.start();

		java.io.SequenceInputStream is = new SequenceInputStream(process.getInputStream(), process.getErrorStream());

		InputStreamReader isr = new InputStreamReader(is);
		BufferedReader br = new BufferedReader(isr);

		String line;
		while ((line = br.readLine()) != null) {
			System.out.println(line);
		}
	}
}
