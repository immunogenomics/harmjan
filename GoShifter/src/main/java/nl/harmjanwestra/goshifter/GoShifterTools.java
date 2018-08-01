package nl.harmjanwestra.goshifter;

/**
 * Created by hwestra on 2/26/15.
 */
public class GoShifterTools {
	public static void main(String[] args) {

		GoShifterTools tools = new GoShifterTools();
		tools.run(args);


	}

	private void run(String[] args) {
		if (args.length < 2) {
			printUsage();
		} else {
			boolean jobRun = false;
			for (int i = 0; i < args.length; i++) {

				if (args[i].equals("--mode")) {

					if (i + 1 < args.length) {
						String val = args[i + 1];

						if (val.equals("jobs")) {
							GoShifterJobCreator jobs = new GoShifterJobCreator();
							jobs.run(args);
							jobRun = true;
							break;
						} else if (val.equals("summarize")) {

							GoShifterResultSummarizer summarizer = new GoShifterResultSummarizer();
							summarizer.run(args);
							jobRun = true;
							break;
						}
					}
				}
			}

			if (!jobRun) {
				System.out.println("Did you specify a mode?");
				printUsage();
			}
		}
	}

	private void printUsage() {
		System.out.println("Usage: GoShifter --mode jobs/summarize");
	}
}
