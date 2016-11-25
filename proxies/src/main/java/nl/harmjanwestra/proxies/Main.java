
package nl.harmjanwestra.proxies;

import java.io.IOException;

public class Main {

	public static void main(String[] args) {
		ProxyFinderOptions options = new ProxyFinderOptions(args);
		try {
			ProxyFinder finder = new ProxyFinder(options);
			if (options.locusld) {
				finder.locusLD();
			} else if (options.pairwise) {
				finder.pairwiseLD();
			} else {
				finder.findProxies();
			}
		} catch (IOException e) {
			System.out.println("Unrecoverable error: ");
			e.printStackTrace();
		}

	}

}