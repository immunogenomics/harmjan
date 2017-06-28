
package nl.harmjanwestra.proxyfinder;

import java.io.IOException;

public class Main {

	public static void main(String[] args) {
		ProxyFinderOptions options = new ProxyFinderOptions(args);
		try {
			ProxyFinder finder = new ProxyFinder(options);
			if (options.mode.equals(ProxyFinderOptions.MODE.LOCUSLD)) {
				finder.locusLD();
			} else if (options.mode.equals(ProxyFinderOptions.MODE.PAIRWISE)) {
				finder.pairwiseLD();
			} else if (options.mode.equals(ProxyFinderOptions.MODE.PROXY)) {
				finder.findProxies();
			}
		} catch (IOException e) {
			System.out.println("Unrecoverable error: ");
			e.printStackTrace();
		}

	}

}