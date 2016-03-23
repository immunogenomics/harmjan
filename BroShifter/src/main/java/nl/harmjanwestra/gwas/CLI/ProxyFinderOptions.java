package nl.harmjanwestra.gwas.CLI;

/**
 * Created by hwestra on 3/22/16.
 */
public class ProxyFinderOptions {
	public String tabixrefprefix;
	public int windowsize = 1000000;
	public double threshold = 0.8;
	public String snpfile;
	public String output;
	public int nrthreads = 1;


	public ProxyFinderOptions(String[] args) {

	}


}
