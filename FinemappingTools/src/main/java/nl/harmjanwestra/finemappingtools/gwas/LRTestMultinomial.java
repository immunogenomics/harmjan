package nl.harmjanwestra.finemappingtools.gwas;

import nl.harmjanwestra.finemappingtools.gwas.CLI.LRTestOptions;

import java.io.IOException;
import java.util.concurrent.Executors;

/**
 * Created by westr on 07/13/17.
 */
public class LRTestMultinomial extends LRTest {
    public LRTestMultinomial(LRTestOptions options) throws IOException {
        super(options);
        exService = Executors.newWorkStealingPool(options.getNrThreads());
        run();
        exService.shutdown();
    }



    public void run() throws IOException {



    }
}
