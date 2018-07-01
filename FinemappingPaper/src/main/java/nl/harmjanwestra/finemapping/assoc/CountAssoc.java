package nl.harmjanwestra.finemapping.assoc;

import nl.harmjanwestra.utilities.association.AssociationFile;
import nl.harmjanwestra.utilities.association.AssociationResult;
import nl.harmjanwestra.utilities.bedfile.BedFileReader;
import nl.harmjanwestra.utilities.features.Feature;

import java.io.IOException;
import java.util.ArrayList;

public class CountAssoc {


    public static void main(String[] args) {

        String[] assocs = new String[]{
                "D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\RA-assoc0.3-COSMO-merged-posterior.txt.gz",
                "D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\T1D-assoc0.3-COSMO-merged-posterior.txt.gz",
                "D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\missp\\META-assoc0.3-COSMO-merged-posterior.txt.gz"
        };

        double t = 7.5E-7;
        String regions = "d:/Sync/SyncThing/Postdoc/2016-03-RAT1D-Finemapping/Data/LocusDefinitions/AllICLoci-overlappingWithImmunobaseT1DOrRALoci.bed";
        for (String f : assocs) {
            CountAssoc c = new CountAssoc();
            try {
                c.run(f, regions, t);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

    }

    public void run(String file, String regionsf, double t) throws IOException {

        AssociationFile f = new AssociationFile();
        ArrayList<AssociationResult> assoc = f.read(file);

        BedFileReader r = new BedFileReader();
        ArrayList<Feature> regions = r.readAsList(regionsf);
        int maxvals = 0;
        for (Feature region : regions) {

            int nrvals = 0;
            boolean sig = false;
            for (AssociationResult rv : assoc) {
                if (rv.getSnp().overlaps(region)) {
                    nrvals++;
                    if (rv.getPval() < t) {
                        sig = true;
                    }
                }
            }
            if (sig) {
                if (nrvals > maxvals) {
                    maxvals = nrvals;
                }
            }
        }

        System.out.println("total: " + assoc.size() + "Max vals: " + maxvals + "\t" + (0.05 / maxvals));
    }
}
