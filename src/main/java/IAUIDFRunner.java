
import java.io.File;

import fr.ign.cogit.simplu3d.io.feature.AttribNames;

public class IAUIDFRunner {
  public static double run(File folder, File folderOut, File parameterFile,
      long seed) {

    double result = 0;
    // AttribNames.setATT_CODE_PARC("ID_Parcell");
    AttribNames.setATT_CODE_PARC("IDPAR");
    // AttribNames.setATT_NOM_RUE("NOM_VOIE_G");
    try {
      result = IAUIDFTaskMod.run(folder, folderOut, parameterFile, seed);
    } catch (Exception e) {
      e.printStackTrace();
    }
    return result;
  }

  public static void main(String[] args) throws Exception {
    // AttribNames.setATT_CODE_PARC("ID_Parcell");
    // File folder = new File("/home/imran/jobData/EstEnsemble/75025410/");
    // File folder = new File("/home/imran/jobData/EstEnsemble/75015601/");
    String numrep = "91026042";
    File folder = new File(
        "/home/imran/.openmole/imran-OptiPlex-9010/webui/projects/dataBasicSimu/idf/"
            + numrep + "/");
    File folderOut = new File("/home/imran/testoss/out/" + numrep + "/");
    File parameterFile = new File(
        "/home/imran/jobData/scenario/parameters_iauidf.xml");
    long seed = 42L;
    double res = run(folder, folderOut, parameterFile, seed);

    System.out.println("parcelle coverage : " + res);

  }

}
