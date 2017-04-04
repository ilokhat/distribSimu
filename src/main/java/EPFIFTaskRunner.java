import java.io.File;

import fr.ign.cogit.simplu3d.experiments.openmole.EPFIFTask;
import fr.ign.cogit.simplu3d.io.feature.AttribNames;

public class EPFIFTaskRunner {
  public static String run(File folder, File folderOut, File parameterFile,
      long seed) {
    AttribNames.setATT_CODE_PARC("IDPAR");
    EPFIFTask.USE_DEMO_SAMPLER = false;
    EPFIFTask.INTERSECTION = true;
    EPFIFTask.FLOOR_SIZE = 3;

    String result = "";
    try {
      result = EPFIFTask.run(folder, folderOut, parameterFile, seed);
    } catch (Exception e) {
      result += e.toString() + "\n";
      for (StackTraceElement s : e.getStackTrace())
        result += s.toString() + "\n";
      System.out.println("el klodo");
      // e.printStackTrace();
    }
    return result;
  }

  public static void main(String[] args) {
    String numrep = "77017278";// "75033409";// "77017278";
    String foldName = "/home/imran/Téléchargements/Test_IAUIDF/Eval_EPF_2/";
    foldName = "/home/imran/.openmole/imran-OptiPlex-9010/webui/projects/dataBasicSimu/idf/";
    File folder = new File(foldName + numrep + "/");
    File folderOut = new File("/home/imran/testoss/out/" + numrep + "/");
    File parameterFile = new File(
        "/home/imran/testoss/EPFIF/parameters_iauidf.xml");
    long seed = 42L;
    String res = "";
    res = run(folder, folderOut, parameterFile, seed);
    System.out.println("result : " + res);
  }
}
