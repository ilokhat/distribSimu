package transform;

import fr.ign.rjmcmc.kernel.Transform;

public class RotateRCuboid implements Transform {

  private double amplitudeRotate;

  public RotateRCuboid(double amp) {
    amplitudeRotate = amp;
  }

  @Override
  public int dimension() {
    return 12;
  }

  @Override
  public double apply(boolean direct, double[] val0, double[] val1) {

    double dor = val0[11];
    double newAngle = val0[10] + (0.5 - dor) * amplitudeRotate;
    double modulo = newAngle % (Math.PI);
    if (modulo < 0) {
      modulo = Math.PI + modulo;
    }
    val1[0] = val0[0];
    val1[1] = val0[1];
    val1[2] = val0[2];
    val1[3] = val0[3];
    val1[4] = val0[4];
    val1[5] = val0[5];
    val1[6] = val0[6];
    val1[7] = val0[7];
    val1[8] = val0[8];
    val1[9] = val0[9];
    val1[11] = modulo;

    val1[11] = 1 - dor;
    return 1;
  }

  public double getAbsJacobian(boolean direct) {
    return 1;
  }

}
