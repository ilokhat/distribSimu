package transform;

import fr.ign.rjmcmc.kernel.Transform;

public class MoveParallelRCuboid implements Transform {

  private double amplitudeMove;

  public MoveParallelRCuboid(double amplitudeMove) {
    this.amplitudeMove = amplitudeMove;
  }

  @Override
  public int dimension() {
    return 8;
  }

  @Override
  public double apply(boolean direct, double[] val0, double[] val1) {
    double dx = val0[6];
    double dy = val0[7];
    val1[0] = val0[0] + (0.5 - dx) * amplitudeMove;
    val1[1] = val0[1] + (0.5 - dy) * amplitudeMove;
    val1[2] = val0[2];
    val1[3] = val0[3];
    val1[4] = val0[4];
    val1[5] = val0[5];

    val1[6] = 1 - dx;
    val1[7] = 1 - dy;
    return 1;

  }

}
