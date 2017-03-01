package sampler;

import org.apache.commons.math3.random.RandomGenerator;

import fr.ign.cogit.simplu3d.rjmcmc.paramshp.geometry.impl.CuboidRoofed;
import fr.ign.mpp.kernel.ObjectBuilder;
import fr.ign.mpp.kernel.ObjectSampler;
import fr.ign.rjmcmc.kernel.Transform;
import fr.ign.rjmcmc.kernel.Variate;

public class MixCuboidRoofedSampler implements ObjectSampler<CuboidRoofed> {

  RandomGenerator engine;
  double p_simple;
  CuboidRoofed object;
  Variate variate;
  Transform transformBand2;
  Transform transformBand1;
  ObjectBuilder<CuboidRoofed> builder1;
  ObjectBuilder<CuboidRoofed> builder2;
  Class<?> c1;
  Class<?> c2;

  public MixCuboidRoofedSampler(RandomGenerator e, double p_simple,
      Transform transformBand1, Transform transformBand2,
      ObjectBuilder<CuboidRoofed> builder1,
      ObjectBuilder<CuboidRoofed> builder2, Class<?> c1, Class<?> c2) {
    this.engine = e;
    this.p_simple = p_simple;
    this.transformBand2 = transformBand2;
    this.transformBand1 = transformBand1;
    this.variate = new Variate(e);
    this.builder1 = builder1;
    this.builder2 = builder2;
    this.c1 = c1;
    this.c2 = c2;
  }

  @Override
  public double sample(RandomGenerator e) {
    double[] val0;
    double[] val1;
    if (engine.nextDouble() < p_simple) {
      val0 = new double[builder2.size()];
      val1 = new double[builder2.size()];
      double phi = this.variate.compute(val0, 0);
      double jacob = this.transformBand2.apply(true, val0, val1);
      this.object = builder2.build(val1);
      return phi / jacob;
    }
    val0 = new double[builder1.size()];
    val1 = new double[builder1.size()];
    double phi = this.variate.compute(val0, 0);
    double jacob = this.transformBand1.apply(true, val0, val1);
    this.object = builder1.build(val1);
    return phi / jacob;
  }

  @Override
  public CuboidRoofed getObject() {
    return this.object;
  }

  @Override
  public double pdf(CuboidRoofed t) {
    if (c2 != null && c2.isInstance(t)) {
      double[] val1 = new double[builder2.size()];
      builder2.setCoordinates(t, val1);
      double[] val0 = new double[builder2.size()];
      double J10 = this.transformBand2.apply(false, val1, val0);
      double pdf = this.variate.pdf(val0, 0);
      return pdf * J10;
    }

    double[] val1 = new double[builder1.size()];
    builder1.setCoordinates(t, val1);
    double[] val0 = new double[builder1.size()];
    double J10 = this.transformBand1.apply(false, val1, val0);
    if (J10 == 0) {
      return 0;
    }
    double pdf = this.variate.pdf(val0, 0);
    return pdf * J10;

  }

}
