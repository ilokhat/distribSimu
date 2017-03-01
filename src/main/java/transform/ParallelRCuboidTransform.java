package transform;

import org.apache.log4j.Logger;

import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;

import fr.ign.cogit.geoxygene.api.spatial.geomroot.IGeometry;
import fr.ign.cogit.geoxygene.util.conversion.AdapterFactory;
import fr.ign.cogit.simplu3d.rjmcmc.cuboid.transformation.birth.ParallelPolygonTransform;
import fr.ign.geometry.transform.PolygonTransform;
import fr.ign.rjmcmc.kernel.Transform;

public class ParallelRCuboidTransform implements Transform {

  /**
   * Logger.
   */
  static Logger LOGGER = Logger
      .getLogger(ParallelPolygonTransform.class.getName());

  private double absJacobian[];
  private PolygonTransform polygonTransform;
  // private MultiLineString limits;
  private GeometryFactory factory = new GeometryFactory();

  private double deltaLength;
  private double deltaHeight;
  private double rangeLength;
  private double rangeHeight;
  private double deltaHeightT;
  private double rangeHeightT;
  private double deltaDFSide;
  private double rangeDFSide;

  public ParallelRCuboidTransform(double[] d, double[] v, IGeometry polygon)
      throws Exception {
    this.rangeLength = d[2];
    this.rangeHeight = d[4];
    this.rangeHeightT = d[6];
    this.rangeDFSide = d[7];
    this.deltaLength = v[2];
    this.deltaHeight = v[4];
    this.deltaHeightT = v[6];
    this.deltaDFSide = v[7];

    double determinant = rangeLength * rangeHeight * rangeHeightT * rangeDFSide;

    Geometry pp = AdapterFactory.toGeometry(factory, polygon);
    this.polygonTransform = new PolygonTransform(pp, 0.1);
    this.absJacobian = new double[2];
    this.absJacobian[0] = Math.abs(determinant)
        * this.polygonTransform.getAbsJacobian(true);
    this.absJacobian[1] = Math.abs(1 / determinant)
        * this.polygonTransform.getAbsJacobian(false);

  }

  @Override
  public double apply(boolean direct, double[] val0, double[] val1) {
    double pt = this.polygonTransform.apply(direct, val0, val1);
    if (direct) {
      // Coordinate p = new Coordinate(val1.get(0), val1.get(1));
      // DistanceOp op = new DistanceOp(this.limits, factory.createPoint(p));
      // Coordinate projected = op.nearestPoints()[0];
      // double distance = op.distance();
      // double orientation = Angle.angle(p, projected);
      val1[2] = val0[2] * rangeLength + deltaLength;
      val1[3] = val0[3] * rangeHeight + deltaHeight;
      val1[4] = val0[4] * rangeHeightT + deltaHeightT;
      val1[5] = val0[5] * rangeDFSide + deltaDFSide;

      // val1[4] = val0[4] * rangeHeight + deltaHeight;
      // val1[6] = val0[6] * rangeHeightT + deltaHeightT;
      // val1[7] = val0[7] * rangeDFSide + deltaDFSide;

      // val1.set(5, orientation + Math.PI / 2);
      return pt * this.absJacobian[0];
    } else {
      val1[2] = (val0[2] - deltaLength) / rangeLength;
      val1[3] = (val0[3] - deltaHeight) / rangeHeight;
      val1[4] = (val0[4] - deltaHeightT) / rangeHeightT;
      val1[5] = (val0[5] - deltaHeight) / rangeHeight;

      // val1[4] = (val0[4] - deltaHeight) / rangeHeight;
      // val1[6] = (val0[7] - deltaHeightT) / rangeHeightT;
      // val1[7] = (val0[7] - deltaHeight) / rangeHeight;
      // var1.set(4, 0.0);
      // var1.set(5, 0.0);
      return pt * this.absJacobian[1];
    }
  }

  public double getAbsJacobian(boolean direct) {
    return this.absJacobian[direct ? 0 : 1];
  }

  @Override
  public int dimension() {
    return 6;
  }

}
