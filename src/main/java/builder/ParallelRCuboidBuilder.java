package builder;

import com.vividsolutions.jts.algorithm.Angle;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.MultiLineString;
import com.vividsolutions.jts.operation.distance.DistanceOp;

import cuboidRoofed.geometry.simple.AbstractParallelCuboidRoofed;
import cuboidRoofed.geometry.simple.ParallelCuboidRoofed;
import cuboidRoofed.geometry.simple.ParallelCuboidRoofed2;
import fr.ign.cogit.geoxygene.api.spatial.geomroot.IGeometry;
import fr.ign.cogit.geoxygene.util.conversion.AdapterFactory;
import fr.ign.cogit.simplu3d.rjmcmc.paramshp.geometry.impl.CuboidRoofed;
import fr.ign.mpp.kernel.ObjectBuilder;

public class ParallelRCuboidBuilder implements ObjectBuilder<CuboidRoofed> {

  private GeometryFactory factory;
  private MultiLineString limits;
  private int bandType;

  public ParallelRCuboidBuilder(IGeometry[] limits, int bandType)
      throws Exception {
    factory = new GeometryFactory();
    LineString[] lineStrings = new LineString[limits.length];
    for (int i = 0; i < limits.length; i++) {
      lineStrings[i] = (LineString) AdapterFactory.toGeometry(factory,
          limits[i]);
    }
    this.limits = factory.createMultiLineString(lineStrings);
    this.bandType = bandType;

  }

  @Override
  public CuboidRoofed build(double[] coordinates) {
    Coordinate p = new Coordinate(coordinates[0], coordinates[1]);
    DistanceOp op = new DistanceOp(this.limits, factory.createPoint(p));
    Coordinate projected = op.nearestPoints()[0];
    double distance = op.distance();
    double orientation = Angle.angle(p, projected);
    AbstractParallelCuboidRoofed result;
    if (bandType == 1) {

      result = new ParallelCuboidRoofed(coordinates[0], coordinates[1],
          coordinates[2], distance * 2, coordinates[3],
          orientation + Math.PI / 2, coordinates[4], coordinates[5]);

    } else {
      result = new ParallelCuboidRoofed2(coordinates[0], coordinates[1],
          coordinates[2], distance * 2, coordinates[3],
          orientation + Math.PI / 2, coordinates[4], coordinates[5]);

    }

    return result;

  }

  @Override
  public void setCoordinates(CuboidRoofed t, double[] coordinates) {
    AbstractParallelCuboidRoofed pc = (AbstractParallelCuboidRoofed) t;
    coordinates[0] = pc.centerx;
    coordinates[1] = pc.centery;
    coordinates[2] = pc.length;
    coordinates[3] = pc.height;
    coordinates[4] = pc.getHeightT();
    coordinates[5] = pc.getDeltaFromSide();
  }

  @Override
  public int size() {
    return 6;
  }

}
