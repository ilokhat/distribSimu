package optimizer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.random.RandomGenerator;

import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;

import builder.CuboidRoofedBuilder2;
import builder.ParallelRCuboidBuilder;
import cuboidRoofed.geometry.impl.CuboidRoofed2;
import cuboidRoofed.geometry.simple.ParallelCuboidRoofed;
import cuboidRoofed.geometry.simple.ParallelCuboidRoofed2;
import fr.ign.cogit.geoxygene.api.spatial.coordgeom.IEnvelope;
import fr.ign.cogit.geoxygene.api.spatial.geomaggr.IMultiCurve;
import fr.ign.cogit.geoxygene.api.spatial.geomprim.IOrientableCurve;
import fr.ign.cogit.geoxygene.api.spatial.geomroot.IGeometry;
import fr.ign.cogit.geoxygene.convert.FromGeomToLineString;
import fr.ign.cogit.geoxygene.spatial.geomaggr.GM_MultiCurve;
import fr.ign.cogit.geoxygene.util.conversion.AdapterFactory;
import fr.ign.cogit.simplu3d.experiments.iauidf.predicate.PredicateIAUIDF;
import fr.ign.cogit.simplu3d.experiments.iauidf.regulation.Regulation;
import fr.ign.cogit.simplu3d.experiments.iauidf.tool.BandProduction;
import fr.ign.cogit.simplu3d.model.BasicPropertyUnit;
import fr.ign.cogit.simplu3d.model.Environnement;
import fr.ign.cogit.simplu3d.model.ParcelBoundary;
import fr.ign.cogit.simplu3d.rjmcmc.cuboid.transformation.birth.TransformToSurface;
import fr.ign.cogit.simplu3d.rjmcmc.generic.energy.IntersectionVolumeBinaryEnergy;
import fr.ign.cogit.simplu3d.rjmcmc.generic.energy.VolumeUnaryEnergy;
import fr.ign.cogit.simplu3d.rjmcmc.generic.optimizer.DefaultSimPLU3DOptimizer;
import fr.ign.cogit.simplu3d.rjmcmc.generic.sampler.GreenSamplerBlockTemperature;
import fr.ign.cogit.simplu3d.rjmcmc.generic.visitor.PrepareVisitors;
import fr.ign.cogit.simplu3d.rjmcmc.paramshp.builder.CuboidRoofedBuilder;
import fr.ign.cogit.simplu3d.rjmcmc.paramshp.geometry.impl.CuboidRoofed;
import fr.ign.cogit.simplu3d.rjmcmc.paramshp.transform.MoveRCuboid;
import fr.ign.mpp.DirectRejectionSampler;
import fr.ign.mpp.DirectSampler;
import fr.ign.mpp.configuration.BirthDeathModification;
import fr.ign.mpp.configuration.GraphConfiguration;
import fr.ign.mpp.kernel.ObjectBuilder;
import fr.ign.mpp.kernel.UniformTypeView;
import fr.ign.parameters.Parameters;
import fr.ign.random.Random;
import fr.ign.rjmcmc.acceptance.Acceptance;
import fr.ign.rjmcmc.acceptance.MetropolisAcceptance;
import fr.ign.rjmcmc.configuration.ConfigurationModificationPredicate;
import fr.ign.rjmcmc.distribution.PoissonDistribution;
import fr.ign.rjmcmc.energy.BinaryEnergy;
import fr.ign.rjmcmc.energy.ConstantEnergy;
import fr.ign.rjmcmc.energy.MinusUnaryEnergy;
import fr.ign.rjmcmc.energy.MultipliesBinaryEnergy;
import fr.ign.rjmcmc.energy.MultipliesUnaryEnergy;
import fr.ign.rjmcmc.energy.UnaryEnergy;
import fr.ign.rjmcmc.kernel.ChangeValue;
import fr.ign.rjmcmc.kernel.Kernel;
import fr.ign.rjmcmc.kernel.NullView;
import fr.ign.rjmcmc.kernel.Transform;
import fr.ign.rjmcmc.kernel.Variate;
import fr.ign.rjmcmc.sampler.Sampler;
import fr.ign.simulatedannealing.SimulatedAnnealing;
import fr.ign.simulatedannealing.endtest.EndTest;
import fr.ign.simulatedannealing.schedule.Schedule;
import fr.ign.simulatedannealing.temperature.SimpleTemperature;
import fr.ign.simulatedannealing.visitor.CompositeVisitor;
import sampler.MixCuboidRoofedSampler;
import transform.MoveParallelRCuboid;
import transform.ParallelRCuboidTransform;
import transform.RotateRCuboid;

public class MultipleBuildingsRCuboid
    extends DefaultSimPLU3DOptimizer<CuboidRoofed> {

  public static boolean ALLOW_INTERSECTING_CUBOID = false;

  public GraphConfiguration<CuboidRoofed> process(BasicPropertyUnit bpu,
      Parameters p, Environnement env,
      PredicateIAUIDF<CuboidRoofed, GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> pred,
      Regulation r1, Regulation r2, BandProduction bP) throws Exception {

    // Géométrie de l'unité foncière sur laquelle porte la génération (on se
    // permet de faire un petit buffer)
    IGeometry geom = bpu.getpol2D().buffer(1);

    // Définition de la fonction d'optimisation (on optimise en décroissant)
    // relative au volume
    GraphConfiguration<CuboidRoofed> conf = null;

    try {
      conf = create_configuration(p, geom, bpu);
    } catch (Exception e) {
      e.printStackTrace();
    }

    // Création de l'échantilloneur
    Sampler<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> samp = create_sampler(
        Random.random(), p, bpu, pred, r1, r2, bP);
    if (samp == null) {
      return null;
    }
    // Température
    Schedule<SimpleTemperature> sch = create_schedule(p);

    EndTest end = create_end_test(p);

    PrepareVisitors<CuboidRoofed> pv = new PrepareVisitors<>(env);
    CompositeVisitor<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> mVisitor = pv
        .prepare(p, bpu.getId());
    /*
     * < This is the way to launch the optimization process. Here, the magic
     * happen... >
     */
    SimulatedAnnealing.optimize(Random.random(), conf, samp, sch, end,
        mVisitor);
    return conf;
  }

  public GraphConfiguration<CuboidRoofed> create_configuration(Parameters p,
      IGeometry geom, BasicPropertyUnit bpu) throws Exception {
    return this.create_configuration(p,
        AdapterFactory.toGeometry(new GeometryFactory(), geom), bpu);
  }

  // Création de la configuration
  /**
   * @param p paramètres importés depuis le fichier XML
   * @param bpu l'unité foncière considérée
   * @return la configuration chargée, c'est à dire la formulation énergétique
   *         prise en compte
   */
  public GraphConfiguration<CuboidRoofed> create_configuration(Parameters p,
      Geometry geom, BasicPropertyUnit bpu) {

    if (ALLOW_INTERSECTING_CUBOID) {
      return create_configuration_intersection(p, geom, bpu);
    }

    return create_configuration_no_inter(p, geom, bpu);

  }

  private GraphConfiguration<CuboidRoofed> create_configuration_intersection(
      Parameters p, Geometry geom, BasicPropertyUnit bpu) {
    // Énergie constante : à la création d'un nouvel objet
    ConstantEnergy<CuboidRoofed, CuboidRoofed> energyCreation = new ConstantEnergy<CuboidRoofed, CuboidRoofed>(
        p.getDouble("energy"));
    // Énergie constante : pondération de l'intersection
    ConstantEnergy<CuboidRoofed, CuboidRoofed> ponderationVolume = new ConstantEnergy<CuboidRoofed, CuboidRoofed>(
        p.getDouble("ponderation_volume"));
    // Énergie unaire : aire dans la parcelle
    UnaryEnergy<CuboidRoofed> energyVolume = new VolumeUnaryEnergy<CuboidRoofed>();
    // Multiplication de l'énergie d'intersection et de l'aire
    UnaryEnergy<CuboidRoofed> energyVolumePondere = new MultipliesUnaryEnergy<CuboidRoofed>(
        ponderationVolume, energyVolume);

    // On retire de l'énergie de création, l'énergie de l'aire
    UnaryEnergy<CuboidRoofed> u3 = new MinusUnaryEnergy<CuboidRoofed>(
        energyCreation, energyVolumePondere);

    // Énergie binaire : intersection entre deux rectangles
    ConstantEnergy<CuboidRoofed, CuboidRoofed> c3 = new ConstantEnergy<CuboidRoofed, CuboidRoofed>(
        p.getDouble("ponderation_volume_inter"));
    BinaryEnergy<CuboidRoofed, CuboidRoofed> b1 = new IntersectionVolumeBinaryEnergy<CuboidRoofed>();
    BinaryEnergy<CuboidRoofed, CuboidRoofed> binaryEnergy = new MultipliesBinaryEnergy<CuboidRoofed, CuboidRoofed>(
        c3, b1);
    // empty initial configuration*/
    GraphConfiguration<CuboidRoofed> conf = new GraphConfiguration<>(u3,
        binaryEnergy);
    return conf;
  }

  private GraphConfiguration<CuboidRoofed> create_configuration_no_inter(
      Parameters p, Geometry geom, BasicPropertyUnit bpu) {
    // Énergie constante : à la création d'un nouvel objet
    double energyCrea = p.getDouble("energy");
    ConstantEnergy<CuboidRoofed, CuboidRoofed> energyCreation = new ConstantEnergy<CuboidRoofed, CuboidRoofed>(
        energyCrea);
    // Énergie constante : pondération de l'intersection
    ConstantEnergy<CuboidRoofed, CuboidRoofed> ponderationVolume = new ConstantEnergy<CuboidRoofed, CuboidRoofed>(
        p.getDouble("ponderation_volume"));
    // Énergie unaire : aire dans la parcelle
    UnaryEnergy<CuboidRoofed> energyVolume = new VolumeUnaryEnergy<CuboidRoofed>();
    // Multiplication de l'énergie d'intersection et de l'aire
    UnaryEnergy<CuboidRoofed> energyVolumePondere = new MultipliesUnaryEnergy<CuboidRoofed>(
        ponderationVolume, energyVolume);
    // On retire de l'énergie de création, l'énergie de l'aire
    UnaryEnergy<CuboidRoofed> u3 = new MinusUnaryEnergy<CuboidRoofed>(
        energyCreation, energyVolumePondere);
    // empty initial configuration*/
    GraphConfiguration<CuboidRoofed> conf = new GraphConfiguration<>(u3,
        new ConstantEnergy<CuboidRoofed, CuboidRoofed>(0));
    return conf;
  }

  public static Sampler<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> create_sampler(
      RandomGenerator rng, Parameters p, BasicPropertyUnit bpU,
      ConfigurationModificationPredicate<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> pred,
      Regulation r1, Regulation r2, BandProduction bP) throws Exception {

    ////////////////////////
    // On prépare les paramètres de tirage aléatoire des boîtes
    ////////////////////////

    // On récupère les intervalles dans lesquels on va tirer aléatoirement
    // les carac des boîtes
    double minlen = p.getDouble("minlen");
    double maxlen = p.getDouble("maxlen");

    double minwid = p.getDouble("minwid");
    double maxwid = p.getDouble("maxwid");

    double minheight = p.getDouble("minheightG");
    double maxheight = p.getDouble("maxheightG");
    double maxheight2 = p.getDouble("maxheightG");

    double minheightT = p.getDouble("minheightT");
    double maxheightT = p.getDouble("maxheightT");
    double maxheightT2 = p.getDouble("maxheightT");

    double minDFS = p.getDouble("minDelta");
    double maxDFS = p.getDouble("maxDelta");

    IEnvelope env = bpU.getGeom().envelope();
    // in multi object situations, we need an object builder for each
    // subtype and a sampler for the supertype (end of file)

    // Est-ce que l'implémentation dans la première bande se fait
    // parallèlement à la voirie ?
    // On regarde si c'est possible avec la présence d'une limite donnant
    // sur la voirie sinon c'est non
    boolean band1Parallel = !(bP.getLineRoad() == null
        || bP.getLineRoad().isEmpty());

    // On prépare le vecteur dans lequel on va tirer aléatoirement les
    // boîtes dans la première bande
    double[] v = new double[] { env.minX(), env.minY(), minlen, minwid,
        minheight, 0., minheightT, minDFS };

    // On regarde si la contrainte de hauteur ne permet pas de réduire
    // l'intervallle des hauteurs
    if (r1 != null && r1.getArt_102() != 99) {
      maxheight = Math.min(maxheight, r1.getArt_102());
    }

    double[] d = new double[] { env.maxX(), env.maxY(), maxlen, maxwid,
        maxheight, Math.PI, maxheightT, maxDFS };

    if (r2 != null && r2.getArt_102() != 99) {
      maxheight2 = Math.min(maxheight2, r2.getArt_102());
    }

    // On répète la même chose pour la seconde
    double[] d2 = new double[] { env.maxX(), env.maxY(), maxlen, maxwid,
        maxheight2, Math.PI, maxheightT2, maxDFS };

    for (int i = 0; i < d.length; i++) {
      d[i] = d[i] - v[i];
    }

    for (int i = 0; i < d.length; i++) {
      d2[i] = d2[i] - v[i];
    }

    System.out.println(Arrays.toString(d));
    System.out.println(Arrays.toString(d2));
    System.out.println(Arrays.toString(v));
    ////////////////////////
    // On prépare les zones dans lesquelles les boîtes seront tirées
    ////////////////////////

    // On récupère la bande numéro 1
    IGeometry geomBand1 = r1.getGeomBande();
    IGeometry geomBand2 = null;

    ////////////////////////
    // On prépare les transforme
    ////////////////////////

    // On calcule la transforme 1 => il n'est pas initialisé s'il n'y a pas
    // de bande 1
    Transform transformBand1 = null;
    ObjectBuilder<CuboidRoofed> builderBand1 = null;
    Class<?> c1 = null;

    if (geomBand1 != null && !geomBand1.isEmpty()) {
      if (band1Parallel) {
        // S'il n'y a qu'une seule bande de constructibilité
        // On peut demander à construire des bâtiments dans la bande
        // derrière
        // le bâtiment aligné à la voirie
        // On va séparer en 2 en fonction de la largeur max du bâtiment
        if (r2 == null) {
          geomBand2 = geomBand1.difference(bP.getLineRoad().buffer(maxwid));

          // Si la bande est toute petite alors, on ne met rien
          if (geomBand2.area() < 5) {
            geomBand2 = null;
          }

        }

        // The center is included in a band equals to half of max
        // allowed width according to alignment line

        geomBand1 = geomBand1.intersection(bP.getLineRoad().buffer(maxwid / 2));

        builderBand1 = new ParallelRCuboidBuilder(bP.getLineRoad().toArray(),
            1);
        transformBand1 = new ParallelRCuboidTransform(d2, v, geomBand1);
        c1 = ParallelCuboidRoofed.class;
      } else {

        geomBand1 = geomBand1.buffer(-minwid / 2);

        transformBand1 = new TransformToSurface(d2, v, geomBand1);
        builderBand1 = new CuboidRoofedBuilder();
        c1 = CuboidRoofed.class;// SimpleCuboid.class;
      }
    }
    ObjectBuilder<CuboidRoofed> builderBand2 = null;

    boolean band2parallel = false;
    Class<?> c2 = null;

    // On calcule la transforme 2 => il n'est pas initialisé s'il n'y a pas
    // de bande 2
    Transform transformBand2 = null;

    // On récupère la seconde bande

    if (r2 != null) {
      geomBand2 = r2.getGeomBande();
    }

    // On calcule la transforme 2 et le builder 2
    if (r2 != null && geomBand2 != null && !geomBand2.isEmpty()) {

      if (r2 != null && r2.getArt_71() == 2) {

        band2parallel = true;

        List<ParcelBoundary> featC = bpU.getCadastralParcels().get(0)
            .getBoundariesBySide(PredicateIAUIDF.RIGHT_OF_LEFT_FOR_ART_71);
        IMultiCurve<IOrientableCurve> ims = new GM_MultiCurve<>();
        for (ParcelBoundary s : featC) {
          ims.addAll(FromGeomToLineString.convert(s.getGeom()));
        }

        // On se colle partout si on peut pas déterminer de côté.
        if (ims.isEmpty()) {

          featC = bpU.getCadastralParcels().get(0).getBoundaries();

          for (ParcelBoundary s : featC) {
            ims.addAll(FromGeomToLineString.convert(s.getGeom()));
          }

        }

        // The center is included in a band equals to half of max
        // allowed width according to alignment line
        geomBand2 = geomBand2.intersection(ims.buffer(maxwid / 2));

        builderBand2 = new ParallelRCuboidBuilder(ims.toArray(), 2);
        transformBand2 = new ParallelRCuboidTransform(d, v, geomBand2);
        c2 = ParallelCuboidRoofed2.class;

      } else {

        geomBand2 = geomBand2.buffer(-minwid / 2);

        builderBand2 = new CuboidRoofedBuilder2();
        transformBand2 = new TransformToSurface(d, v, geomBand2);
        c2 = CuboidRoofed2.class;// SimpleCuboid2.class;
      }
    }

    // Cas où il n'y a qu'une seule bande, mais qu'on implante des bâtiments
    // derrière
    if (r2 == null && geomBand2 != null) {
      builderBand2 = new CuboidRoofedBuilder();
      transformBand2 = new TransformToSurface(d, v, geomBand2);
      c2 = CuboidRoofed.class;
    }

    // System.out.println("geomband1 " + geomBand1);
    // System.out.println("geomband2 " + geomBand2);

    ////////////////////////
    // Préparation des noyaux de modification
    ////////////////////////

    // Probabilité de s'implenter en bande 1 ou 2 (si c'est inférieur c'est
    // 2 et supérieur c'est 1)
    double p_simple = 0.5;

    Variate variate = new Variate(rng);
    // Probabilité de naissance-morts modifications
    List<Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>>> kernels = new ArrayList<>(
        3);
    NullView<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> nullView = new NullView<>();

    // Noyau pour la bande 1
    if (transformBand1 != null) {
      List<Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>>> lKernelsBand1 = new ArrayList<>();
      lKernelsBand1 = getBande1Kernels(variate, nullView, p, transformBand1,
          builderBand1, band1Parallel);
      kernels.addAll(lKernelsBand1);
    } else {
      p_simple = 1; // pas de transform on ne sera jamais dans la bande 1
    }

    System.out.println("kernel size before r2" + kernels.size());
    // Noyaux pour la bande 2
    if (r2 != null && transformBand2 != null) {
      List<Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>>> lKernelsBand2 = new ArrayList<>();
      lKernelsBand2 = getBande2Kernels(variate, nullView, p, transformBand2,
          builderBand2, band2parallel);
      kernels.addAll(lKernelsBand2);
    } else if (r2 == null && transformBand2 != null) { // Cas une seule
                                                       // bande et on bâtie
                                                       // derrière

      List<Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>>> lKernelsBand2 = new ArrayList<>();
      lKernelsBand2 = getBande1Kernels(variate, nullView, p, transformBand2,
          builderBand2, false);
      kernels.addAll(lKernelsBand2);

    } else {
      p_simple = 0; // pas de transform on ne sera jamais dans la bande 2
    }
    System.out.println("kernel size after r2" + kernels.size());

    // Si on ne peut pas construire dans la deuxième bande ni dans la
    // première ça sert à rien de continue
    if (kernels.isEmpty()) {
      return null;
    }

    // When direct sampling (solomon, etc.), what is the prob to choose a
    // simple cuboid

    // System.out.println("************* transformBand1 " + transformBand1);
    // System.out.println("************* transformBand2 " + transformBand2);
    // CuboidSampler objectSampler = new CuboidSampler(rng, p_simple,
    // transformSimple, transformParallel);
    MixCuboidRoofedSampler objectSampler = new MixCuboidRoofedSampler(rng,
        p_simple, transformBand1, transformBand2, builderBand1, builderBand2,
        c1, c2);

    // poisson distribution
    PoissonDistribution distribution = new PoissonDistribution(rng,
        p.getDouble("poisson"));
    DirectSampler<CuboidRoofed, GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> ds = new DirectRejectionSampler<>(
        distribution, objectSampler, pred);

    Acceptance<SimpleTemperature> acceptance = new MetropolisAcceptance<>();
    Sampler<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> s = new GreenSamplerBlockTemperature<>(
        rng, ds, acceptance, kernels);
    return s;
    // return null;
  }

  private static List<Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>>> getBande1Kernels(
      Variate variate,
      NullView<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> nullView,
      Parameters p, Transform transform, ObjectBuilder<CuboidRoofed> pbuilder,
      boolean parallel) {
    List<Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>>> kernels = new ArrayList<>();
    // Kernel de naissance
    // TODO Use a KernelProposalRatio to propose only birth when size is 0
    UniformTypeView<CuboidRoofed, GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> pView = null;

    if (parallel) {
      pView = new UniformTypeView<>(ParallelCuboidRoofed.class, pbuilder);

      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> kernel2 = new Kernel<>(
          nullView, pView, variate, variate, transform,
          p.getDouble("pbirthdeath"), p.getDouble("pbirth"), "Parallel");
      kernels.add(kernel2);

      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> parallelMovekernel = new Kernel<>(
          pView, pView, variate, variate,
          new MoveParallelRCuboid(p.getDouble("amplitudeMove")), 0.2, 1.0,
          "ParallelMoveP2");
      kernels.add(parallelMovekernel);

      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> parallelLength = new Kernel<>(
          pView, pView, variate, variate,
          new ChangeValue(p.getDouble("amplitudeMaxDim"), 6, 2), 0.2, 1.0,
          "ChgLengthP2");
      kernels.add(parallelLength);

      if ((p.getDouble("maxheightG") - p.getDouble("minheightG")) > 0.2) {
        double amplitudeHeight = p.getDouble("amplitudeHeight");
        Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> parallelHeight = new Kernel<>(
            pView, pView, variate, variate,
            new ChangeValue(amplitudeHeight, 6, 3), 0.2, 1.0, "ChgHeightP");
        kernels.add(parallelHeight);
      }
      if ((p.getDouble("maxheightT") - p.getDouble("minheightT")) > 0.2) {
        double amplitudeHeight = p.getDouble("amplitudeHeight");
        Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> parallelHeightT = new Kernel<>(
            pView, pView, variate, variate,
            new ChangeValue(amplitudeHeight, 6, 4), 0.2, 1.0, "ChgHeightTP");
        kernels.add(parallelHeightT);
      }
      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> changeDFSkernel = new Kernel<>(
          pView, pView, variate, variate,
          new ChangeValue(p.getDouble("amplitudeMaxDim"), 6, 5), 0.2, 1.0,
          "ChgDFS");
      kernels.add(changeDFSkernel);
    } else {
      pView = new UniformTypeView<>(CuboidRoofed.class, pbuilder); // faut il un
      // simpleCuboidRoofed
      // ?

      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> kernel1 = new Kernel<>(
          nullView, pView, variate, variate, transform,
          p.getDouble("pbirthdeath"), p.getDouble("pbirth"), "Simple");
      kernels.add(kernel1);

      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> simpleMovekernel = new Kernel<>(
          pView, pView, variate, variate,
          new MoveRCuboid(p.getDouble("amplitudeMove")), 0.2, 1.0,
          "SimpleMoveP2");
      kernels.add(simpleMovekernel);

      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> simpleLength = new Kernel<>(
          pView, pView, variate, variate,
          new ChangeValue(p.getDouble("amplitudeMaxDim"), 8, 2), 0.2, 1.0,
          "ChgLengthP2");
      kernels.add(simpleLength);
      if ((p.getDouble("maxheightG") - p.getDouble("minheightG")) > 0.2) {
        double amplitudeHeight = p.getDouble("amplitudeHeight");
        Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> simplelHeight = new Kernel<>(
            pView, pView, variate, variate,
            new ChangeValue(amplitudeHeight, 8, 4), 0.2, 1.0, "ChgHeightP");
        kernels.add(simplelHeight);
      }
      if ((p.getDouble("maxheightT") - p.getDouble("minheightT")) > 0.2) {
        double amplitudeHeight = p.getDouble("amplitudeHeight");
        Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> simpleHeightT = new Kernel<>(
            pView, pView, variate, variate,
            new ChangeValue(amplitudeHeight, 8, 6), 0.2, 1.0, "ChgHeightTP");
        kernels.add(simpleHeightT);
      }
      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> changeDFSkernel = new Kernel<>(
          pView, pView, variate, variate,
          new ChangeValue(p.getDouble("amplitudeMaxDim"), 8, 7), 0.2, 1.0,
          "ChgDFS");
      kernels.add(changeDFSkernel);
      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> simpleWidth = new Kernel<>(
          pView, pView, variate, variate,
          new ChangeValue(p.getDouble("amplitudeMaxDim"), 8, 3), 0.2, 1.0,
          "ChgWidth");
      kernels.add(simpleWidth);
      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> simpleRotatekernel = new Kernel<>(
          pView, pView, variate, variate,
          new RotateRCuboid(p.getDouble("amplitudeRotate") * Math.PI / 180),
          0.2, 1.0, "RotateS");
      kernels.add(simpleRotatekernel);
    }
    return kernels;
  }

  private static List<Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>>> getBande2Kernels(
      Variate variate,
      NullView<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> nullView,
      Parameters p, Transform transform, ObjectBuilder<CuboidRoofed> pbuilder,
      boolean parallel) {
    List<Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>>> kernels = new ArrayList<>();
    // Kernel de naissance
    // TODO Use a KernelProposalRatio to propose only birth when size is 0
    UniformTypeView<CuboidRoofed, GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> pView = null;

    if (parallel) {
      pView = new UniformTypeView<>(ParallelCuboidRoofed2.class, pbuilder);

      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> kernel2 = new Kernel<>(
          nullView, pView, variate, variate, transform,
          p.getDouble("pbirthdeath"), p.getDouble("pbirth"), "Parallel");
      kernels.add(kernel2);

      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> parallelMovekernel = new Kernel<>(
          pView, pView, variate, variate,
          new MoveParallelRCuboid(p.getDouble("amplitudeMove")), 0.2, 1.0,
          "ParallelMoveP2");
      kernels.add(parallelMovekernel);

      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> parallelLength = new Kernel<>(
          pView, pView, variate, variate,
          new ChangeValue(p.getDouble("amplitudeMaxDim"), 6, 2), 0.2, 1.0,
          "ChgLengthP2");
      kernels.add(parallelLength);

      if ((p.getDouble("maxheightG") - p.getDouble("minheightG")) > 0.2) {
        double amplitudeHeight = p.getDouble("amplitudeHeight");
        Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> parallelHeight = new Kernel<>(
            pView, pView, variate, variate,
            new ChangeValue(amplitudeHeight, 6, 3), 0.2, 1.0, "ChgHeightP");
        kernels.add(parallelHeight);
      }
      if ((p.getDouble("maxheightT") - p.getDouble("minheightT")) > 0.2) {
        double amplitudeHeight = p.getDouble("amplitudeHeight");
        Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> parallelHeightT = new Kernel<>(
            pView, pView, variate, variate,
            new ChangeValue(amplitudeHeight, 6, 4), 0.2, 1.0, "ChgHeightTP");
        kernels.add(parallelHeightT);
      }
      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> changeDFSkernel = new Kernel<>(
          pView, pView, variate, variate,
          new ChangeValue(p.getDouble("amplitudeMaxDim"), 6, 5), 0.2, 1.0,
          "ChgDFS");
      kernels.add(changeDFSkernel);
    } else {
      pView = new UniformTypeView<>(CuboidRoofed2.class, pbuilder); // faut il
                                                                    // un
      // simpleCuboidRoofed 1 et 2
      // ?

      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> kernel1 = new Kernel<>(
          nullView, pView, variate, variate, transform,
          p.getDouble("pbirthdeath"), p.getDouble("pbirth"), "Simple");
      kernels.add(kernel1);

      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> simpleMovekernel = new Kernel<>(
          pView, pView, variate, variate,
          new MoveRCuboid(p.getDouble("amplitudeMove")), 0.2, 1.0,
          "SimpleMoveP2");
      kernels.add(simpleMovekernel);

      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> simpleLength = new Kernel<>(
          pView, pView, variate, variate,
          new ChangeValue(p.getDouble("amplitudeMaxDim"), 8, 2), 0.2, 1.0,
          "ChgLengthP2");
      kernels.add(simpleLength);
      if ((p.getDouble("maxheightG") - p.getDouble("minheightG")) > 0.2) {
        double amplitudeHeight = p.getDouble("amplitudeHeight");
        Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> simplelHeight = new Kernel<>(
            pView, pView, variate, variate,
            new ChangeValue(amplitudeHeight, 8, 4), 0.2, 1.0, "ChgHeightP");
        kernels.add(simplelHeight);
      }
      if ((p.getDouble("maxheightT") - p.getDouble("minheightT")) > 0.2) {
        double amplitudeHeight = p.getDouble("amplitudeHeight");
        Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> simpleHeightT = new Kernel<>(
            pView, pView, variate, variate,
            new ChangeValue(amplitudeHeight, 8, 6), 0.2, 1.0, "ChgHeightTP");
        kernels.add(simpleHeightT);
      }
      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> changeDFSkernel = new Kernel<>(
          pView, pView, variate, variate,
          new ChangeValue(p.getDouble("amplitudeMaxDim"), 8, 7), 0.2, 1.0,
          "ChgDFS");
      kernels.add(changeDFSkernel);

      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> simpleWidth = new Kernel<>(
          pView, pView, variate, variate,
          new ChangeValue(p.getDouble("amplitudeMaxDim"), 8, 3), 0.2, 1.0,
          "ChgWidth");
      kernels.add(simpleWidth);
      Kernel<GraphConfiguration<CuboidRoofed>, BirthDeathModification<CuboidRoofed>> simpleRotatekernel = new Kernel<>(
          pView, pView, variate, variate,
          new RotateRCuboid(p.getDouble("amplitudeRotate") * Math.PI / 180),
          0.2, 1.0, "RotateS");
      kernels.add(simpleRotatekernel);
    }
    return kernels;
  }
}
