/**
 * 			Unit :
 *			time t = ms = 12.9 [-]
 * 			length l = mm
 * 			mass m = g
 *			density rho = g * (mm)^(-3)
 *			Pressure pa = g * (mm)^(-1) * (ms)^(-2)
 *			diffusion d = (mm)^(2) * (ms)^(-2)
 */

#include "CLI11.hpp"
#include "binary_writer.h"
#include "sphinxsys.h"
#include "stimulus.h"
#include "models.h"

using namespace SPH;

int main(int argc, char *argv[])
{
  // ===========================================================================
  /// Parameters
  // ===========================================================================

  /// Full path to stl model file.
  std::string full_path_to_stl_file = "./input/heart-new.stl";

  /// Scaling factors
  Real length_scale = 1.0;
  Real time_scale = 1.0 / 12.9;
  Real stress_scale = 1.0e-6;

  /// System domain bounds
  Vec3d domain_lower_bound(-55.0 * length_scale, -75.0 * length_scale, -35.0 * length_scale);
  Vec3d domain_upper_bound(35.0 * length_scale, 5.0 * length_scale, 35.0 * length_scale);
  BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);

  /// Initial particle spacing
  Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 45.0;

  /// Density of myocardium muscle
  Real rho0_s = 1.06e-3;

  /// Fiber (f^0) and sheet (s^0) directions.
  Vec3d fiber_direction(1.0, 0.0, 0.0);
  Vec3d sheet_direction(0.0, 1.0, 0.0);

  /// Capacitance of cell membrane.
  Real c_m = 1.0;

  /// Conductivity parameters (d^iso and d^ani ?).
  Real diffusion_coff = 1.2;
  Real bias_coff = 0.0;

  /// Aliev-Panfilov parameters.
  Real k = 8.0;
  Real a = 0.01;
  Real b = 0.15;
  Real epsilon = 0.002;
  Real mu_1 = 0.2;
  Real mu_2 = 0.3;
  Real chaos = 35.0;
  Real scaling = 1.0;

  /// Passive stress parameters.
  Real a0[4] = {496.0 * stress_scale, 15196.0 * stress_scale, 3283.0 * stress_scale,
                662.0 * stress_scale};
  Real b0[4] = {7.209, 20.417, 11.176, 9.466};

  // Compressibility parameters (weakly compressible).
  Real poisson = 0.4995;
  Real bulk_modulus = 2.0 * a0[0] * (1.0 + poisson) / (3.0 * (1.0 - 2.0 * poisson));

  /// Active stress parameter.
  Real k_a = 100 * stress_scale;

  /// Simulation parameters
  bool run_particle_relaxation = false;
  bool reload_particles = true;
  bool add_observers = false;

  Real end_time = 300;

  Real stimulus_radius = 10;
  int extra_stimulus_count = 1;

  int relax_steps = 1000;
  int relax_diffusion_steps = 100;

  /// Data collection
  int run_id = 0;
  int sampling_interval = 1;
  int time_step = 1;

  // ===========================================================================
  /// CLI
  // ===========================================================================

  CLI::App app{"Computes SPH simulations of a 3D heart"};

  app.set_help_all_flag("--help-all", "Expand all help");
  app.set_config("--config", "sph_heart.toml", "config file");

  app.add_flag("-r,--relax", run_particle_relaxation, "Particle relaxation.");
  app.add_flag("--reload", reload_particles, "Reload particles from previous simulation.");
  app.add_flag("--add_observers", add_observers);

  app.add_option("--stl", full_path_to_stl_file, "Heart STL file")->group("Simulation Parameters");
  app.add_option("--end_time", end_time, "End time of simulation.")->group("Simulation Parameters");

  app.add_option("--radius", stimulus_radius)->group("Stimulus Parameters");
  app.add_option("--count", extra_stimulus_count)->group("Stimulus Parameters");

  app.add_option("--iso", diffusion_coff, "Isotrophic conductivity parameter")
      ->group("Conductivity Parameters");
  app.add_option("--aniso", bias_coff, "Anisotropic conductivity parameter")
      ->group("Conductivity Parameters");

  app.add_option("--ep_k", k)->group("Electrophysiology Parameters");
  app.add_option("--ep_a", a)->group("Electrophysiology Parameters");
  app.add_option("--ep_b", b)->group("Electrophysiology Parameters");
  app.add_option("--ep_mu_1", mu_1)->group("Electrophysiology Parameters");
  app.add_option("--ep_mu_2", mu_2)->group("Electrophysiology Parameters");
  app.add_option("--ep_epsilon", epsilon)->group("Electrophysiology Parameters");
  app.add_option("--chaos", chaos)->group("Electrophysiology Parameters");
  app.add_option("--chaos_scaling", scaling)->group("Electrophysiology Parameters");

  app.add_option("--relax_steps", relax_steps, "Relaxation steps")->group("Particle Relaxation");
  app.add_option("--relax_diffusion_steps", relax_diffusion_steps, "Relaxation Diffusion steps")
      ->group("Particle Relaxation");

  app.add_option("--id", run_id, "Run ID")->group("Data Collection");
  app.add_option("--interval", sampling_interval, "Interval between output samples")
      ->group("Data Collection");
  app.add_option("--step", time_step, "Time step")->group("Data Collection");

  try
  {
    app.parse(argc, argv);
  }
  catch (const CLI::ParseError &e)
  {
    return app.exit(e);
  }

  std::cout << app.config_to_str(false, true) << std::endl;
  std::cout << "################ END OF CONFIG ################\n\n"
            << std::endl;

  // ===========================================================================
  /// SPH system
  // ===========================================================================

  SPHSystem system(system_domain_bounds, dp_0);
  GlobalStaticVariables::physical_time_ = 0.0;
  system.run_particle_relaxation_ = run_particle_relaxation;
  system.reload_particles_ = reload_particles;
  /// Tag for computation from restart files. 0: not from restart files
  system.restart_step_ = 0;

  In_Output in_output(system);
  if (!run_particle_relaxation && reload_particles && !fs::exists(in_output.reload_folder_))
  {
    cout
        << "ERROR: Reload direcory '" << in_output.reload_folder_
        << "' does not exist. First run the program with '--relax' to perform particle relaxation.\n"
        << endl;
    return app.exit(CLI::Error("CallForHelp", "", CLI::ExitCodes::FileError));
  }

  // ===========================================================================
  /// SPH data
  // ===========================================================================

  /// Create SPH physiology body, material, and particles
  HeartBody physiology_body(system, "ExcitationHeart", full_path_to_stl_file, length_scale);
  SharedPtr<ParticleGenerator> physiology_body_particle_generator = makeShared<ParticleGeneratorLattice>();
  if (!system.run_particle_relaxation_ && system.reload_particles_)
    physiology_body_particle_generator = makeShared<ParticleGeneratorReload>(in_output, physiology_body.getBodyName());

  AnisotropicAlievPanfilowModel muscle_reaction_model(k_a, c_m, k, a, b, mu_1, mu_2, epsilon);

  SharedPtr<LocalMonoFieldElectroPhysiology> myocardium_excitation =
      makeShared<LocalMonoFieldElectroPhysiology>(muscle_reaction_model, diffusion_coff, bias_coff, fiber_direction);
  ElectroPhysiologyParticles physiology_particles(physiology_body, myocardium_excitation, physiology_body_particle_generator);

  muscle_reaction_model.generate_heterogeneities(physiology_particles, chaos, scaling);

  /// Create SPH mechanical body, material, and particles
  HeartBody mechanics_body(system, "ContractionHeart", full_path_to_stl_file, length_scale);
  SharedPtr<ParticleGenerator> mechanics_body_particle_generator = makeShared<ParticleGeneratorLattice>();
  if (!system.run_particle_relaxation_ && system.reload_particles_)
    mechanics_body_particle_generator = makeShared<ParticleGeneratorReload>(in_output, mechanics_body.getBodyName());
  SharedPtr<ActiveMuscle<LocallyOrthotropicMuscle>> myocardium_muscle =
      makeShared<ActiveMuscle<LocallyOrthotropicMuscle>>(rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0);
  ActiveMuscleParticles mechanics_particles(mechanics_body, myocardium_muscle, mechanics_body_particle_generator);

  /** check whether reload material properties. */
  if (!system.run_particle_relaxation_ && system.reload_particles_)
  {
    ReloadMaterialParameterIO read_muscle_fiber_and_sheet(in_output, myocardium_muscle);
    ReloadMaterialParameterIO read_myocardium_excitation_fiber(in_output, myocardium_excitation, myocardium_muscle->LocalParametersName());
    read_muscle_fiber_and_sheet.readFromFile();
    read_myocardium_excitation_fiber.readFromFile();
  }

  // ===========================================================================
  /// Particle relaxation
  // ===========================================================================

  if (system.run_particle_relaxation_)
  {
    HeartBody relax_body(system, "RelaxationHeart", full_path_to_stl_file, length_scale);
    StdVec<std::string> species_name_list{"Phi"};
    SharedPtr<DiffusionReaction<ElasticSolidParticles, LocallyOrthotropicMuscle>> relax_body_material =
        makeShared<DiffusionReaction<ElasticSolidParticles, LocallyOrthotropicMuscle>>(
            species_name_list, rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0);
    relax_body_material->initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi", diffusion_coff);
    DiffusionReactionParticles<ElasticSolidParticles, LocallyOrthotropicMuscle> diffusion_particles(relax_body, relax_body_material);

    /** topology */
    BodyRelationInner relax_body_inner(relax_body);
    /** Random reset the relax solid particle position. */
    RandomizePartilePosition random_particles(relax_body);

    /**
     * @brief 	Algorithms for particle relaxation.
     */
    /** A  Physics relaxation step. */
    relax_dynamics::RelaxationStepInner relaxation_step_inner(relax_body_inner);
    /**
     * Diffusion process.
     */
    /** Time step for diffusion. */
    GetDiffusionTimeStepSize<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle> get_time_step_size(relax_body);
    /** Diffusion process for diffusion body. */
    DiffusionRelaxation diffusion_relaxation(relax_body_inner);
    /** Compute the fiber and sheet after diffusion. */
    ComputeFiberandSheetDirections compute_fiber_sheet(relax_body);
    /** Write the particle reload files. */
    ReloadParticleIO write_particle_reload_files(in_output, {&relax_body, &relax_body}, {physiology_body.getBodyName(), mechanics_body.getBodyName()});
    /** Write material property to xml file. */
    ReloadMaterialParameterIO write_material_property(in_output, relax_body_material, myocardium_muscle->LocalParametersName());
    /**
     * @brief 	Physics relaxation starts here.
     */
    /** Relax the elastic structure. */
    random_particles.parallel_exec(0.25);
    relaxation_step_inner.surface_bounding_.parallel_exec();
    /**
     * From here the time stepping begins.
     * Set the starting time.
     */
    int ite = 0;
    int relax_step = 1000;
    int diffusion_step = 100;
    while (ite < relax_step)
    {
      relaxation_step_inner.parallel_exec();
      ite++;
      if (ite % 100 == 0)
      {
        std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
      }
    }

    BodySurface surface_part(relax_body);
    /** constraint boundary condition for diffusion. */
    DiffusionBCs impose_diffusion_bc(relax_body, surface_part);
    impose_diffusion_bc.parallel_exec();

    Real dt = get_time_step_size.parallel_exec();
    while (ite <= diffusion_step + relax_step)
    {
      diffusion_relaxation.parallel_exec(dt);
      impose_diffusion_bc.parallel_exec();
      if (ite % 10 == 0)
      {
        std::cout << "Diffusion steps N=" << ite - relax_step << "	dt: " << dt << "\n";
      }
      ite++;
    }
    compute_fiber_sheet.exec();
    ite++;
    compute_fiber_sheet.parallel_exec();
    write_material_property.writeToFile(0);
    write_particle_reload_files.writeToFile(0);

    fs::remove_all(in_output.output_folder_);
    return 0;
  }

  // ===========================================================================
  /// Body relation/topology
  // ===========================================================================

  BodyRelationInner physiology_body_inner(physiology_body);
  BodyRelationInner mechanics_body_inner(mechanics_body);
  BodyRelationContact physiology_body_contact(physiology_body, {&mechanics_body});
  BodyRelationContact mechanics_body_contact(mechanics_body, {&physiology_body});

  // ===========================================================================
  /// SPH methods
  // ===========================================================================

  // Corrected configuration.
  solid_dynamics::CorrectConfiguration correct_configuration_excitation(physiology_body_inner);
  // Time step size calculation.
  electro_physiology::GetElectroPhysiologyTimeStepSize get_physiology_time_step(physiology_body);
  // Diffusion process for diffusion body.
  electro_physiology::ElectroPhysiologyDiffusionRelaxationInner diffusion_relaxation(physiology_body_inner);
  // Solvers for ODE system.
  electro_physiology::ElectroPhysiologyReactionRelaxationForward reaction_relaxation_forward(physiology_body);
  electro_physiology::ElectroPhysiologyReactionRelaxationBackward reaction_relaxation_backward(physiology_body);
  // Active mechanics.
  solid_dynamics::CorrectConfiguration correct_configuration_contraction(mechanics_body_inner);
  observer_dynamics::CorrectInterpolationKernelWeights correct_kernel_weights_for_interpolation(mechanics_body_contact);
  /** Interpolate the active contract stress from electrophysiology body. */
  observer_dynamics::InterpolatingAQuantity<Real> active_stress_interpolation(mechanics_body_contact, "ActiveContractionStress", "ActiveContractionStress");
  /** Interpolate the particle position in physiology_body  from mechanics_body. */
  observer_dynamics::InterpolatingAQuantity<Vecd> interpolation_particle_position(physiology_body_contact, "Position", "Position");
  /** Time step size calculation. */
  solid_dynamics::AcousticTimeStepSize get_mechanics_time_step(mechanics_body);
  /** active and passive stress relaxation. */
  solid_dynamics::StressRelaxationFirstHalf stress_relaxation_first_half(mechanics_body_inner);
  solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(mechanics_body_inner);
  /** Constrain region of the inserted body. */
  MuscleBaseShapeParameters muscle_base_parameters(domain_upper_bound, domain_lower_bound, length_scale, dp_0);
  TriangleMeshShapeBrick muscle_base_shape(muscle_base_parameters);
  BodyRegionByParticle muscle_base(mechanics_body, "Holder", muscle_base_shape);
  solid_dynamics::ConstrainSolidBodyRegion constrain_holder(mechanics_body, muscle_base);

  // ===========================================================================
  /// Random Engine
  // ===========================================================================

  std::random_device r;
  std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
  std::mt19937 rng(seed);
  std::uniform_int_distribution<size_t> random_particle_index(
      0, physiology_body.base_particles_->pos_n_.size() - 1);
  std::uniform_real_distribution<Real> random_start_time(0, end_time);

  // ===========================================================================
  /// Stimulation
  // ===========================================================================

  RandomPulseTrain rpt(&physiology_body, random_particle_index(rng), stimulus_radius);
  for (int i = 0; i < extra_stimulus_count; i++)
  {
    rpt.register_extra_stimulus(random_particle_index(rng));
  }

  // ===========================================================================
  /// Output
  // ===========================================================================

  BodyStatesRecordingToNpz write_states(in_output, system.real_bodies_, run_id);

  // ===========================================================================
  /// Pre-simulation
  // ===========================================================================

  system.initializeSystemCellLinkedLists();
  system.initializeSystemConfigurations();
  correct_configuration_excitation.parallel_exec();
  correct_configuration_contraction.parallel_exec();
  correct_kernel_weights_for_interpolation.parallel_exec();

  // ===========================================================================
  /// Main loop
  // ===========================================================================

  int screen_output_interval = 10;
  int ite = 0;
  int reaction_step = 2;
  Real Ouput_T = end_time / 200.0;
  Real Observer_time = 0.01 * Ouput_T;
  Real dt = 0.0;   /**< Default acoustic time step sizes for physiology. */
  Real dt_s = 0.0; /**< Default acoustic time step sizes for mechanics. */
  /// Statistics for computing time
  tick_count t1 = tick_count::now();
  tick_count::interval_t interval;
  std::cout << "Main Loop Starts Here : " << std::endl;
  /// Track sample recording status
  int write_tracker = 0;
  int iter_cnt = 0;
  /** Main loop starts here. */
  while (GlobalStaticVariables::physical_time_ < end_time)
  {
    Real integration_time = 0.0;
    while (integration_time < Ouput_T)
    {
      Real relaxation_time = 0.0;
      while (relaxation_time < Observer_time)
      {
        if (ite % screen_output_interval == 0)
        {
          std::cout << std::fixed << std::setprecision(9) << "N=" << ite << "	Time = "
                    << GlobalStaticVariables::physical_time_
                    << "	dt = " << dt
                    << "	dt_s = " << dt_s << "\n";
        }

        /// Apply stimulation
        rpt.generate_and_execute(dt);

        /**Strong splitting method. */
        // forward reaction
        int ite_forward = 0;
        while (ite_forward < reaction_step)
        {
          reaction_relaxation_forward.parallel_exec(0.5 * dt / Real(reaction_step));
          ite_forward++;
        }
        /** 2nd Runge-Kutta scheme for diffusion. */
        diffusion_relaxation.parallel_exec(dt);

        // backward reaction
        int ite_backward = 0;
        while (ite_backward < reaction_step)
        {
          reaction_relaxation_backward.parallel_exec(0.5 * dt / Real(reaction_step));
          ite_backward++;
        }

        active_stress_interpolation.parallel_exec();

        Real dt_s_sum = 0.0;
        while (dt_s_sum < dt)
        {
          dt_s = get_mechanics_time_step.parallel_exec();
          if (dt - dt_s_sum < dt_s)
            dt_s = dt - dt_s_sum;
          stress_relaxation_first_half.parallel_exec(dt_s);
          constrain_holder.parallel_exec(dt_s);
          stress_relaxation_second_half.parallel_exec(dt_s);
          dt_s_sum += dt_s;
        }

        ite++;
        dt = get_physiology_time_step.parallel_exec();

        relaxation_time += dt;
        integration_time += dt;
        GlobalStaticVariables::physical_time_ += dt;
      }
    }
    tick_count t2 = tick_count::now();
    interpolation_particle_position.parallel_exec();
    if (write_tracker % sampling_interval <= time_step)
    {
      write_states.writeToFile(iter_cnt);
    }
    write_tracker++;
    iter_cnt++;
    tick_count t3 = tick_count::now();
    interval += t3 - t2;
  }
  tick_count t4 = tick_count::now();

  tick_count::interval_t tt;
  tt = t4 - t1 - interval;
  std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
  return 0;
}

/// @}