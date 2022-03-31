#pragma once

#include "sphinxsys.h"
#include <unordered_map>

#define DB_PERLIN_IMPL
#include "db_perlin.hpp"

using namespace SPH;

class AnisotropicAlievPanfilowModel : public AlievPanfilowModel
{
protected:
  unordered_map<size_t, Real> heterogeneities_;

  virtual Real getProductionRateIonicCurrent(StdVec<StdLargeVec<Real>> &species, size_t particle_i) override
  {
    Real voltage = species[voltage_][particle_i];
    return -k_ * voltage * (voltage * voltage - a_ * voltage - voltage) / c_m_;
  }
  virtual Real getLossRateIonicCurrent(StdVec<StdLargeVec<Real>> &species, size_t particle_i) override
  {
    Real gate_variable = species[gate_variable_][particle_i];
    return (k_ * a_ + gate_variable) / c_m_;
  }
  virtual Real getProductionRateGateVariable(StdVec<StdLargeVec<Real>> &species, size_t particle_i) override
  {
    Real voltage = species[voltage_][particle_i];
    Real gate_variable = species[gate_variable_][particle_i];
    Real temp = epsilon_ + mu_1_ * gate_variable / (mu_2_ + voltage + Eps);
    return -temp * k_ * voltage * (voltage - (b_ * heterogeneities_[particle_i]) - 1.0);
  }
  virtual Real getLossRateGateVariable(StdVec<StdLargeVec<Real>> &species, size_t particle_i) override
  {
    Real voltage = species[voltage_][particle_i];
    Real gate_variable = species[gate_variable_][particle_i];
    return epsilon_ + mu_1_ * gate_variable / (mu_2_ + voltage + Eps);
  }

public:
  explicit AnisotropicAlievPanfilowModel(Real k_a, Real c_m, Real k, Real a, Real b, Real mu_1, Real mu_2, Real epsilon)
      : AlievPanfilowModel(k_a, c_m, k, a, b, mu_1, mu_2, epsilon){};
  virtual ~AnisotropicAlievPanfilowModel(){};

  void generate_heterogeneities(ElectroPhysiologyParticles &particles, Real chaos, Real scaling)
  {
    for (size_t i = 0; i < particles.pos_0_.size(); i++)
    {
      Vec3d pos = particles.pos_0_[i];
      heterogeneities_[i] = (Real)1.0 + scaling * db::perlin(pos[0]/chaos, pos[1]/chaos, pos[2]/chaos);
    }
  }
};

/** Define geometry and initial conditions of SPH bodies. */
class HeartBody : public SolidBody
{
public:
  HeartBody(SPHSystem &system, const std::string &body_name, const std::string &path_to_file, Real length_scale)
      : SolidBody(system, body_name)
  {
    Vecd translation(-53.5 * length_scale, -70.0 * length_scale, -32.5 * length_scale);
    TriangleMeshShapeSTL triangle_mesh_heart_shape(path_to_file, translation, length_scale);
    body_shape_.add<LevelSetShape>(this, triangle_mesh_heart_shape);
  }
};

/** Set diffusion relaxation method. */
class DiffusionRelaxation
    : public RelaxationOfAllDiffusionSpeciesRK2<
          SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle,
          RelaxationOfAllDiffussionSpeciesInner<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle>,
          BodyRelationInner>
{
public:
  explicit DiffusionRelaxation(BodyRelationInner &body_inner_relation)
      : RelaxationOfAllDiffusionSpeciesRK2(body_inner_relation){};
  virtual ~DiffusionRelaxation(){};
};

/** Imposing diffusion boundary condition */
class DiffusionBCs
    : public ConstrainDiffusionBodyRegion<SolidBody, ElasticSolidParticles, BodySurface, LocallyOrthotropicMuscle>
{
protected:
  size_t phi_;
  virtual void Update(size_t index_i, Real dt = 0.0) override
  {
    Vecd dist_2_face = body_->body_shape_.findNormalDirection(pos_n_[index_i]);
    Vecd face_norm = dist_2_face / (dist_2_face.norm() + 1.0e-15);

    Vecd center_norm = pos_n_[index_i] / (pos_n_[index_i].norm() + 1.0e-15);

    Real angle = dot(face_norm, center_norm);
    if (angle >= 0.0)
    {
      species_n_[phi_][index_i] = 1.0;
    }
    else
    {
      if (pos_n_[index_i][1] < -body_->sph_adaptation_->ReferenceSpacing())
        species_n_[phi_][index_i] = 0.0;
    }
  };

public:
  DiffusionBCs(SolidBody &body, BodySurface &body_part)
      : ConstrainDiffusionBodyRegion<SolidBody, ElasticSolidParticles, BodySurface, LocallyOrthotropicMuscle>(body, body_part)
  {
    phi_ = material_->SpeciesIndexMap()["Phi"];
  };
  virtual ~DiffusionBCs(){};
};

/** Compute Fiber and Sheet direction after diffusion */
class ComputeFiberandSheetDirections
    : public DiffusionBasedMapping<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle>
{
protected:
  size_t phi_;
  Real beta_epi_, beta_endo_;
  /** We define the centerline vector, which is parallel to the ventricular centerline and pointing  apex-to-base.*/
  Vecd center_line_;
  virtual void Update(size_t index_i, Real dt = 0.0) override
  {
    /**
     * Ref: original doi.org/10.1016/j.euromechsol.2013.10.009
     * 		Present  doi.org/10.1016/j.cma.2016.05.031
     */
    /** Probe the face norm from Levelset field. */
    Vecd dist_2_face = body_->body_shape_.findNormalDirection(pos_n_[index_i]);
    Vecd face_norm = dist_2_face / (dist_2_face.norm() + 1.0e-15);
    Vecd center_norm = pos_n_[index_i] / (pos_n_[index_i].norm() + 1.0e-15);
    if (dot(face_norm, center_norm) <= 0.0)
    {
      face_norm = -face_norm;
    }
    /** Compute the centerline's projection on the plane orthogonal to face norm. */
    Vecd circumferential_direction = SimTK::cross(center_line_, face_norm);
    Vecd cd_norm = circumferential_direction / (circumferential_direction.norm() + 1.0e-15);
    /** The rotation angle is given by beta = (beta_epi - beta_endo) phi + beta_endo */
    Real beta = (beta_epi_ - beta_endo_) * species_n_[phi_][index_i] + beta_endo_;
    /** Compute the rotation matrix through Rodrigues rotation formulation. */
    Vecd f_0 = cos(beta) * cd_norm + sin(beta) * SimTK::cross(face_norm, cd_norm) +
               dot(face_norm, cd_norm) * (1.0 - cos(beta)) * face_norm;

    if (pos_n_[index_i][1] < -body_->sph_adaptation_->ReferenceSpacing())
    {
      material_->local_f0_[index_i] = f_0 / (f_0.norm() + 1.0e-15);
      material_->local_s0_[index_i] = face_norm;
    }
    else
    {
      material_->local_f0_[index_i] = Vecd(0);
      material_->local_s0_[index_i] = Vecd(0);
    }
  };

public:
  explicit ComputeFiberandSheetDirections(SolidBody &body)
      : DiffusionBasedMapping<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle>(body)
  {
    phi_ = material_->SpeciesIndexMap()["Phi"];
    center_line_ = Vecd(0.0, 1.0, 0.0);
    beta_epi_ = -(70.0 / 180.0) * M_PI;
    beta_endo_ = (80.0 / 180.0) * M_PI;
  };
  virtual ~ComputeFiberandSheetDirections(){};
};

//	define shape parameters which will be used for the constrained body part.
class MuscleBaseShapeParameters : public TriangleMeshShapeBrick::ShapeParameters
{
public:
  MuscleBaseShapeParameters(const Vec3d &domain_upper_bound, const Vec3d &domain_lower_bound, Real length_scale, Real dp_0) : TriangleMeshShapeBrick::ShapeParameters()
  {
    Real l = domain_upper_bound[0] - domain_lower_bound[0];
    Real w = domain_upper_bound[2] - domain_lower_bound[2];
    halfsize_ = Vec3d(0.5 * l, 1.0 * dp_0, 0.5 * w);
    resolution_ = 20;
    translation_ = Vec3d(-10.0 * length_scale, -1.0 * dp_0, 0.0);
  }
};