#pragma once

#include "sphinxsys.h"

using namespace SPH;

class MyocardiumPhysiology : public LocalMonoFieldElectroPhysiology {
 public:
  MyocardiumPhysiology(ElectroPhysiologyReaction* electro_physiology_reaction,
                       Real diff_cf,
                       Real bias_diff_cf,
                       Vec3d& bias_direction)
      : LocalMonoFieldElectroPhysiology(electro_physiology_reaction) {

    diff_cf_        = diff_cf;
    bias_diff_cf_   = bias_diff_cf;
    bias_direction_ = bias_direction;

    assignDerivedMaterialParameters();
    initializeDiffusion();
  }
};

class MyocardiumMuscle : public ActiveMuscle<LocallyOrthotropicMuscle> {
 public:
  MyocardiumMuscle(Real rho0, Real bulk_modulus, Vec3d& f0, Vec3d& s0, Real* a_0, Real* b_0)
      : ActiveMuscle<LocallyOrthotropicMuscle>() {

    rho0_         = rho0;
    bulk_modulus_ = bulk_modulus;
    f0_           = f0;
    s0_           = s0;
    std::copy(a_0, a_0 + 4, a0_);
    std::copy(b_0, b_0 + 4, b0_);

    assignDerivedMaterialParameters();
  }
};

class MuscleReactionModel : public AlievPanfilowModel {
 public:
  MuscleReactionModel(Real k_a, Real c_m, Real k, Real a, Real b, Real mu_1, Real mu_2, Real epsilon)
      : AlievPanfilowModel() {

    k_a_     = k_a;
    c_m_     = c_m;
    k_       = k;
    a_       = a;
    b_       = b;
    mu_1_    = mu_1;
    mu_2_    = mu_2;
    epsilon_ = epsilon;

    assignDerivedReactionParameters();
  }
};

class HeartBody : public SolidBody {
 public:
  HeartBody(SPHSystem& system, string body_name, std::string& path_to_file, Real length_scale)
      : SolidBody(system, body_name) {

    ComplexShape original_body_shape;
    original_body_shape.addTriangleMeshShape(CreateHeart(path_to_file, length_scale),
                                             ShapeBooleanOps::add);
    body_shape_ = new LevelSetComplexShape(this, original_body_shape);
  }

 private:
  TriangleMeshShape* CreateHeart(std::string& path_to_file, Real length_scale) {

    Vecd translation(-53.5 * length_scale, -70.0 * length_scale, -32.5 * length_scale);
    TriangleMeshShape* geometry_myocardium =
        new TriangleMeshShape(path_to_file, translation, length_scale);
    return geometry_myocardium;
  }
};

/// Setup diffusion material properties for mapping the fiber direction
class DiffusionMaterial
    : public DiffusionReactionMaterial<ElasticSolidParticles, LocallyOrthotropicMuscle> {

 public:
  DiffusionMaterial()
      : DiffusionReactionMaterial<ElasticSolidParticles, LocallyOrthotropicMuscle>() {

    insertASpecies("Phi");
    assignDerivedMaterialParameters();
    initializeDiffusion();
  }

  virtual void initializeDiffusion() override {

    IsotropicDiffusion* phi_diffusion =
        new IsotropicDiffusion(species_indexes_map_["Phi"], species_indexes_map_["Phi"]);
    species_diffusion_.push_back(phi_diffusion);
  };
};

/// Set diffusion relaxation
class DiffusionRelaxation : public RelaxationOfAllDiffusionSpeciesRK2<
                                SolidBody,
                                ElasticSolidParticles,
                                LocallyOrthotropicMuscle,
                                RelaxationOfAllDiffussionSpeciesInner<SolidBody,
                                                                      ElasticSolidParticles,
                                                                      LocallyOrthotropicMuscle>,
                                BodyRelationInner> {
 public:
  DiffusionRelaxation(BodyRelationInner* body_inner_relation)
      : RelaxationOfAllDiffusionSpeciesRK2(body_inner_relation){};

  virtual ~DiffusionRelaxation(){};
};

/// Impose diffusion boundary conditions
class DiffusionBCs : public ConstrainDiffusionBodyRegion<SolidBody,
                                                         ElasticSolidParticles,
                                                         ShapeSurface,
                                                         LocallyOrthotropicMuscle> {

 protected:
  size_t phi_;

  virtual void Update(size_t index_i, Real dt = 0.0) override {

    Vecd dist_2_face = body_->body_shape_->findNormalDirection(pos_n_[index_i]);
    Vecd face_norm   = dist_2_face / (dist_2_face.norm() + 1.0e-15);
    Vecd center_norm = pos_n_[index_i] / (pos_n_[index_i].norm() + 1.0e-15);
    Real angle       = dot(face_norm, center_norm);
    if (angle >= 0.0) {
      species_n_[phi_][index_i] = 1.0;
    } else {
      if (pos_n_[index_i][1] < -body_->particle_adaptation_->ReferenceSpacing())
        species_n_[phi_][index_i] = 0.0;
    }
  };

 public:
  DiffusionBCs(SolidBody* body, ShapeSurface* body_part)
      : ConstrainDiffusionBodyRegion<SolidBody,
                                     ElasticSolidParticles,
                                     ShapeSurface,
                                     LocallyOrthotropicMuscle>(body, body_part) {

    phi_ = material_->SpeciesIndexMap()["Phi"];
  };

  virtual ~DiffusionBCs(){};
};

/// Update fiber and sheet direction after diffusion
class ComputeFiberandSheetDirections
    : public DiffusionBasedMapping<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle> {
 protected:
  size_t phi_;
  Real beta_epi_, beta_endo_;
  /// The centerline vector that is parallel to the ventricular centerline and pointing apex-to-base
  Vecd center_line_;

  virtual void Update(size_t index_i, Real dt = 0.0) override {

    /// Probe the face norm from levelset field.
    Vecd dist_2_face = body_->body_shape_->findNormalDirection(pos_n_[index_i]);
    Vecd face_norm   = dist_2_face / (dist_2_face.norm() + 1.0e-15);
    Vecd center_norm = pos_n_[index_i] / (pos_n_[index_i].norm() + 1.0e-15);
    if (dot(face_norm, center_norm) <= 0.0) {
      face_norm = -face_norm;
    }
    /// Compute the centerline's projection on the plane orthogonal to face norm
    Vecd circumferential_direction = SimTK::cross(center_line_, face_norm);
    Vecd cd_norm = circumferential_direction / (circumferential_direction.norm() + 1.0e-15);
    /// The rotation angle is given by beta = (beta_epi - beta_endo) phi + beta_endo
    Real beta = (beta_epi_ - beta_endo_) * species_n_[phi_][index_i] + beta_endo_;
    /// Compute the rotation matrix through Rodrigues rotation formulation
    Vecd f_0 = cos(beta) * cd_norm + sin(beta) * SimTK::cross(face_norm, cd_norm) +
               dot(face_norm, cd_norm) * (1.0 - cos(beta)) * face_norm;

    if (pos_n_[index_i][1] < -body_->particle_adaptation_->ReferenceSpacing()) {
      material_->local_f0_[index_i] = f_0 / (f_0.norm() + 1.0e-15);
      material_->local_s0_[index_i] = face_norm;
    } else {
      material_->local_f0_[index_i] = Vecd(0);
      material_->local_s0_[index_i] = Vecd(0);
    }
  };

 public:
  ComputeFiberandSheetDirections(SolidBody* body)
      : DiffusionBasedMapping<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle>(body) {

    phi_         = material_->SpeciesIndexMap()["Phi"];
    center_line_ = Vecd(0.0, 1.0, 0.0);
    beta_epi_    = -(70.0 / 180.0) * M_PI;
    beta_endo_   = (80.0 / 180.0) * M_PI;
  };

  virtual ~ComputeFiberandSheetDirections(){};
};

/// Define the beam base which will be constrained.
/// Must be instantiated after all body particles have been generated.
class MuscleBase : public BodyPartByParticle {
 public:
  MuscleBase(SolidBody* solid_body,
             std::string constrained_region_name,
             BoundingBox& domain_bounds,
             Real dp_0,
             Real length_scale)
      : BodyPartByParticle(solid_body, constrained_region_name) {

    body_part_shape_ = new ComplexShape(constrained_region_name);
    body_part_shape_->addTriangleMeshShape(CreateBaseShape(domain_bounds, dp_0, length_scale),
                                           ShapeBooleanOps::add);

    /// Tag the constrained particles to the base for constraint
    tagBodyPart();
  }

 private:
  TriangleMeshShape* CreateBaseShape(BoundingBox& domain_bounds, Real dp_0, Real length_scale) {

    Real l = domain_bounds.second[0] - domain_bounds.first[0];
    Real w = domain_bounds.second[2] - domain_bounds.first[2];
    Vecd halfsize_shape(0.5 * l, 1.0 * dp_0, 0.5 * w);
    Vecd translation_shape(-10.0 * length_scale, -1.0 * dp_0, 0.0);
    TriangleMeshShape* geometry = new TriangleMeshShape(halfsize_shape, 20, translation_shape);

    return geometry;
  }
};