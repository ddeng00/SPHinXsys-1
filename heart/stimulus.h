#pragma once

#include "sphinxsys.h"
#include <unordered_map>

using namespace SPH;

class StimulusCurrent : public electro_physiology::ElectroPhysiologyInitialCondition {
 protected:
  size_t voltage_;
  size_t gate_variable_;

 public:
  StimulusCurrent(SolidBody* muscle)
      : electro_physiology::ElectroPhysiologyInitialCondition(*muscle) {

    voltage_       = material_->SpeciesIndexMap()["Voltage"];
    gate_variable_ = material_->SpeciesIndexMap()["GateVariable"];
  }
};

class IndexSphereStimulus : public StimulusCurrent {
 protected:
  const size_t center_particle_index_i_;
  const Real radius_;
  const Real voltage_target_;

  void Update(size_t index_i, Real dt) override {

    Real distance = (pos_n_[index_i] - pos_n_[center_particle_index_i_]).norm();
    if (distance < radius_) {
      species_n_[voltage_][index_i] = voltage_target_;
    }
  };

 public:
  IndexSphereStimulus(SolidBody* muscle,
                      size_t center_particle_index_i,
                      Real radius,
                      Real voltage_target = 0.92)
      : StimulusCurrent(muscle),
        center_particle_index_i_(center_particle_index_i),
        radius_(radius),
        voltage_target_(voltage_target) {}
};

class RandomPulseTrain {
 private:
  const size_t s1_index_;
  list<size_t> extra_stimuli_indexes_;
  unordered_map<size_t, bool> ready_map_;

  SolidBody* muscle_;
  SPH::ElectroPhysiologyParticles* particles_;
  size_t voltage_;
  size_t gate_variable_;

  const Real radius_;
  const Real voltage_target_;

 public:
  RandomPulseTrain(SolidBody* muscle, size_t s1_index, Real radius, Real voltage_target = 0.95)
      : s1_index_(s1_index), radius_(radius), voltage_target_(voltage_target) {
    muscle_        = muscle;
    particles_     = reinterpret_cast<SPH::ElectroPhysiologyParticles*>(muscle->base_particles_);
    voltage_       = particles_->SpeciesIndexMap()["Voltage"];
    gate_variable_ = particles_->SpeciesIndexMap()["GateVariable"];
  }

  void register_extra_stimulus(size_t particle_index) {
    extra_stimuli_indexes_.push_back(particle_index);
    ready_map_[particle_index] = false;
  }

  void generate_and_execute(Real dt) {
    if (0 <= GlobalStaticVariables::physical_time_ && GlobalStaticVariables::physical_time_ <= dt) {
      IndexSphereStimulus s1(muscle_, s1_index_, radius_, voltage_target_);
      s1.parallel_exec(dt);
    }

    auto it = extra_stimuli_indexes_.begin();
    while (it != extra_stimuli_indexes_.end()) {
      if (!ready_map_[*it] && particles_->species_n_[gate_variable_][*it] >= 0.8) {
        ready_map_[*it] = true;
      }
      if (ready_map_[*it] && particles_->species_n_[gate_variable_][*it] <= 0.3) {
        IndexSphereStimulus stimulus(muscle_, *it, radius_, voltage_target_);
        stimulus.parallel_exec(dt);
        ready_map_.erase(*it);
        extra_stimuli_indexes_.erase(it++);
        continue;
      }
      it++;
    }
  }
};