#pragma once
#include "sphinxsys.h"
#include "cnpy.h"

using namespace SPH;

class BodyStatesRecordingToNpz : public BodyStatesRecording {
 private:
  const int run_id_;
  std::string out_directory_;

 public:
  BodyStatesRecordingToNpz(In_Output& in_output, SPHBodyVector bodies, int run_id)
      : BodyStatesRecording(in_output, bodies), run_id_(run_id) {
    out_directory_ = in_output_.output_folder_ + "/" + std::to_string(run_id_);
    if (!fs::exists(out_directory_)) {
      fs::create_directory(out_directory_);
    }
  };
  virtual ~BodyStatesRecordingToNpz(){};

  virtual void writeWithFileName(const std::string& sequence) override {
    for (SPHBody* body : bodies_) {
      if (body->checkNewlyUpdated()) {
        if (body->getBodyName() == "ExcitationHeart") {
          std::string filefullpath = out_directory_ + "/" + sequence + ".npz";
          if (fs::exists(filefullpath)) {
            fs::remove(filefullpath);
          }

          auto* particles =
              reinterpret_cast<SPH::ElectroPhysiologyParticles*>(body->base_particles_);
          size_t total_real_particles = particles->total_real_particles_;
          if (!particles) {
            throw std::runtime_error("particles is nullptr");
          }

          using RealOut = float;

          cnpy::npz_save<RealOut>(
              filefullpath, "XYZ",
              [total_real_particles, &particles](char* out) {
                for (size_t i = 0; i != total_real_particles; ++i) {
                  Vec3d pos = particles->pos_n_[i];

                  RealOut x;
                  x = pos[0];
                  std::memcpy(out, &x, sizeof(RealOut));
                  x = pos[1];
                  std::memcpy(out + sizeof(RealOut), &x, sizeof(RealOut));
                  x = pos[2];
                  std::memcpy(out + 2 * sizeof(RealOut), &x, sizeof(RealOut));
                  out += 3 * sizeof(RealOut);
                }
              },
              {total_real_particles, 3}, "w");

          size_t voltage = particles->SpeciesIndexMap()["Voltage"];
          cnpy::npz_save<RealOut>(
              filefullpath, "Voltage",
              [&particles, k = voltage](char* out) {
                auto out_f = reinterpret_cast<RealOut*>(out);
                std::copy(particles->species_n_[k].begin(), particles->species_n_[k].end(), out_f);
              },
              {total_real_particles, 1}, "a");
        }
      }
      body->setNotNewlyUpdated();
    }
  }
};

// class ReloadParticleIO_NPZ : public ReloadParticleIO {
//  public:
//   ReloadParticleIO_NPZ(In_Output& in_output, SPHBodyVector bodies)
//       : ReloadParticleIO(in_output, bodies) {
//     for (size_t i = 0; i != bodies.size(); ++i) {
//       file_paths_[i] =
//           in_output.reload_folder_ + "/SPHBody_" + bodies[i]->getBodyName() + "_rld.npz";
//     }
//   }
//   ReloadParticleIO_NPZ(In_Output& in_output,
//                        SPHBodyVector bodies,
//                        StdVec<std::string> given_body_names)
//       : ReloadParticleIO(in_output, bodies) {
//     if (bodies.size() != given_body_names.size()) {
//       throw std::logic_error("bodies.size() != given_body_names.size()");
//     }
//     for (size_t i = 0; i != bodies.size(); ++i) {
//       file_paths_[i] = in_output.reload_folder_ + "/SPHBody_" + given_body_names[i] + "_rld.npz";
//     }
//   }
//   ~ReloadParticleIO_NPZ() override = default;

//   void writeToFile(size_t iteration_step = 0) override {
//     for (size_t i = 0; i < bodies_.size(); ++i) {
//       std::string filefullpath = file_paths_[i];

//       if (fs::exists(filefullpath)) {
//         fs::remove(filefullpath);
//       }
//       bodies_[i]->writeToXmlForReloadParticle(filefullpath);
//       auto particles            = bodies_[i]->base_particles_;
//       auto total_real_particles = particles->total_real_particles_;
//       using RealOut             = Real;
//       cnpy::npz_save<RealOut>(
//           filefullpath, "Position",
//           [total_real_particles, &particles](char* out) {
//             for (size_t i = 0; i != total_real_particles; ++i) {
//               Vec3d pos = particles->pos_n_[i];

//               RealOut x;
//               x = pos[0];
//               std::memcpy(out, &x, sizeof(RealOut));
//               x = pos[1];
//               std::memcpy(out + sizeof(RealOut), &x, sizeof(RealOut));
//               x = pos[2];
//               std::memcpy(out + 2 * sizeof(RealOut), &x, sizeof(RealOut));
//               out += 3 * sizeof(RealOut);
//             }
//           },
//           {total_real_particles, 3}, "w");
//       cnpy::npz_save<RealOut>(
//           filefullpath, "Volume",
//           [&particles](char* out) {
//             auto out_f = reinterpret_cast<RealOut*>(out);
//             std::copy(particles->Vol_.begin(), particles->Vol_.end(), out_f);
//           },
//           {total_real_particles, 1}, "a");
//     }
//   };
//   void readFromFile(size_t iteration_step = 0) override {
//     std::cout << "\n Reloading particles from files." << std::endl;
//     for (size_t i = 0; i < bodies_.size(); ++i) {
//       std::string filefullpath = file_paths_[i];

//       if (!fs::exists(filefullpath)) {
//         std::cout << "\n Error: the input file:" << filefullpath << " does not exist" << std::endl;
//         std::cout << __FILE__ << ':' << __LINE__ << std::endl;
//         exit(1);
//       }

//       auto pos      = cnpy::npz_load(filefullpath, "Position");
//       auto pos_data = pos.data<Vec3d>();
//       std::copy(pos_data, pos_data + pos.shape[0], bodies_[i]->base_particles_->pos_n_.begin());

//       auto vol      = cnpy::npz_load(filefullpath, "Volume");
//       auto vol_data = vol.data<Real>();
//       std::copy(vol_data, vol_data + vol.shape[0], bodies_[i]->base_particles_->Vol_.begin());
//     }
//   };
// };

// class ParticleGeneratorReloadNPZ : public ParticleGeneratorReload {
//   std::string file_path_;
// class ReloadParticleIO_NPZ : public ReloadParticleIO {
//  public:
//   ReloadParticleIO_NPZ(In_Output& in_output, SPHBodyVector bodies)
//       : ReloadParticleIO(in_output, bodies) {
//     for (size_t i = 0; i != bodies.size(); ++i) {
//       file_paths_[i] =
//           in_output.reload_folder_ + "/SPHBody_" + bodies[i]->getBodyName() + "_rld.npz";
//     }
//   }
//   ReloadParticleIO_NPZ(In_Output& in_output,
//                        SPHBodyVector bodies,
//                        StdVec<std::string> given_body_names)
//       : ReloadParticleIO(in_output, bodies) {
//     if (bodies.size() != given_body_names.size()) {
//       throw std::logic_error("bodies.size() != given_body_names.size()");
//     }
//     for (size_t i = 0; i != bodies.size(); ++i) {
//       file_paths_[i] = in_output.reload_folder_ + "/SPHBody_" + given_body_names[i] + "_rld.npz";
//     }
//   }
//   ~ReloadParticleIO_NPZ() override = default;

//   void writeToFile(size_t iteration_step = 0) override {
//     for (size_t i = 0; i < bodies_.size(); ++i) {
//       std::string filefullpath = file_paths_[i];

//       if (fs::exists(filefullpath)) {
//         fs::remove(filefullpath);
//       }
//       bodies_[i]->writeToXmlForReloadParticle(filefullpath);
//       auto particles            = bodies_[i]->base_particles_;
//       auto total_real_particles = particles->total_real_particles_;
//       using RealOut             = Real;
//       cnpy::npz_save<RealOut>(
//           filefullpath, "Position",
//           [total_real_particles, &particles](char* out) {
//             for (size_t i = 0; i != total_real_particles; ++i) {
//               Vec3d pos = particles->pos_n_[i];

//               RealOut x;
//               x = pos[0];
//               std::memcpy(out, &x, sizeof(RealOut));
//               x = pos[1];
//               std::memcpy(out + sizeof(RealOut), &x, sizeof(RealOut));
//               x = pos[2];
//               std::memcpy(out + 2 * sizeof(RealOut), &x, sizeof(RealOut));
//               out += 3 * sizeof(RealOut);
//             }
//           },
//           {total_real_particles, 3}, "w");
//       cnpy::npz_save<RealOut>(
//           filefullpath, "Volume",
//           [&particles](char* out) {
//             auto out_f = reinterpret_cast<RealOut*>(out);
//             std::copy(particles->Vol_.begin(), particles->Vol_.end(), out_f);
//           },
//           {total_real_particles, 1}, "a");
//     }
//   };
//   void readFromFile(size_t iteration_step = 0) override {
//     std::cout << "\n Reloading particles from files." << std::endl;
//     for (size_t i = 0; i < bodies_.size(); ++i) {
//       std::string filefullpath = file_paths_[i];

//       if (!fs::exists(filefullpath)) {
//  public:
//   ParticleGeneratorReloadNPZ(In_Output* in_output, std::string reload_body_name)
//       : ParticleGeneratorReload(in_output, reload_body_name) {
//     file_path_ = in_output->reload_folder_ + "/SPHBody_" + reload_body_name + "_rld.npz";
//   }
//   ~ParticleGeneratorReloadNPZ() override = default;

//   void createBaseParticles(BaseParticles* base_particles) override {
//     if (!fs::exists(file_path_)) {
//       std::cout << "\n Error: the reload file:" << file_path_ << " does not exist" << std::endl;
//       std::cout << __FILE__ << ':' << __LINE__ << std::endl;
//       std::exit(1);
//     }
//     auto pos      = cnpy::npz_load(file_path_, "Position");
//     auto pos_data = pos.data<Vec3d>();
//     auto vol      = cnpy::npz_load(file_path_, "Volume");
//     auto vol_data = vol.data<Real>();

//     for (size_t i = 0; i < pos.shape[0]; i++) {
//       base_particles->initializeABaseParticle(pos_data[i], vol_data[i]);
//     }
//   }
// };