#include "solid_particles.h"
#include "solid_particles_variable.h"
#include "base_body.h"

#include <iterator>

namespace SPH
{
	//=============================================================================================//
	void SolidParticles::ParticleTranslationAndRotation(Transformd &transform)
	{
		std::cout << "\n Error: the function ParticleTranslationAndRotation in 3d is not defined!" << std::endl;
		std::cout << __FILE__ << ':' << __LINE__ << std::endl;
		exit(1);
	}
	//=================================================================================================//
	Matd ElasticSolidParticles::get_GreenLagrange_strain(size_t particle_i)
	{
		Mat3d F = F_[particle_i];
		return 0.5 * (~F * F - Matd(1.0)); // calculation of the Green-Lagrange strain tensor
	}
	//=================================================================================================//
	Vecd ElasticSolidParticles::get_Principal_strains(size_t particle_i)
	{
		Mat3d epsilon = get_GreenLagrange_strain(particle_i); // calculation of the Green-Lagrange strain tensor
		return getPrincipalValuesFromMatrix(epsilon);
	}
	//=================================================================================================//
	Real ElasticSolidParticles::von_Mises_strain_static(size_t particle_i)
	{
		Mat3d epsilon = get_GreenLagrange_strain(particle_i); // calculation of the Green-Lagrange strain tensor

		Real epsilonxx = epsilon(0, 0);
		Real epsilonyy = epsilon(1, 1);
		Real epsilonzz = epsilon(2, 2);
		Real epsilonxy = epsilon(0, 1);
		Real epsilonxz = epsilon(0, 2);
		Real epsilonyz = epsilon(1, 2);

		return sqrt( (1.0 / 3.0) * (powerN(epsilonxx - epsilonyy, 2) + powerN(epsilonyy - epsilonzz, 2) + powerN(epsilonzz - epsilonxx, 2))
		 + 2.0 * (powerN(epsilonxy, 2) + powerN(epsilonyz, 2) + powerN(epsilonxz, 2)));
	}
	//=================================================================================================//
	Real ElasticSolidParticles::von_Mises_strain_dynamic(size_t particle_i, Real poisson)
	{
		Mat3d epsilon = get_GreenLagrange_strain(particle_i); // calculation of the Green-Lagrange strain tensor
		
		Vec3d principal_strains = getPrincipalValuesFromMatrix(epsilon);
		Real eps_1 = principal_strains[0];
		Real eps_2 = principal_strains[1];
		Real eps_3 = principal_strains[2];

		return 1.0/(1.0 + poisson) * std::sqrt(0.5 * (powerN(eps_1 - eps_2, 2) + powerN(eps_2 - eps_3, 2) + powerN(eps_3 - eps_1, 2)));
	}
	//=================================================================================================//
	Matd ElasticSolidParticles::get_Cauchy_stress(size_t particle_i)
	{
		Mat3d F = F_[particle_i];
		Real J = det(F);
		Mat3d stress = stress_PK1_[particle_i];

		return ( 1.0 / J ) * F * stress * ~F; // Cauchy stress 
	}
	//=================================================================================================//
	Matd ElasticSolidParticles::get_PK2_stress(size_t particle_i)
	{
		Mat3d F = F_[particle_i];
		Mat3d stress = stress_PK1_[particle_i];

		return SimTK::inverse(F) * stress; // Second Piola-Kirchhof stress
	}
	//=================================================================================================//
	Vecd ElasticSolidParticles::get_Principal_stresses(size_t particle_i)
	{
		Mat3d sigma;
		if (stress_measure_ == "Cauchy") {
			sigma = get_Cauchy_stress(particle_i); // Cauchy stress
		} else if (stress_measure_ == "PK2") {
			sigma = get_PK2_stress(particle_i); // Second Piola-Kirchhof stress
		} else {
			throw std::runtime_error("get_Principal_stresses: wrong input");
		}

		return getPrincipalValuesFromMatrix(sigma);
	}
	//=================================================================================================//
	Real ElasticSolidParticles::get_von_Mises_stress(size_t particle_i)
	{
		Mat3d sigma;
		if (stress_measure_ == "Cauchy") {
			sigma = get_Cauchy_stress(particle_i); // Cauchy stress
		} else if (stress_measure_ == "PK2") {
			sigma = get_PK2_stress(particle_i); // Second Piola-Kirchhof stress
		} else {
			throw std::runtime_error("get_von_Mises_stress: wrong input");
		}

		return getVonMisesStressFromMatrix(sigma);
	}
	//=================================================================================================//
	Real ElasticSolidParticles::von_Mises_strain(size_t particle_i)
	{

		Mat3d F = F_[particle_i];
		Mat3d epsilon = 0.5 * (~F * F - Matd(1.0)); //calculation of the Green-Lagrange strain tensor

		Real epsilonxx = epsilon(0, 0);
		Real epsilonyy = epsilon(1, 1);
		Real epsilonzz = epsilon(2, 2);
		Real epsilonxy = epsilon(0, 1);
		Real epsilonxz = epsilon(0, 2);
		Real epsilonyz = epsilon(1, 2);

		return sqrt((1.0 / 3.0) * (std::pow(epsilonxx - epsilonyy, 2.0) + std::pow(epsilonyy - epsilonzz, 2.0) + std::pow(epsilonzz - epsilonxx, 2.0)) + 2.0 * (std::pow(epsilonxy, 2.0) + std::pow(epsilonyz, 2.0) + std::pow(epsilonxz, 2.0)));
	}
	//=============================================================================================//
	void VonMisesStress::operator()(size_t index_i, Real dt)
	{
		Real J = rho0_ / rho_n_[index_i];
		Mat3d F = F_[index_i];
		Mat3d stress = stress_PK1_[index_i];
		Mat3d sigma = (stress * ~F) / J;

		Real sigmaxx = sigma(0, 0);
		Real sigmayy = sigma(1, 1);
		Real sigmazz = sigma(2, 2);
		Real sigmaxy = sigma(0, 1);
		Real sigmaxz = sigma(0, 2);
		Real sigmayz = sigma(1, 2);

		derived_variable_[index_i] =
			sqrt(sigmaxx * sigmaxx + sigmayy * sigmayy + sigmazz * sigmazz -
				 sigmaxx * sigmayy - sigmaxx * sigmazz - sigmayy * sigmazz +
				 3.0 * (sigmaxy * sigmaxy + sigmaxz * sigmaxz + sigmayz * sigmayz));
	}
	//=============================================================================================//
	void VonMisesStrain::operator()(size_t index_i, Real dt)
	{
		Mat3d F = F_[index_i];
		Mat3d epsilon = 0.5 * (~F * F - Matd(1.0)); //calculation of the Green-Lagrange strain tensor

		Real epsilonxx = epsilon(0, 0);
		Real epsilonyy = epsilon(1, 1);
		Real epsilonzz = epsilon(2, 2);
		Real epsilonxy = epsilon(0, 1);
		Real epsilonxz = epsilon(0, 2);
		Real epsilonyz = epsilon(1, 2);

		derived_variable_[index_i] =
			sqrt((1.0 / 3.0) * (powerN(epsilonxx - epsilonyy, 2) + powerN(epsilonyy - epsilonzz, 2) +
								powerN(epsilonzz - epsilonxx, 2)) +
				 2.0 * (powerN(epsilonxy, 2) + powerN(epsilonyz, 2) + powerN(epsilonxz, 2)));
	}
	//=================================================================================================//
}
