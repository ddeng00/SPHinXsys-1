#pragma once

#include "sphinxsys.h"

using namespace SPH;

class VoltageObserver : public FictitiousBody {
 public:
  VoltageObserver(SPHSystem& system, string body_name, Real length_scale)
      : FictitiousBody(system, body_name) {

    body_input_points_volumes_.push_back(
        std::make_pair(Vecd(-45.0 * length_scale, -30.0 * length_scale, 0.0), 0.0));
    body_input_points_volumes_.push_back(
        std::make_pair(Vecd(0.0, -30.0 * length_scale, 26.0 * length_scale), 0.0));
    body_input_points_volumes_.push_back(
        std::make_pair(Vecd(-30.0 * length_scale, -50.0 * length_scale, 0.0), 0.0));
    body_input_points_volumes_.push_back(
        std::make_pair(Vecd(0.0, -50.0 * length_scale, 20.0 * length_scale), 0.0));
    body_input_points_volumes_.push_back(std::make_pair(Vecd(0.0, -70.0 * length_scale, 0.0), 0.0));
  }
};

class MyocardiumObserver : public FictitiousBody {
 public:
  MyocardiumObserver(SPHSystem& system, std::string body_name, Real length_scale)
      : FictitiousBody(system, body_name) {

    body_input_points_volumes_.push_back(
        std::make_pair(Vecd(-45.0 * length_scale, -30.0 * length_scale, 0.0), 0.0));
    body_input_points_volumes_.push_back(
        std::make_pair(Vecd(0.0, -30.0 * length_scale, 26.0 * length_scale), 0.0));
    body_input_points_volumes_.push_back(
        std::make_pair(Vecd(-30.0 * length_scale, -50.0 * length_scale, 0.0), 0.0));
    body_input_points_volumes_.push_back(
        std::make_pair(Vecd(0.0, -50.0 * length_scale, 20.0 * length_scale), 0.0));
    body_input_points_volumes_.push_back(std::make_pair(Vecd(0.0, -70.0 * length_scale, 0.0), 0.0));
  }
};