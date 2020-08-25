/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"
//#include "multiv_gauss.h"
#include <cmath>

using std::string;
using std::vector;

using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  
  std::default_random_engine gen;
  
  num_particles = 10;  // TODO: Set the number of particles
  
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  //std::cout << num_particles << std::endl;
  
  
  for (int i = 0; i < num_particles; ++i) {
    double sample_x, sample_y, sample_theta;
    //weights.push_back(1./100.);
    weights.push_back(1.);
    
    sample_x = dist_x(gen);
    sample_y = dist_y(gen);
    sample_theta = dist_theta(gen);
    
    //Particle particle = {i, sample_x, sample_y, sample_theta, 1./100.};
    Particle particle;
    particle.id = i;
    particle.x = sample_x;
    particle.y = sample_y;
    particle.theta = sample_theta;
    particle.weight = 1.;
    
    particles.push_back(particle);
    
  }
  
  is_initialized = true;
  
  //std::cout << "ParticleFilter::init" << std::endl;
  
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  std::default_random_engine gen;
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
  
  for (int i = 0; i < num_particles; ++i) {
  
    Particle particle = particles[i];
    double x = particle.x;
    double y = particle.y;
    double theta = particle.theta;
    
    double yaw = yaw_rate * delta_t;
    
    
    
    if (fabs(yaw_rate) < 0.01 ) { // Drive stright, when yaw not change
      x += velocity * delta_t * cos(theta);
      y += velocity * delta_t * sin(theta);
    } 
    else 
    {
      x += ((velocity/yaw_rate)*(sin(theta+yaw)-sin(theta)));
      y += ((velocity/yaw_rate)*(-cos(theta+yaw)+cos(theta)));
      theta += yaw;
    }
    
    double sample_x = x+dist_x(gen);
    double sample_y = y+dist_y(gen);
    double sample_theta = theta+dist_theta(gen);
    
    particle.x = sample_x;
    particle.y = sample_y;
    particle.theta = sample_theta;
    
    particles[i] = particle;
  
  }
  
//std::cout << "ParticleFilter::prediction" << std::endl;
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  for (int i = 0; i < observations.size(); ++i) {
    
    double obs_x = observations[i].x;
    double obs_y = observations[i].y;
    
    double lowest_distance = 1000000000.;
    int selected_id = -1;
    
    double distance;
    for (int k = 0; k < predicted.size(); ++k) {
      double pred_x = predicted[k].x;
      double pred_y = predicted[k].y;
      int pred_id = predicted[k].id;
      
      //distance = sqrt(pow((obs_x-pred_x), 2.0) + pow((obs_y-pred_y), 2.0));
      distance = dist(obs_x, obs_y, pred_x, pred_y);
      /**std::cout << obs_x << '; ' << obs_y << std::endl;
      std::cout << pred_x << '; ' << pred_y << std::endl;
      std::cout << distance << std::endl;**/
      if (distance < lowest_distance){
        lowest_distance = distance;
        selected_id = pred_id;
        //observations[i].id = k;
      }
    }
    /**std::cout << lowest_distance << std::endl;
    std::cout << selected_id << std::endl;
    std::cout << std::endl;**/
    observations[i].id = selected_id;
  }  
  //std::cout << "ParticleFilter::dataAssociation" << std::endl;
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations   are given in the VEHICLE'S coordinate system. 
   *   Yo particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
    struct single_landmark_s {
      int id_i ; // Landmark ID
      float x_f; // Landmark x-position in the map (global coordinates)
      float y_f; // Landmark y-position in the map (global coordinates)
    };
  
  double x_part, y_part, x_obs, y_obs, theta, distance;
  
  double sig_x = std_landmark[0];
  double sig_y = std_landmark[1];
  
  
  /**for each particle we need to convert all observations to map coordinates
  #afterwards we need to check which of the map landmarks in the sensor range of that particle and include them in a vector of LandmarkObs objects
  #now we are able to use dataAssociation function**/
  for (int i = 0; i < num_particles; ++i) {
  
    Particle particle = particles[i];
   
    x_part = particle.x;
    y_part = particle.y;
    theta = particle.theta;
    
    vector<LandmarkObs> observations_transformed;
    
    for (int k = 0; k < observations.size(); ++k) {
      LandmarkObs obs = observations[k];
      x_obs = obs.x;
      y_obs = obs.y;
      
      // transform to map x coordinate
      double x_map;
      x_map = x_part + (cos(theta) * x_obs) - (sin(theta) * y_obs);

      // transform to map y coordinate
      double y_map;
      y_map = y_part + (sin(theta) * x_obs) + (cos(theta) * y_obs);
      
      obs.x = x_map;
      obs.y = y_map;
      
      observations_transformed.push_back(obs);
      
    }
    
    vector<LandmarkObs> landmarks_in_range;
    LandmarkObs lm;
    
    //check which landmarks are in sensor range
    for (int k = 0; k < map_landmarks.landmark_list.size(); ++k) {
      //single_landmark_s landmark = map_landmarks.landmark_list[k];
      double x_lm = map_landmarks.landmark_list[k].x_f;
      double y_lm = map_landmarks.landmark_list[k].y_f;
      
      //distance = sqrt(pow((x_lm-x_part), 2.0) + pow((y_lm-y_part), 2.0));
      distance = dist(x_lm, y_lm, x_part, y_part);
      
      if (distance < sensor_range){
        lm.x = x_lm;
        lm.y = y_lm;
        lm.id = map_landmarks.landmark_list[k].id_i;
        landmarks_in_range.push_back(lm);
      } 
    }
    
    if (landmarks_in_range.size() ==0){
      //std::cout << "particle " << i << std::endl;
      //std::cout << particle.x << particle.y << std::endl;
      /**double numbers[3] {0};
      numbers[0] = 0.3;
      numbers[1] = 0.3;
      numbers[2] = 0.01;
      init(best_particle_.x, best_particle_.y, best_particle_.theta, numbers);
      
      for (int k = 0; k < map_landmarks.landmark_list.size(); ++k) {
        double x_lm = map_landmarks.landmark_list[k].x_f;
        double y_lm = map_landmarks.landmark_list[k].y_f;
        distance = dist(x_lm, y_lm, x_part, y_part);
      
        if (distance < sensor_range){
        
          lm.x = x_lm;
          lm.y = y_lm;
          lm.id = map_landmarks.landmark_list[k].id_i;
          landmarks_in_range.push_back(lm);
        } 
      }
      **/
    }
    
    
    
    //use data association function with observations_transformed and landmarks_in_range
    dataAssociation(landmarks_in_range, observations_transformed);
    //std::cout << "particle" << i << std::endl;
    //std::cout << particle.x << particle.y << std::endl;
    
    
    LandmarkObs obs;
    particle.weight = 1.;
    for (int k = 0; k < landmarks_in_range.size(); ++k) {
      lm = landmarks_in_range[k];
      for (int q = 0; q < observations_transformed.size(); ++q) {
        obs = observations_transformed[q];
        if(lm.id == obs.id){
          particle.weight *= multiv_prob(sig_x, sig_y, obs.x, obs.y, lm.x, lm.y);
        }
      }
    }
    weights[i] = particle.weight;
    particles[i].weight = particle.weight;
  }
  
  /**
  double normaliser = 0.0;
  for (int i = 0; i < num_particles; ++i) {
    normaliser += weights[i];
  }
  for (int i = 0; i < num_particles; ++i) {
    weights[i] = weights[i]/normaliser;
    particles[i].weight = weights[i]/normaliser;
  }**/
  
  //std::cout << "ParticleFilter::updateWeights" << std::endl;
}


void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  std::cout << "particle 0" << std::endl;
  std::cout << particles[0].x << particles[0].y << std::endl;
  
  
  std::vector<Particle> newparticles;
  std::vector<double> newweights;
  
  std::discrete_distribution<> distribution(weights.begin(), weights.end());

  std::default_random_engine gen;
  for (int i = 0; i < num_particles; ++i) {
    
    int idx = distribution(gen);
    std::cout << idx << std::endl;
    std::cout << particles[idx].x << particles[idx].y << std::endl;
    std::cout << std::endl;
    
    newparticles.push_back(particles[idx]);
    newweights.push_back(newparticles[i].weight);
  }
  
  particles = newparticles;
  weights = newweights;
  
  //std::cout << "particle 0" << std::endl;
  //std::cout << particles[0].x << particles[0].y << std::endl;
  
  /**for (int w=0; w<newweights.size() ; ++w){
    std::cout << newweights[w] << std::endl;
  }**/
  
  
  //std::cout << "ParticleFilter::resample" << std::endl;
  
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

double ParticleFilter::multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}
