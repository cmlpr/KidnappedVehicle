/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"
#include "helper_functions.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	// Number of particles
	// num_particles = 200;
	num_particles = 200;

	// Generate gaussian normal distributin 
	default_random_engine gen;
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	// Initialize all particles
    for (int i = 0; i < num_particles; ++i) {

		// create a particle instance and update its id, xloc, yloc, theta, weight
	    Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);           // Sample from defined normal distributions
		particle.y = dist_y(gen);           // Sample from defined normal distributions
		particle.theta = dist_theta(gen);   // Sample from defined normal distributions
		particle.weight = 1.0;

		// add the partile and the weight to the arrays
		particles.push_back(particle);
		weights.push_back(1.0);
	}
	// complete initialization so that it doesn't come back to here 
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;

	// for each particle 
	for (int i=0; i < num_particles; i++) 
	{	
		// compute the new location and angle based on the controls data 
		// if yaw rate is zero we have a simplier equation to update the x,y,theta values
		double new_x, new_y, new_theta;
		if (yaw_rate == 0.0) {
			new_x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
			new_y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
			new_theta = particles[i].theta;
		} 
		// when yaw_rate is not zero
		else {
			new_x = particles[i].x + (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
			new_y = particles[i].y + (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
			new_theta = particles[i].theta + yaw_rate * delta_t;
		}

		// add some gaussian noise to controls
		normal_distribution<double> dist_x(new_x, std_pos[0]);
		normal_distribution<double> dist_y(new_y, std_pos[1]);
		normal_distribution<double> dist_theta(new_theta, std_pos[2]);

		// update particles x,y,theta variables with the updated values
		particles[i].x = dist_x(gen);           // Sample from defined normal distributions
		particles[i].y = dist_y(gen);           // Sample from defined normal distributions
		particles[i].theta = dist_theta(gen);   // Sample from defined normal distributions
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	for (int p = 0; p < num_particles; p++) {
		vector<int> associations;
		vector<double> sense_x;
		vector<double> sense_y;

		// First observations are needed to be transformed into map coordinates
		vector<LandmarkObs> trans_observations;
		LandmarkObs obs;
		for (int i = 0; i < observations.size(); i++) {
			
			LandmarkObs trans_obs;
			obs = observations[i];

			// perform space transformation from vehicle to map coordinates 
			// rotation and translation 
			trans_obs.x = particles[p].x + (obs.x * cos(particles[p].theta) - obs.y * sin(particles[p].theta));
			trans_obs.y = particles[p].y + (obs.x * sin(particles[p].theta) + obs.y * cos(particles[p].theta));
			trans_obs.id = obs.id;
			trans_observations.push_back(trans_obs);
		}

		// Initialize the weight to 1
		particles[p].weight = 1.0;

		// For each observation we need to find the closest landmark and get a new weight from the measurements
		for (int i = 0; i < trans_observations.size(); i++) {

			// Assign the closest distance to the sensor range
			double closest_dist = sensor_range;
			// Initialize association id
			int association = 0;
			// Loop thru all the ladmarks and find the closest one
			for (int j=0; j < map_landmarks.landmark_list.size(); j++)
			{	
				// get landmarks coordinates
				double landmark_x = map_landmarks.landmark_list[j].x_f;
				double landmark_y = map_landmarks.landmark_list[j].y_f;

				// calculate the distance between the observation and the landmark
				double calc_dist = dist(trans_observations[i].x, trans_observations[i].y, landmark_x, landmark_y);
				// if the distance is smaller than the closest distance then 
				// update the closest distance and assing the landmark position in the list to association
				if (calc_dist < closest_dist)
				{
					closest_dist = calc_dist;
					association = j;
				}
			}

			// If there is a landmark associated with the observation then multiply weight with the multivariate gaussian probability
			if (association!=0) {
				double meas_x = trans_observations[i].x;
				double meas_y = trans_observations[i].y;
				double mu_x = map_landmarks.landmark_list[association].x_f;
				double mu_y = map_landmarks.landmark_list[association].y_f;

				double first_term = 1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]);
				double x_div = 2.0 * std_landmark[0] * std_landmark[0];
				double y_div = 2.0 * std_landmark[1] * std_landmark[1];
				long double multiplier = first_term * exp(-(pow(meas_x - mu_x, 2.0) / x_div + pow(meas_y - mu_y, 2.0) / y_div));
				
				if (multiplier > 0) {
					particles[p].weight *= multiplier;
				}
			}

			associations.push_back(association+1);
			sense_x.push_back(trans_observations[i].x);
			sense_y.push_back(trans_observations[i].y);
		}

		// Send this for debugging purpose
		particles[p] = SetAssociations(particles[p], associations, sense_x, sense_y);
		weights[p] = particles[p].weight;
	
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	discrete_distribution<int> distribution(weights.begin(), weights.end());

	vector<Particle> resample_particles;

	for (int i = 0; i < num_particles; i++)
	{
		resample_particles.push_back(particles[distribution(gen)]);
	}

	particles = resample_particles;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
