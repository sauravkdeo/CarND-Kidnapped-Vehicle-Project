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
#include "map.h"

using std::string;
using std::vector;
//struct single_landmark_s;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	/**
	 * TODO: Set the number of particles. Initialize all particles to
	 *   first position (based on estimates of x, y, theta and their uncertainties
	 *   from GPS) and all weights to 1.
	 * TODO: Add random Gaussian noise to each particle.
	 * NOTE: Consult particle_filter.h for more information about this method
	 *   (and others in this file).
	 */
	if(is_initialized) {
		return;
	}

	// Initializing the number of particles
	num_particles = 100;
	// Extracting standard deviations
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	// Creating normal distributions
	std::normal_distribution<double> dist_x(x, std_x);
	std::normal_distribution<double> dist_y(y, std_y);
	std::normal_distribution<double> dist_theta(theta, std_theta);

	// Generate particles with normal distribution with mean on GPS values.
	for (int i = 0; i < num_particles; i++) {

		Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1.0;

		particles.push_back(particle);
	}

	// The filter is now initialized.
	is_initialized = true;
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

	// Extracting standard deviations
	double std_x = std_pos[0];
	double std_y = std_pos[1];
	double std_theta = std_pos[2];

	// Creating normal distributions
	std::normal_distribution<double> dist_x(0, std_x);
	std::normal_distribution<double> dist_y(0, std_y);
	std::normal_distribution<double> dist_theta(0, std_theta);

	// Calculate new state.
	for(std::vector<Particle> ::iterator  p_it = particles.begin(); p_it != particles.end(); ++p_it) {

		double theta = p_it->theta;

		if ( fabs(yaw_rate) < 0.0001 ) { // When yaw is not changing.
			p_it->x += velocity * delta_t * cos( theta );
			p_it->y += velocity * delta_t * sin( theta );
			// yaw continue to be the same.
		} else {
			p_it->x += velocity / yaw_rate * ( sin( theta + yaw_rate * delta_t ) - sin( theta ) );
			p_it->y += velocity / yaw_rate * ( cos( theta ) - cos( theta + yaw_rate * delta_t ) );
			p_it->theta += yaw_rate * delta_t;
		}

		// Adding noise.
		p_it->x += dist_x(gen);
		p_it->y += dist_y(gen);
		p_it->theta += dist_theta(gen);
	}
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


	for (vector<LandmarkObs>::iterator obs_it = observations.begin();
			obs_it != observations.end(); ++obs_it) { // For each observation

		// Initialize min distance as a really big number.
		double minDistance = std::numeric_limits<double>::max();

		// Initialize the found map in something not possible.
		int mapId = -1;

		for (vector<LandmarkObs>::iterator pred_it = predicted.begin();
				pred_it != predicted.end(); ++pred_it) {

			double xDistance = obs_it->x - pred_it->x;
			double yDistance = obs_it->y - pred_it->y;

			double distance = xDistance * xDistance + yDistance * yDistance;

			// If the "distance" is less than min, stored the id and update min.
			if ( distance < minDistance ) {
				minDistance = distance;
				mapId = pred_it->id;
			}
		}

		// Update the observation identifier.
		obs_it->id = mapId;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const vector<LandmarkObs> &observations,
		const Map &map_landmarks) {
	/**
	 * TODO: Update the weights of each particle using a mult-variate Gaussian
	 *   distribution. You can read more about this distribution here:
	 *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	 * NOTE: The observations are given in the VEHICLE'S coordinate system.
	 *   Your particles are located according to the MAP'S coordinate system.
	 *   You will need to transform between the two systems. Keep in mind that
	 *   this transformation requires both rotation AND translation (but no scaling).
	 *   The following is a good resource for the theory:
	 *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	 *   and the following is a good resource for the actual equation to implement
	 *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
	 */
	double stdLandmarkRange = std_landmark[0];
	double stdLandmarkBearing = std_landmark[1];

	for (std::vector<Particle>::iterator Particle_it = particles.begin();
			Particle_it != particles.end(); ++Particle_it) {

		double x = Particle_it->x;
		double y = Particle_it->y;
		double theta = Particle_it->theta;
		// Find landmarks in particle's range.
		double sensor_range_2 = sensor_range * sensor_range;
		vector<LandmarkObs> inRangeLandmarks;
		//		Map::single_landmark_s landmark_it;
		for(auto landmark_it = map_landmarks.landmark_list.begin();
				landmark_it != map_landmarks.landmark_list.end();++landmark_it)	{
			float landmarkX = landmark_it->x_f;
			float landmarkY = landmark_it->y_f;
			int id = landmark_it->id_i;
			double dX = x - landmarkX;
			double dY = y - landmarkY;
			if ( dX*dX + dY*dY <= sensor_range_2 ) {
				inRangeLandmarks.push_back(LandmarkObs{ id, landmarkX, landmarkY });
			}
		}

		// Transform observation coordinates.
		vector<LandmarkObs> mappedObservations;
		for(vector<LandmarkObs>::const_iterator obs_it = observations.begin();
				obs_it != observations.end();++obs_it)	{
			double xx = cos(theta)*obs_it->x - sin(theta)*obs_it->y + x;
			double yy = sin(theta)*obs_it->x + cos(theta)*obs_it->y + y;
			mappedObservations.push_back(LandmarkObs{ obs_it->id, xx, yy });
		}

		// Observation association to landmark.
		dataAssociation(inRangeLandmarks, mappedObservations);

		// Reseting weight.
		Particle_it->weight = 1.0;
		// Calculate weights.
		for(vector<LandmarkObs>::iterator m_obs_it = mappedObservations.begin();
				m_obs_it != mappedObservations.end();++m_obs_it)	{
			double observationX = m_obs_it->x;
			double observationY = m_obs_it->y;

			int landmarkId = m_obs_it->id;

			double landmarkX, landmarkY;
			unsigned int k = 0;
			unsigned int nLandmarks = inRangeLandmarks.size();
			bool found = false;
			while( !found && k < nLandmarks ) {
				if ( inRangeLandmarks[k].id == landmarkId) {
					found = true;
					landmarkX = inRangeLandmarks[k].x;
					landmarkY = inRangeLandmarks[k].y;
				}
				k++;
			}

			// Calculating weight.
			double dX = observationX - landmarkX;
			double dY = observationY - landmarkY;

			double weight = ( 1/(2*M_PI*stdLandmarkRange*stdLandmarkBearing)) * exp( -( dX*dX/(2*stdLandmarkRange*stdLandmarkRange) + (dY*dY/(2*stdLandmarkBearing*stdLandmarkBearing)) ) );
			if (weight == 0) {
				Particle_it->weight *= 0.0001;
			} else {
				Particle_it->weight *= weight;
			}
		}
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	// Get weights and max weight.
	vector<double> weights;
	double maxWeight = std::numeric_limits<double>::min();
	for(std::vector<Particle> ::iterator  p_it = particles.begin(); p_it != particles.end(); ++p_it) {
		weights.push_back(p_it->weight);
		if ( p_it->weight > maxWeight ) {
			maxWeight = p_it->weight;
		}
	}

	// Creating distributions.
	std::uniform_real_distribution<double> distDouble(0.0, maxWeight);
	std::uniform_int_distribution<int> distInt(0, num_particles - 1);

	// Generating index.
	int index = distInt(gen);

	double beta = 0.0;

	// the wheel
	vector<Particle> resampledParticles;

	for(std::vector<Particle> ::iterator  p_it = particles.begin(); p_it != particles.end(); ++p_it) {
		beta += distDouble(gen) * 2.0;
		while( beta > weights[index]) {
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}
		resampledParticles.push_back(particles[index]);
	}

	particles = resampledParticles;
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
