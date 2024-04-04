#include "models.hpp"
#include <random>

namespace simulation {
	namespace primatives {
		//No functions for structs, yet...
	}// namespace primatives

	namespace models {
		//Unless changed in main, this will call only once before the window is created
		BoidsModel::BoidsModel()
			: boid_geometry(givr::geometry::Mesh(givr::geometry::Filename("./models/dart.obj")))
			, boid_style(givr::style::Colour(1.f, 1.f, 0.f), givr::style::LightPosition(100.f, 100.f, 100.f))
		{
			// Reset Dynamic elements
			reset();

			// Render
			boid_render = givr::createInstancedRenderable(boid_geometry, boid_style);
		}

		void BoidsModel::reset() {
			// static so they are persistent
			static std::random_device random_device;
			static std::mt19937 generator(random_device());
			static std::uniform_real_distribution<double> position_distribution(-10., 10.);
			static std::uniform_real_distribution<double> theta_distribution(-180., 180.);
			static std::uniform_real_distribution<double> phi_distribution(-90., 90.);
			static std::uniform_real_distribution<double> speed_distribution(5., 25.);

			boids.resize(n_boids);
			for (primatives::boid& boid : boids) {
				//Random Position in 20 x 20 x 20 cube
				boid.p.x = position_distribution(generator);
				boid.p.y = position_distribution(generator);
				boid.p.z = position_distribution(generator);

				//Random heading (two angles) and Random speed
				double theta = glm::radians(theta_distribution(generator));
				double phi = glm::radians(phi_distribution(generator));
				double speed = speed_distribution(generator);

				// https://stackoverflow.com/questions/30011741/3d-vector-defined-by-2-angles
				boid.v.x = speed * std::cos(theta) * std::cos(phi);
				boid.v.y = speed * std::sin(phi);
				boid.v.z = speed * std::sin(theta) * std::cos(phi);

			}

		}

		void BoidsModel::step(float dt) {
			reset(); // Comment this out
			//TODO: apply forces and compute integration step
		}

		void BoidsModel::render(const ModelViewContext& view) {
			//Add Mass render
			for (const primatives::boid& boid : boids)
				givr::addInstance(boid_render, glm::translate(glm::mat4(1.f), boid.p)); //NEED TO FRAME!!!

			//Render
			givr::style::draw(boid_render, view);
		}
	} // namespace models
} // namespace simulation