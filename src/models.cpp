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
            // Iterate through each boid
            for (primatives::boid& bi : boids) {
                bi.f = glm::vec3(0.f);
                bi.g = g;

                // Iterate through each other boid
                for (const primatives::boid& bj : boids) {
                    if (&bi == &bj) continue; // Skip self-comparison

                    glm::vec3 DeltaXij = bj.p - bi.p;
                    float d = glm::length(DeltaXij);
//                    std::cout << d << "\n";
                    float alpha = glm::dot(glm::normalize(DeltaXij), glm::normalize(bi.v));

                    if (d < r_s && alpha > cos(theta_s)) {
                        bi.f += separationForce(bi, bj);
                    } else if (d < r_a && alpha > cos(theta_a)) {
                        bi.f += alignmentForce(bi, bj);
                    } else if (d < r_c && alpha > cos(theta_c)) {
                        bi.f += cohesionForce(bi, bj);
                    }
                }

                // Integrate forward velocity and position for bi
//                std::cout << bi.f.x << ", " << bi.f.y << ", " << bi.f.z << "\n";
                bi.integrate(dt);
                bi.orientate();
            }

//            boids[0].v = glm::normalize(glm::vec3(boids[0].p.y * -1.f, boids[0].p.x, 0.f)) * 10.f;
//            boids[1].v = glm::vec3(-5.f, 0.f, 0.f);
//            boids[1].v = glm::vec3(10.f, 0.f, 0.f);
//            boids[2].v = glm::vec3(10.f, 0.f, 0.f);
//            boids[0].v = 10.f * glm::vec3((rand() - RAND_MAX / 1.f) / (RAND_MAX + 1.0), (rand() - RAND_MAX / 1.f) / (RAND_MAX + 1.0), (rand() - RAND_MAX / 1.f) / (RAND_MAX + 1.0));
//            boids[0].integrate(dt);
        }


        void BoidsModel::render(const ModelViewContext& view) {
			//Add Mass render
			for (const primatives::boid& boid : boids)
                givr::addInstance(boid_render, calculateTransformMatrix(boid.p, boid.t, boid.n, boid.b));
//				givr::addInstance(boid_render, glm::translate(glm::mat4(1.f), boid.p)); //NEED TO FRAME!!!

			//Render
			givr::style::draw(boid_render, view);
		}

        glm::vec3 BoidsModel::separationForce(const primatives::boid& bi, const primatives::boid& bj) const {
            glm::vec3 DeltaXij = bj.p - bi.p;
            float d = glm::length(DeltaXij);
            return -1.f * k_s * DeltaXij / (d * d); // Inverse Linear Force
        }

        // Function to calculate alignment force
        glm::vec3 BoidsModel::alignmentForce(const primatives::boid& bi, const primatives::boid& bj) const {
            return k_a * (bj.v - bi.v); // Linear force on velocity
        }

        // Function to calculate cohesion force
        glm::vec3 BoidsModel::cohesionForce(const primatives::boid& bi, const primatives::boid& bj) const {
            return k_c * (bj.p - bi.p); // Linear force
        }

        glm::mat4 BoidsModel::calculateTransformMatrix(glm::vec3 position, glm::vec3 tangent, glm::vec3 normal, glm::vec3 binormal) {
            // Calculate the rotation matrix based on the tangent and normal
            glm::mat4 rotationMatrix = glm::mat4(1.0f);
            rotationMatrix[0] = glm::vec4(binormal, 0.0f);
            rotationMatrix[1] = glm::vec4(normal, 0.0f);
            rotationMatrix[2] = glm::vec4(tangent, 0.0f);
            rotationMatrix[3] = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);

            // Construct the model matrix with the combined rotation
            glm::mat4 TransformMatrix = glm::mat4(1.0f);
            TransformMatrix[0] = rotationMatrix[0];
            TransformMatrix[1] = rotationMatrix[1];
            TransformMatrix[2] = rotationMatrix[2];
            TransformMatrix[3] = glm::vec4(position, 1.0f);

            return TransformMatrix;
        }
	} // namespace models
} // namespace simulation