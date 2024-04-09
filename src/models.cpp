#include "models.hpp"
#include <random>
#include <omp.h>

namespace simulation {
	namespace primatives {
		//No functions for structs, yet...
	}// namespace primatives

	namespace models {
		//Unless changed in main, this will call only once before the window is created
		BoidsModel::BoidsModel()
			: boid_geometry(givr::geometry::Mesh(givr::geometry::Filename("./models/dart.obj")))
			, boid_style(givr::style::Colour(1.f, 1.f, 0.f), givr::style::LightPosition(100.f, 100.f, 100.f), givr::style::AmbientFactor(0.25f))
            , wall_geometry()
            , wall_style(givr::style::Colour(1.f, 1.f, 1.f), givr::style::LightPosition(100.f, 100.f, 100.f))
            , sphere_geometry()
            , sphere_style(givr::style::Colour(0.f, 1.f, 1.f), givr::style::LightPosition(100.f, 100.f, 100.f), givr::style::AmbientFactor(0.25f))
		{
            int plane_count = 6;
            planes.resize(plane_count);
            for (int i = 0; i < plane_count; ++i) {
                primatives::plane pl = primatives::plane();
                if (i < 2) {
                    pl.origin = glm::vec3((i - 0.5f) * grid_w, 0.f, 0.f);
                } else if (i < 4) {
                    pl.origin = glm::vec3(0.f, (i - 2.5f) * grid_h, 0.f);
                } else {
                    pl.origin = glm::vec3(0.f, 0.f, (i - 4.5f) * grid_d);
                }
//                std::cout << pl.origin.x << ", " << pl.origin.y << ", " << pl.origin.z << "\n";
                pl.normal = glm::normalize(glm::vec3(0.f) - pl.origin);
                planes[i] = pl;
            }

            int sphere_count = 1;
            spheres.resize(sphere_count);
            for (int i = 0; i < sphere_count; ++i) {
                primatives::sphere sp = primatives::sphere();
                sp.origin = glm::vec3(-0.f, -25.f, -25.f);
                sp.radius = 15.f;
                sp.epsilon = 2.f;
                sp.collision_avoid_distance = 5.f;
                spheres[i] = sp;
            }

			// Reset Dynamic elements
			reset();

			// Render
			boid_render = givr::createInstancedRenderable(boid_geometry, boid_style);
//            wall_render = givr::createInstancedRenderable(wall_geometry, wall_style);
            sphere_render = givr::createInstancedRenderable(sphere_geometry, sphere_style);
		}

		void BoidsModel::reset() {
            grid.resize_grid(grid_w, grid_h, grid_d, glm::vec3(-0.5f * grid_w, -0.5f * grid_h, -0.5f * grid_d), cell_size);

            // static so they are persistent
			static std::random_device random_device;
			static std::mt19937 generator(random_device());
			static std::uniform_real_distribution<double> position_distribution(-0.5f * std::min(std::min(grid_w, grid_h), grid_d),0.5f * std::min(std::min(grid_w, grid_h), grid_d));
			static std::uniform_real_distribution<double> theta_distribution(-180., 180.);
			static std::uniform_real_distribution<double> phi_distribution(-90., 90.);
			static std::uniform_real_distribution<double> speed_distribution(5., 25.);

			boids.resize(n_boids);
            int boid_i = 0;
//#pragma omp parallel for private(boid_i) shared(boids)
            for (primatives::boid& boid : boids) {
//                int tid = omp_get_thread_num(); // Get thread ID
//                std::cout << "Thread " << tid << " is initializing boid " << boid_i << std::endl;


                //Random Position in 20 x 20 x 20 cube
                boid.index = boid_i;
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

                boid_i++;

                grid.add_boid(boid);
			}

		}

        void BoidsModel::step(float dt) {
            // Iterate through each boid
            grid.t_g++;
//#pragma omp parallel for shared(boids, grid)
            for (primatives::boid& bi : boids) {
                bi.p = clampPositionInGrid(bi.p);
                grid.add_boid(bi);
            }
//            std::cout << "2nd";
#pragma omp parallel for shared(boids, grid)
            for (primatives::boid& bi : boids) {
                std::vector<primatives::cell> n_i = grid.neighbourhood(bi.p, cell_size);

                bi.f = glm::vec3(0.f);
                bi.g = g;

                for (const primatives::cell& cell : n_i) {
                    for (int b_index : cell.boids) {
                        primatives::boid bj = boids[b_index];
                        if (&bi == &bj) continue; // Skip self-comparison
                        glm::vec3 DeltaXij = bj.p - bi.p;
                        float d = glm::length(DeltaXij);
                        float alpha = glm::dot(glm::normalize(DeltaXij), glm::normalize(bi.v));

                        if (d < r_s && alpha > cos(theta_s)) {
                            bi.f += separationForce(bi, bj);
                        } else if (d < r_a && alpha > cos(theta_a)) {
                            bi.f += alignmentForce(bi, bj);
                        } else if (d < r_c && alpha > cos(theta_c)) {
                            bi.f += cohesionForce(bi, bj);
                        }
                    }
                }

                for (const auto & plane : planes) {
                    bi.f += planeAvoidanceForce(bi, plane);
                }

                for (const auto & sphere : spheres) {
                    bi.f += sphereAvoidanceForce(bi, sphere);
                }

                // Integrate forward velocity and position for bi
                bi.integrate(dt);

                // Clamping the velocity
                float bi_speed = glm::length(bi.v);
                bi_speed = std::clamp(bi_speed, min_boid_v, max_boid_v);
                if (glm::length(bi.v) > 0.01f) {
                    bi.v = glm::normalize(bi.v) * bi_speed;
                }

                bi.orientate();
            }
        }


        void BoidsModel::render(const ModelViewContext& view) {
			//Add Mass render
			for (const primatives::boid& boid : boids)
                givr::addInstance(boid_render, calculateTransformMatrix(boid.p, boid.t, boid.n, boid.b));
//				givr::addInstance(boid_render, glm::translate(glm::mat4(1.f), boid.p)); //NEED TO FRAME!!!

            for (const primatives::sphere &sphere: spheres) {
                givr::addInstance(sphere_render, glm::scale(glm::translate(glm::mat4(1.f), sphere.origin), glm::vec3(sphere.radius)));
//                givr::addInstance(sphere_render, glm::translate(glm::scale(glm::mat4(1.f), glm::vec3(sphere.radius)), sphere.origin));
            }

			//Render
			givr::style::draw(boid_render, view);
            givr::style::draw(sphere_render, view);
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

        glm::vec3 BoidsModel::planeAvoidanceForce(const primatives::boid& bi, const primatives::plane& plane) const {
            glm::vec3 u = bi.p - plane.origin;
            glm::vec3 n_hat = glm::normalize(plane.normal);
            glm::vec3 v_hat = glm::normalize(bi.v);
            if (glm::dot(u, n_hat) < plane.epsilon) {
                float s_f = (glm::dot(u, n_hat) - plane.epsilon) * plane.wall_k_s * (-1.f);
                float d_f = glm::dot(bi.v, n_hat) * plane.wall_k_d;
                return (s_f + d_f) * n_hat;
            }
            else if (glm::dot(v_hat, n_hat) < 0.f && glm::dot(u, n_hat) < plane.collision_avoid_distance + plane.epsilon) {
                float dot_n_v = glm::dot(n_hat, v_hat);
                if (dot_n_v > 0.99f || dot_n_v < -0.99f) {
                    dot_n_v = std::clamp(dot_n_v, -0.99f, 0.99f);
                }
                glm::vec3 a_c_hat = glm::normalize(n_hat - dot_n_v * v_hat);
                float dot_a_n = glm::dot(a_c_hat, n_hat);
                if (dot_a_n > 0.99f) {
                    dot_a_n = 0.99f;
                }
                float r = (glm::dot(u, n_hat) - plane.epsilon) / (1 - dot_a_n);
//                if (std::isnan(r)) {
//                    std::cout << "WOOOW: " << glm::dot(a_c_hat, n_hat) << "\n";
//                }
                glm::vec3 a_c = (glm::dot(bi.v, bi.v) / r) * a_c_hat;
                return a_c;
            }
            return glm::vec3(0.f);
        }

        glm::vec3 BoidsModel::sphereAvoidanceForce(const primatives::boid& bi, const primatives::sphere& sphere) const {
            glm::vec3 u = bi.p - sphere.origin;
            glm::vec3 n_hat = glm::normalize(u);
            float u_length = glm::length(u);
            if (u_length <= sphere.radius + sphere.epsilon) {
                float s_f = (sphere.radius + sphere.epsilon - u_length) * sphere.wall_k_s;
                float d_f = glm::dot(bi.v, n_hat) * sphere.wall_k_d * (-1.f);
                return (s_f + d_f) * n_hat;
            }

            float v_length = glm::length(bi.v);
            if (v_length < 0.01f) {
                return glm::vec3(0.f);
            }

            glm::vec3 v_hat = glm::normalize(bi.v);
            float future_distance = std::sqrt((u_length * u_length) - (sphere.radius * sphere.radius));
            glm::vec3 future_pos = v_hat * future_distance;
//            if (glm::distance(future_pos, sphere.origin) <= sphere.radius + sphere.epsilon && u_length - sphere.radius - sphere.epsilon < sphere.collision_avoid_distance) {
            if (glm::dot(n_hat, bi.v) < 0.f && u_length - sphere.radius - sphere.epsilon < sphere.collision_avoid_distance) {
                float dot_n_v = glm::dot(n_hat, v_hat);
                if (dot_n_v > 0.99f || dot_n_v < -0.99f) {
                    dot_n_v = std::clamp(dot_n_v, -0.99f, 0.99f);
                }
                glm::vec3 a_c_hat = glm::normalize(n_hat - dot_n_v * v_hat);
                float dot_u_a = glm::dot(u, a_c_hat);
                if (dot_u_a > sphere.radius + sphere.epsilon - 0.01f) {
                    dot_u_a = sphere.radius + sphere.epsilon - 0.01f;
                }
                float r = ((u_length * u_length) - ((sphere.radius + sphere.epsilon) * (sphere.radius + sphere.epsilon))) / (2.f * (sphere.radius + sphere.epsilon - dot_u_a));
                glm::vec3 a_c = (glm::dot(bi.v, bi.v) / r) * a_c_hat;
//                std::cout << glm::length(a_c) << "\n";
                if (r == 0) {
                    std::cout << dot_u_a << "\n";
                }
                return a_c;
            }

            return glm::vec3(0.f);
        }

        glm::vec3 BoidsModel::clampPositionInGrid(const glm::vec3 p) const {
            glm::vec3 new_p = p;
            if (new_p.x >= grid.w + grid.origin.x) {
                new_p.x = grid.w + grid.origin.x - 0.01f;
            }
            else if (new_p.x <= grid.origin.x) {
                new_p.x = grid.origin.x + 0.01f;
            }

            if (new_p.y >= grid.h + grid.origin.y) {
                new_p.y = grid.h + grid.origin.y - 0.01f;
            }
            else if (new_p.y <= grid.origin.y) {
                new_p.y = grid.origin.y + 0.01f;
            }

            if (new_p.z >= grid.d + grid.origin.z) {
                new_p.z = grid.d + grid.origin.z - 0.01f;
            }
            else if (new_p.z <= grid.origin.z) {
                new_p.z = grid.origin.z + 0.01f;
            }

            return new_p;
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