#pragma once

#include <vector>
#include <givr.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <cmath>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/compatibility.hpp> // lerp

namespace simulation {
	namespace primatives {
		//Boids in simulation
		struct boid {
			glm::vec3 p = glm::vec3(0.f);
			glm::vec3 v = glm::vec3(0.f);
            glm::vec3 f = glm::vec3(0.f);
            glm::vec3 g = glm::vec3(0.f);

            // for orientation
            glm::vec3 t = glm::vec3(1.f, 0.f, 0.f);
            glm::vec3 b = glm::vec3(0.f, 0.f, 1.f);
            glm::vec3 n = glm::vec3(0.f, 1.f, 0.f);
			//Note: Mass is implicitly 1 (our choice) so Force = Acceleration
			//TO-DO: Modify this class to include certain desired quantities (mass, force, ...)
			//May even add functions! Such as integration ...

            // Integration function
            void integrate(float dt) {
                glm::vec3 a = f + g;
                v += a * dt;
                p += v * dt;
//                std::cout << (v).x << ", " << (v).y << ", " << (v).z << "\n";
            }

            void orientate() {
                t = glm::normalize(v);
                b = glm::normalize(glm::cross(t, f - g));
                n = glm::normalize(glm::cross(b, t));
            }
        };

		//Walls for avoidences and simulation boundry
		struct plane {
			//TO-DO: Define a splne for collision avoidence and container purposes
            glm::vec3 origin = glm::vec3(0.f);
            glm::vec3 normal = glm::vec3(0.f, 1.f, 0.f);
            float epsilon = 5.f;
            float collision_avoid_distance = 200.f;
            float wall_k_s = 100.f, wall_k_d = 0.2f;
		};

		//Spheres for avoidences
		struct sphere {
			//TO-DO: Define a sphere for collision avoidence purposes
            glm::vec3 origin = glm::vec3(0.f);
            float radius = 1.f;
		};
	} // namespace primatives

	namespace models {
		//If you want to use a different view, change this and the one in main
		using ModelViewContext = givr::camera::ViewContext<givr::camera::TurnTableCamera, givr::camera::PerspectiveProjection>;
		// Abstract class used by all models
		class GenericModel {
		public:
			virtual void reset() = 0;
			virtual void step(float dt) = 0;
			virtual void render(const ModelViewContext& view) = 0;
		};

		//Model constructing a single spring
		class BoidsModel : public GenericModel {
		public:
			BoidsModel();
			void reset();
			void step(float dt);
			void render(const ModelViewContext& view);
            glm::vec3 separationForce(const primatives::boid& bi, const primatives::boid& bj) const;
            glm::vec3 cohesionForce(const primatives::boid& bi, const primatives::boid& bj) const;
            glm::vec3 alignmentForce(const primatives::boid& bi, const primatives::boid& bj) const;
            glm::vec3 planeAvoidanceForce(const primatives::boid& bi, const primatives::plane& plane) const;
            glm::vec3 sphereAvoidanceForce(const primatives::boid& bi, const primatives::sphere& sphere) const;
            static glm::mat4 calculateTransformMatrix(glm::vec3 position, glm::vec3 tangent, glm::vec3 normal, glm::vec3 binormal);

			//Simulation Constants (you can re-assign values here from imgui)
			glm::vec3 g = { 0.f, -2.81f, 0.f };
			size_t n_boids = 100; //need alot more eventually for full assignment
            float r_s = 6.f, r_a = 8.f, r_c = 10.f;
            float theta_s = 175.f * M_PI / 180.f, theta_a = 140.f * M_PI / 180.f, theta_c = 120.f * M_PI / 180.f;
            float k_s = 15.f, k_a = 0.4f, k_c = 0.2f;
            float min_boid_v = 1.f, max_boid_v = 25.f;

            //Collisions
//            std::vector<primatives::plane> planes;
//            std::vector<primatives::sphere> spheres;

		private:
			//Simulation Parts
			std::vector<primatives::boid> boids;
			std::vector<primatives::plane> planes;
			std::vector<primatives::sphere> spheres;

			//Render
			givr::geometry::Mesh boid_geometry;
			givr::style::Phong boid_style;
			givr::InstancedRenderContext<givr::geometry::Mesh, givr::style::Phong> boid_render;

			//givr::geometry::TriangleSoup wall_geometry;
			//givr::style::Phong wall_style;
			// Maybe changed this (below) to instanced render????????
			//givr::RenderContext<givr::geometry::TriangleSoup, givr::style::Phong> wall_render;

			//givr::geometry::Sphere sphere_geometry;
			//givr::style::Phong sphere_style;
			//givr::InstancedRenderContext<givr::geometry::TriangleSoup, givr::style::Phong> sphere_render;
		};
	} // namespace models
} // namespace simulation