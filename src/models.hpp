#pragma once

#include <vector>
#include <givr.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/compatibility.hpp> // lerp

namespace simulation {
	namespace primatives {
		//Boids in simulation
		struct boid {
			glm::vec3 p = glm::vec3(0.f);
			glm::vec3 v = glm::vec3(0.f);
            float m = 1.f;
            glm::vec3 f = glm::vec3(0.f);
			//Note: Mass is implicitly 1 (our choice) so Force = Acceleration
			//TO-DO: Modify this class to include certain desired quantities (mass, force, ...)
			//May even add functions! Such as integration ...

            // Integration function
            void integrate(float dt) {
                glm::vec3 a = f / m;
                v += a * dt;
                p += v * dt;
            }
        };

		//Walls for avoidences and simulation boundry
		struct plane {
			//TO-DO: Define a splne for collision avoidence and container purposes
		};

		//Spheres for avoidences
		struct sphere {
			//TO-DO: Define a sphere for collision avoidence purposes
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

			//Simulation Constants (you can re-assign values here from imgui)
			glm::vec3 g = { 0.f, -9.81f, 0.f };
			size_t n_boids = 100; //need alot more eventually for full assignment

		private:
			//Simulation Parts
			std::vector<primatives::boid> boids;
			//std::vector<primatives::plane> planes;
			//std::vector<primatives::sphere> spheres;

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