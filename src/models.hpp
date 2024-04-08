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
            int index;

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
            float epsilon = 2.f;
            float collision_avoid_distance = 20.f;
            float wall_k_s = 100.f, wall_k_d = 0.2f;
		};

		//Spheres for avoidences
		struct sphere {
			//TO-DO: Define a sphere for collision avoidence purposes
            glm::vec3 origin = glm::vec3(0.f);
            float radius = 10.f;
            float epsilon = 1.f;
            float collision_avoid_distance = 10.f;
            float wall_k_s = 100.f, wall_k_d = 0.2f;
		};

        struct cell {
            std::vector<int> boids;
            size_t t_c = 0;

            void clear_boids() {
                boids.clear();
            }

            void add_boid(int boid_i) {
                boids.push_back(boid_i);
//                std::cout << boids.size() << "\n";
            }
        };

        struct grid {
            std::vector<cell> cells;
            size_t t_g = 0;
            float w, h, d;
            glm::vec3 origin = glm::vec3(0.f);
            float cell_size = 1.f;
            int n_w, n_h, n_d;

            void init_cells() {
                int cell_count = n_w * n_h * n_d;
                cells.clear();
                cells.resize(cell_count);

                for (int i = 0; i < cell_count; ++i) {
                    cells[i] = cell();
                    cells[i].clear_boids();
                    cells[i].t_c = t_g;
                }
            }

            cell& get_cell(int i_w, int i_h, int i_d) {
                int i_a = i_w + n_w * i_h + n_w * n_h * i_d;
                cell& c = cells[i_a];
                if (t_g != c.t_c) {
                    c.clear_boids();
                    c.t_c = t_g;
                }
                return c;
            }

            cell& get_cell(glm::vec3 p) {
                int i_w = std::floor((p.x - origin.x) / cell_size);
                int i_h = std::floor((p.y - origin.y) / cell_size);
                int i_d = std::floor((p.z - origin.z) / cell_size);
                if (i_w < 0 || i_w >= n_w || i_h < 0 || i_h >= n_h || i_d < 0 || i_d >= n_d) {
                    std::cout << "WTF: " << p.x << ", " << p.y << ", " << p.z << "\n";
                }
                return get_cell(i_w, i_h, i_d);
            }

            void resize_grid(float new_w, float new_h, float new_d, glm::vec3 new_origin, float new_cell_size) {
                w = new_w;
                h = new_h;
                d = new_d;
                cell_size = new_cell_size;
                origin = new_origin;
                n_w = std::ceil(w / cell_size);
                n_h = std::ceil(h / cell_size);
                n_d = std::ceil(d / cell_size);

                init_cells();
            }

            void add_boid(boid& b) {
                cell& c = get_cell(b.p);
                c.add_boid(b.index);
//                std::cout << c.boids.size() << "\n";
            }

            std::vector<cell> neighbourhood(glm::vec3 p, float r) {
                std::vector<cell> neighborhoodCells;

                // Calculate the indices of the cell containing the position 'p'
                int i_w = std::floor((p.x - origin.x) / cell_size);
                int i_h = std::floor((p.y - origin.y) / cell_size);
                int i_d = std::floor((p.z - origin.z) / cell_size);

                // Add neighboring cells to the result
                int search_size = std::ceil(r / cell_size);
                for (int dx = -1 * search_size; dx <= search_size; ++dx) {
                    for (int dy = -1 * search_size; dy <= search_size; ++dy) {
                        for (int dz = -1 * search_size; dz <= search_size; ++dz) {
                            int neighbor_i_w = i_w + dx;
                            int neighbor_i_h = i_h + dy;
                            int neighbor_i_d = i_d + dz;

                            if (neighbor_i_w >= 0 && neighbor_i_w < n_w &&
                                neighbor_i_h >= 0 && neighbor_i_h < n_h &&
                                neighbor_i_d >= 0 && neighbor_i_d < n_d) {
                                neighborhoodCells.push_back(get_cell(neighbor_i_w, neighbor_i_h, neighbor_i_d));
//                                if (!get_cell(neighbor_i_w, neighbor_i_h, neighbor_i_d).boids.empty()) {
//                                    std::cout << "yay?" << "\n";
//                                }
                            }
                        }
                    }
                }
//                std::cout << neighborhoodCells.size() << "\n";
                return neighborhoodCells;
            }
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
            glm::vec3 clampPositionInGrid(const glm::vec3 p) const;
            static glm::mat4 calculateTransformMatrix(glm::vec3 position, glm::vec3 tangent, glm::vec3 normal, glm::vec3 binormal);

			//Simulation Constants (you can re-assign values here from imgui)
			glm::vec3 g = { 0.f, -2.81f, 0.f };
			size_t n_boids = 500; //need alot more eventually for full assignment
            float r_s = 5.f, r_a = 6.f, r_c = 8.f;
            float theta_s = 175.f * M_PI / 180.f, theta_a = 140.f * M_PI / 180.f, theta_c = 120.f * M_PI / 180.f;
            float k_s = 15.f, k_a = 0.4f, k_c = 0.2f;
            float min_boid_v = 10.f, max_boid_v = 25.f;

            float grid_w = 150.f, grid_h = 150.f, grid_d = 150.f;
            float cell_size = std::max(r_s, std::max(r_c, r_a));

		private:
			//Simulation Parts
			std::vector<primatives::boid> boids;
			std::vector<primatives::plane> planes;
			std::vector<primatives::sphere> spheres;

            primatives::grid grid;

			//Render
			givr::geometry::Mesh boid_geometry;
			givr::style::Phong boid_style;
			givr::InstancedRenderContext<givr::geometry::Mesh, givr::style::Phong> boid_render;

			givr::geometry::TriangleSoup wall_geometry;
			givr::style::Phong wall_style;
			// Maybe changed this (below) to instanced render????????
			givr::RenderContext<givr::geometry::TriangleSoup, givr::style::Phong> wall_render;

			givr::geometry::Sphere sphere_geometry;
			givr::style::Phong sphere_style;
			givr::InstancedRenderContext<givr::geometry::Sphere, givr::style::Phong> sphere_render;
		};
	} // namespace models
} // namespace simulation