#include "imgui_panel.hpp"

namespace imgui_panel {
	// default values
	bool showPanel = true;
	ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
	bool reset_view = false;
	

	//Simulation settings
	int number_of_iterations_per_frame = 1;
	bool play_simulation = false;
	bool reset_simulation = false;
	bool step_simulation = false;
	float dt_simulation = 0.01f;

    int number_of_boids = 500;

    float k_s = 15.f;
    float k_a = 0.4f;
    float k_c = 0.2f;

    float r_s = 6.f;
    float r_a = 8.f;
    float r_c = 10.f;

    float th_s = 175.f;
    float th_a = 140.f;
    float th_c = 120.f;

    bool show_generic = true;
    bool show_k_values = false;
    bool show_radii = false;
    bool show_view_angles = false;

	std::function<void(void)> draw = [](void) {
		if (showPanel && ImGui::Begin("Panel", &showPanel, ImGuiWindowFlags_MenuBar)) {
			ImGui::Spacing();
			ImGui::Separator();

			ImGui::ColorEdit3("Clear color", (float*)&clear_color);
			reset_view = ImGui::Button("Reset View");

			ImGui::Spacing();
			ImGui::Separator();

			ImGui::SliderInt("Iterations Per Frame", &number_of_iterations_per_frame, 1, 100);
			ImGui::Checkbox("Play Simulation", &play_simulation);
			reset_simulation = ImGui::Button("Reset Simulation");
			if (!play_simulation) {
				step_simulation = ImGui::Button("Step Simulation");
			}
			ImGui::DragFloat("Simulation dt", &dt_simulation, 1.e-5f, 1.e-5f, 1.f, "%.6e");

			ImGui::Separator();

			float frame_rate = ImGui::GetIO().Framerate;
			ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
				1000.0f / frame_rate, frame_rate);

			ImGui::Spacing();
			ImGui::Separator();

            if (ImGui::BeginTabBar("1")) {
                if (ImGui::TabItemButton("Generic")) {
                    show_generic = true;
                    show_k_values = false;
                    show_radii = false;
                    show_view_angles = false;
                }

                if (show_generic) {
                    ImGui::DragInt("Number of Boids", &number_of_boids, 5, 1, 5000);
                }

                if (ImGui::TabItemButton("K Values")) {
                    show_generic = false;
                    show_k_values = true;
                    show_radii = false;
                    show_view_angles = false;
                }

                if (show_k_values) {
                    ImGui::DragFloat("Ks", &k_s, 0.05f, 0.f, 100.f);
                    ImGui::DragFloat("Ka", &k_a, 0.05f, 0.f, 100.f);
                    ImGui::DragFloat("Kc", &k_c, 0.05f, 0.f, 100.f);
                }

                if (ImGui::TabItemButton("Radii")) {
                    show_generic = false;
                    show_k_values = false;
                    show_radii = true;
                    show_view_angles = false;
                }

                if (show_radii) {
                    ImGui::DragFloat("Rs", &r_s, 0.5f, 0.f, 100.f);
                    ImGui::DragFloat("Ra", &r_a, 0.05f, 0.f, 100.f);
                    ImGui::DragFloat("Rc", &r_c, 0.05f, 0.f, 100.f);
                }

                if (ImGui::TabItemButton("View Angles")) {
                    show_generic = false;
                    show_k_values = false;
                    show_radii = false;
                    show_view_angles = true;
                }

                if (show_view_angles) {
                    ImGui::DragFloat("Theta_s", &th_s, 1.f, 0.f, 360.f);
                    ImGui::DragFloat("Theta_a", &th_a, 1.f, 0.f, 360.f);
                    ImGui::DragFloat("Theta_c", &th_c, 1.f, 0.f, 360.f);
                }

                ImGui::EndTabBar();
            }

            ImGui::End();
		}
	};
} // namespace panel