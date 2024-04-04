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
	float dt_simulation = 0.001f;

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

			ImGui::End();
		}
	};
} // namespace panel