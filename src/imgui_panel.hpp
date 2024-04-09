#pragma once

#include <givio.h>
#include <givr.h>
#include <imgui/imgui.h>

namespace imgui_panel {
	extern bool showPanel;
	extern ImVec4 clear_color;
	extern bool reset_view;

	//Simulation settings
	extern int number_of_iterations_per_frame;
	extern bool play_simulation;
	extern bool reset_simulation;
	extern bool step_simulation;
	extern float dt_simulation;

    extern int number_of_boids;

    extern float k_s;
    extern float k_a;
    extern float k_c;

    extern float r_s;
    extern float r_a;
    extern float r_c;

    extern float th_s;
    extern float th_a;
    extern float th_c;

	// lambda function
	extern std::function<void(void)> draw;
} // namespace panel