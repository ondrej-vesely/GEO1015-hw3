/*
  GEO1015.2020
  hw03
  --
  [YOUR NAME]
  [YOUR STUDENT NUMBER]
  [YOUR NAME]
  [YOUR STUDENT NUMBER]
*/


#include <iostream>
#include <fstream>
#include <chrono>
#include <filesystem>
namespace fs = std::filesystem;

//-- to read the params.json file
#include <nlohmann-json/json.hpp>
using json = nlohmann::json;

//-- our PlaneDetector class
#include "PlaneDetector.h"

//-- declarations of the runViewer function that launches the Viewer
int runViewer(PlaneDetector& detector, int argc, char** argv);


/*
!!! DO NOT MODIFY main() !!!
*/
int main(int argc, char** argv)
{
	//-- setup path to params.json file
	std::string json_path = JSON_PARAMS_PATH;
	//-- take the json path from command line if supplied so
	if (argc == 2) {
		json_path = argv[1];
	}
	
	//-- set current working directory to the parent directory of json_path. As a result the
	//-- filepaths in the json file will be read relative to the location of the json file
	fs::path working_directory = fs::path(json_path).parent_path();
	fs::current_path(working_directory);
	std::cout << "Active working directory: " << working_directory << std::endl;

	//-- read the params.json file
	std::ifstream json_file(json_path);
	if (!json_file) {
		std::cerr << "JSON file " << json_path << " not found.\n";
		return 1;
	}
	json j; json_file >> j;
	int n_planes = j["n_planes"];
	int k = j["k"];
	int min_score = j["min_score"];
	float epsilon = j["epsilon"];
	std::string input_file = j["input_file"];
	std::string output_file = j["output_file"];

	PlaneDetector detector;

	//-- read point cloud from input .ply file, exit if it fails
	if (!detector.read_ply(input_file)) {
		return 1;
	}

	//-- perform plane detection and time how long it takes
  auto start = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < n_planes; ++i) {
		detector.detect_plane(epsilon, min_score, k);
	}

  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << std::fixed << std::setprecision(3) << "--- Plane detection took " << elapsed.count() << " seconds ---" << std::endl;


	//-- open the viewer
	runViewer(detector, argc, argv);

	//-- write the detection result to output .ply file (after the viewer is closed)
	detector.write_ply(output_file);

	//-- we're done, return 0 to say all went fine
	return 0;
}
