/*
  GEO1015.2020
  hw03
  --
  Ondrej Vesely
  5162130
  Guilherme Spinoza Andreo
  5383994
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

	//-- I extend params.json with some addtional parameters
	//-- If they are missing in the file, defaults defined in PlaneDetector are used.
	try {
		detector.chunk_size = j["chunk_size"];
		detector.chunk_extrapolate = j["chunk_extrapolate"];
		detector.check_normals = j["check_normals"];
		detector.normal_tol = j["normal_tol"];
		detector.normal_radius = j["normal_radius"];
	}
	catch (...) {
		std::cout << "WARN: Your params.json seems to be missing some values.\n"
				  << "		See params below for reference.\n";
	}
	std::cout
		<< "PARAMS: " << "\n"
		<< "	n_planes: " << n_planes << "\n"
		<< "	k: " << k << "\n"
		<< "	min_score: " << min_score << "\n"
		<< "	epsilon: " << epsilon << "\n"
		<< "	chunk_size: " << detector.chunk_size << "\n"
		<< "	chunk_extrapolate: " << std::boolalpha << detector.chunk_extrapolate << "\n"
		<< "	check_normals: " << std::boolalpha << detector.check_normals << "\n"
		<< "	normal_tol: " << detector.normal_tol << "\n"
		<< "	normal_radius: " << detector.normal_radius << "\n";

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
