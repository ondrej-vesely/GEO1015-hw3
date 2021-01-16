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
#include <sstream>
#include <cmath>
#include <iterator>
#include <algorithm>

#include "PlaneDetector.h"

/*
!!! TO BE COMPLETED !!!

Function that implements the RANSAC algorithm to detect one plane in the point cloud that is
accessible from _input_points variable (which contains the point read from the input ply file).

NOTICE that this function can be called multiple times! On each call a *new* plane should be generated
that does not contain any inliers that belong to a previously detected plane. Each plane should
have a unique segment_id (an int; '1' for the first plane, '2' for the second plane etc.).

Input:
	epsilon:          maximum distance from an inlier to corresponding plane.
	min_score:        minimum score (= number of inlier) of a plane. Planes with a lower score
					  should not be detected.
	k:                number of times a new minimal set and a consencus set is computed.

Output:
	Updated .segment_id's on the inliers in the newly detected plane. The inliers contain both the
	consensus set and minimal set.
*/
void PlaneDetector::detect_plane(double epsilon, int min_score, int k) {

	//-- TIP
	//-- access the input points from the _input_points:

	//   double x_of_point_i = _input_points[i].x;
	//   int segment_id_of_point_i = _input_points[i].segment_id;

	//-- set the segment_id of point i:

	//   _input_points[i].segment_id = 1;


	//-- TIP
	//-- Generating random numbers between 0 and 100

	// std::uniform_int_distribution<int> distrib(0, 100);
	// int my_random_number = distrib(_rand);

	//-- see https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution for more info



	// Indices of pts in best result
	std::vector<int> best = {};

	// Do k attempts to find best result
	for (int _k = 0; _k < k; _k++) {

		// Sample 3 unique points
		std::vector<Point> pts = sample(3);
		Point 
			pt1 = pts[0], 
			pt2 = pts[1], 
			pt3 = pts[2];

		// Get vectors of sampled plane
		double3 v1 = { pt2.x - pt1.x , pt2.y - pt1.y , pt2.z - pt1.z };
		double3 v2 = { pt3.x - pt1.x , pt3.y - pt1.y , pt3.z - pt1.z };
		double3 norm = normalize(cross(v1, v2));

		// Get params of sampled plane
		double
			A = norm[0],
			B = norm[1],
			C = norm[2],
			D = -A * pt1.x - B * pt1.y - C * pt1.z;

		// Check distance to plane to find inliers
		std::vector<int> inliers = {};

		for (int i = 0; i < _input_points.size(); i++) {
			Point p = _input_points[i];
			if (p.segment_id != 0) continue;

			double dist_pow = abs(A * p.x + B * p.y + C * p.z + D) / (A * A + B * B + C * C);
			if (dist_pow < epsilon * epsilon) {
				inliers.push_back(i);
			}
		}
		// Found new best result?
		if (inliers.size() > min_score && inliers.size() > best.size()) {
			best = inliers;
		}
	}

	// If the best better than minimal score, segment it
	if (best.size() > min_score) {
		_plane_count++;
		for (int i = 0; i < best.size(); i++) {
			_input_points[best[i]].segment_id = _plane_count;
		}
	}
}


/*
Function that calculates A,B,C,D params of plane defined by 3 points.
Input:
	pts:		Vector of 3 points.
Output:
	Vector of 4 floats defining plane as Ax+By+Cz+D=0.
*/
std::vector<float> _plane(std::vector<double3> pts) {

	double3
		pt1 = pts[0],
		pt2 = pts[1],
		pt3 = pts[2];

	double3 v1 = pt2 - pt1;
	double3 v2 = pt3 - pt1;
	double3 norm = normalize(cross(v1, v2));

	float
		A = norm[0],
		B = norm[1],
		C = norm[2],
		D = -A * pt1.x - B * pt1.y - C * pt1.z;

	std::vector<float> result = { A, B, C, D };
	return result;
}

/*
Function that samples n unique unsegmented points from the _input_points
Input: 
	n:		Number of samples to find
Output:
	std::vector of sampled <PlaneDetector::Point>s
*/
std::vector<PlaneDetector::Point> PlaneDetector::sample(int n) {

	std::vector<int> samples(n);
	std::vector<Point> result(n);
	std::uniform_int_distribution<int> distrib(0, _input_points.size()-1);

	if (n > _input_points.size()) return _input_points;

	for (int i = 0; i < n; i++) {
		int rand = distrib(_rand);
		while (std::count(samples.begin(), samples.end(), rand) 
			|| _input_points[rand].segment_id != 0) 
		{
			rand = distrib(_rand);
		}
		samples[i] = rand;
		result[i] = _input_points[rand];
	}
	return result;
}


// PLY I/O

/*
Function that writes the entire point cloud including the segment_id of each point to a .ply file

Input:
   filepath:  path of the .ply file to write the points with segment id

Template:
	ply
	format ascii 1.0
	element vertex 10000
	property float x
	property float y
	property float z
	property int segment_id
	end_header
	85176.77 446742.10 1.49 0
	85175.69 446742.58 1.52 1
	85174.45 446741.94 9.52 1
*/
void PlaneDetector::write_ply(std::string filepath) {

	std::ofstream file(filepath);
	char buffer[50];

	file << "ply" << "\n";
	file << "format ascii 1.0" << "\n";
	file << "element vertex " << _input_points.size() << "\n";
	file << "property float x" << "\n";
	file << "property float y" << "\n";
	file << "property float z" << "\n";
	file << "property int segment_id" << "\n";
	file << "end_header" << "\n";

	for (int i = 0; i < _input_points.size(); i++) {
		sprintf(buffer, "%f %f %f %d", 
			_input_points[i].x,
			_input_points[i].y,
			_input_points[i].z,
			_input_points[i].segment_id
		);
		file << buffer << "\n";
	}
	file.close();
}

/*
!!! DO NOT MODIFY read_ply() !!!

This function is already implemented.
*/
bool PlaneDetector::read_ply(std::string filepath) {

	std::cout << "Reading file: " << filepath << std::endl;
	std::ifstream infile(filepath.c_str(), std::ifstream::in);
	if (!infile)
	{
		std::cerr << "Input file not found.\n";
		return false;
	}
	std::string cursor;

	// start reading the header
	std::getline(infile, cursor);
	if (cursor != "ply") {
		std::cerr << "Magic ply keyword not found\n";
		return false;
	};

	std::getline(infile, cursor);
	if (cursor != "format ascii 1.0") {
		std::cerr << "Incorrect ply format\n";
		return false;
	};

	// read the remainder of the header
	std::string line = "";
	int vertex_count = 0;
	bool expectVertexProp = false, foundNonVertexElement = false;
	int property_count = 0;
	std::vector<std::string> property_names;
	int pos_x = -1, pos_y = -1, pos_z = -1, pos_segment_id = -1;

	while (line != "end_header") {
		std::getline(infile, line);
		std::istringstream linestream(line);

		linestream >> cursor;

		// read vertex element and properties
		if (cursor == "element") {
			linestream >> cursor;
			if (cursor == "vertex") {
				// check if this is the first element defined in the file. If not exit the function
				if (foundNonVertexElement) {
					std::cerr << "vertex element is not the first element\n";
					return false;
				};

				linestream >> vertex_count;
				expectVertexProp = true;
			}
			else {
				foundNonVertexElement = true;
			}
		}
		else if (expectVertexProp) {
			if (cursor != "property") {
				expectVertexProp = false;
			}
			else {
				// read property type
				linestream >> cursor;
				if (cursor.find("float") != std::string::npos || cursor == "double") {
					// read property name
					linestream >> cursor;
					if (cursor == "x") {
						pos_x = property_count;
					}
					else if (cursor == "y") {
						pos_y = property_count;
					}
					else if (cursor == "z") {
						pos_z = property_count;
					}
					++property_count;
				}
				else if (cursor.find("uint") != std::string::npos || cursor == "int") {
					// read property name
					linestream >> cursor;
					if (cursor == "segment_id") {
						pos_segment_id = property_count;
					}
					++property_count;
				}
			}

		}
	}

	// check if we were able to locate all the coordinate properties
	if (pos_x == -1 || pos_y == -1 || pos_z == -1) {
		std::cerr << "Unable to locate x, y and z vertex property positions\n";
		return false;
	};

	// read the vertex properties
	for (int vi = 0; vi < vertex_count; ++vi) {
		std::getline(infile, line);
		std::istringstream linestream(line);

		double x{}, y{}, z{};
		int sid{};
		for (int pi = 0; pi < property_count; ++pi) {
			linestream >> cursor;
			if (pi == pos_x) {
				x = std::stod(cursor);
			}
			else if (pi == pos_y) {
				y = std::stod(cursor);
			}
			else if (pi == pos_z) {
				z = std::stod(cursor);
			}
			else if (pi == pos_segment_id) {
				sid = std::stoi(cursor);
			}
		}
		auto p = Point{ x, y, z };
		if (pos_segment_id != -1) {
			p.segment_id = sid;
		}
		_input_points.push_back(p);
	}

	std::cout << "Number of points read from .ply file: " << _input_points.size() << std::endl;

	return true;
}