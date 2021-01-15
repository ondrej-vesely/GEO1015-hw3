/*
  GEO1015.2020
  hw03 
  --
  [YOUR NAME] 
  [YOUR STUDENT NUMBER] 
  [YOUR NAME] 
  [YOUR STUDENT NUMBER] 
*/


#include <string>
#include <vector>
#include <random>

//-- Simple linear algebra library that you can use for distance computations etc
//-- See https://github.com/sgorsten/linalg
#include <linalg/linalg.h>
using double3 = linalg::aliases::double3;

/*
The main class of this assignment. You are allowed to add functions and member variables for your own use.
*/
class PlaneDetector {


  //-- you can add your own variables and functions here


  /*
  !!! DO NOT MODIFY below this line !!!
  */
  public:

  /*
  We define a Point struct that inherits from the double3 struct that is defined in linalg.h. 
  This means you can use a Point as a double3, ie. with all the linear algebra functions defined for double3 in linalg.h.
  The only thing added here is the segment_id, which we use to indicate to what plane this point is assigned. 

  NOTICE that the segment_id==0 (the default value) means that the point is not assigned to any plane.
  */
  struct Point : double3 {
    using double3::double3;
    int segment_id{0};
  };

  //-- The main plane detection function where you need to implement the RANSAC algorithm (in the PlaneDetector.cpp file)
  void detect_plane(double epsilon, int min_score, int k);

  //-- .PLY reading (already implemented for you)
  bool read_ply(std::string filepath);
  //-- .PLY writing (you need to implement in PlaneDetector.cpp)
  void write_ply(std::string filepath);
  
  //-- point cloud access (for the viewer application)
  const std::vector<Point>& get_input_points() {
    return _input_points;
  };

  private:

  //-- This variable holds the entire input point cloud after calling read_ply()
  std::vector<Point> _input_points;

  //-- random number generator to generate your random numbers (important for RANSAC!)
  std::mt19937 _rand{ std::mt19937{std::random_device{}()} };

};