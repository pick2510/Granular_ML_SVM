#ifndef DISPLAYSPHERES_INCLUDED
#define DISPLAYSPHERES_INCLUDED

#include <vector>
#include <functional>
#include <string>

//default: false. If set to true, DisplaySpheres() and DisplaySpheres_Movie() will not actually display anything, but will just output 1 or more povray scripts and exit.
extern bool SavePovrayAndExit;
extern std::string SavePovray_Prefix;
extern bool OrthoProjection;
extern long DisplaySpheresMovie_CurrentFrame;

struct sphere
{
	double x, y, z;//coordinates, on the order of 1
	               //the rotation center is (0, 0, 0) 
	double radius;

	double red, green, blue;//color, from 0 to 1, 0:dark, 1:bright
	double transparency;//from 0 to 1, 0:invisible, 1:opaque
};
struct line
{
	double x1, y1, z1;
	double x2, y2, z2;

	double width;// on the order of 1

	double red, green, blue;//color, from 0 to 1, 0:dark, 1:bright
	double transparency;//from 0 to 1, 0:invisible, 1:opaque

	//if the line has an arrow head at the end
	//currently, only PovRay output supports this
	//OpenGL output simply ignores this.
	bool arrowHead;
	line() : arrowHead(false) {}
};

//should only be called inside the Extra2DStuff parameters of DisplaySpheres functions.
void DisplayText(const std::string & text, double x=0.0, double y=0.0);

// Display all spheres in the vector
// Use mouse to drag spheres
// Use keyboard arrows to move the camera
// Use keys '+' and '-' to increase or decrease the precision. Spheres looks bad if precision is too low. Computer is laggy if precision is too high.
// Use key 'q' to return from this function
// Use key 'o' to switch to orthogonal projection
// Use key 'p' to switch to perspective projection (default)
// Use key 'v' to write a Pov-ray source file of current configuration and view angle
// Don't click the Close Window icon in the upper right corner of the title bar. That will kill the program.
// In order for this function to work, define USE_OPENGL when compiling, and link with OpenGL, GLU, and GLUT
void DisplaySpheres(const std::vector<sphere> & VSpheres, const std::vector<line> & VLines, std::function<void(size_t, size_t)> Extra2DStuff=[](size_t, size_t)->void{});


//Similar to DisplaySpheres, except renders a movie
// Use key 'a' to accelerate the movie
// Use key 's' to deccelerate the movie
// Use key 'd' to toggle pause/play (pause at the beginning)
// Use key 'r' to toggle forward/backward (forward by default)
void DisplaySpheres_Movie(std::function<void(std::vector<sphere> & VSpheres, std::vector<line> & VLines, long FrameNumber)> GetFrameFunc, long NumFrames);

#endif