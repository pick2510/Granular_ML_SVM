// DisplaySpheres.cpp
// adapted from SphereWorld.cpp found in OpenGL SuperBible, by Richard S. Wright Jr.

#include "DisplaySpheres.h"

bool SavePovrayAndExit = false;
std::string SavePovray_Prefix = "";

#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <algorithm>
long DisplaySpheresMovie_CurrentFrame;
bool OrthoProjection=false;


#ifdef USE_OPENGL
#include "etc.h"
#include "math3d.h"
#include "glframe.h"
//#include <GL/glut.h>
#include <GL/freeglut.h>



static std::vector<GLFrame> spheres;
static GLFrame    frameCamera;
static GLFrame		AllSpheres;

static std::vector<sphere> Spheres;
static std::vector<size_t> RenderOrder;
static std::vector<line> Lines;
static std::function<void(std::vector<sphere> & VSpheres, std::vector<line> & VLines, long FrameNumber)> * pGetFrameFunc;
static std::function<void(size_t, size_t)> * pExtra2DStuffFunc;
// Light and material Data
static GLfloat fLightPos[4]   = { -100.0f, 100.0f, 100.0f, 1.0f };  // Point source
static GLfloat fNoLight[] = { 0.0f, 0.0f, 0.0f, 0.0f };
static GLfloat fAmbientLight[] = { 0.19225f, 0.19225f, 0.19225f, 1.0f };
static GLfloat fDiffuseLight[] = { 0.50754f, 0.50754f, 0.50754f, 1.0f };
static GLfloat fSpecularLight[] = { 0.508273f, 0.508273f, 0.508273f, 1.0f };
static GLfloat ClearColor[] = { 1.0f, 1.0f, 1.0f, 1.0f };
static GLfloat shine=0.4;

static size_t Precision = 16;
static GLfloat fAspect;
static size_t WindowWidth=800;
static size_t WindowHeight=600;

static long FrameRate = 1000000; //0.1 second in Windows
static bool Playing = false;
static bool Reverse = false;
static long MaxFrame;
static long LastFrameTime;

static int window_handle;

void PrepareFrame (long FrameNumber)
{
	(*pGetFrameFunc)(::Spheres, ::Lines, FrameNumber);
	DisplaySpheresMovie_CurrentFrame=FrameNumber;
	LastFrameTime=GetPreciseClock();

	spheres.resize(Spheres.size());
	for(int iSphere = 0; iSphere < Spheres.size(); iSphere++)
	{
		spheres[iSphere].SetOrigin(static_cast<float>(Spheres.at(iSphere).x), static_cast<float>(Spheres.at(iSphere).y), static_cast<float>(Spheres.at(iSphere).z));
	}
	RenderOrder.clear();
	for(size_t i=0; i<Spheres.size(); i++)
		RenderOrder.push_back(i);
}

//////////////////////////////////////////////////////////////////
// This function does any needed initialization on the rendering
// context. 

void PlaceItems()
{
	// place the spheres
	AllSpheres.SetOrigin(0.0, 0.0, -0.0);
	frameCamera.SetOrigin(0.0, 0.0, 3.0);
}
void SetupRC()
{
	// Grayish background
	glClearColor(ClearColor[0], ClearColor[1], ClearColor[2], ClearColor[3]);

	// Cull backs of polygons
	//glCullFace(GL_BACK);
	//glFrontFace(GL_CCW);
	//glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);

	// Setup light parameters
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, fNoLight);
	glLightfv(GL_LIGHT0, GL_AMBIENT, fAmbientLight);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, fDiffuseLight);
	glLightfv(GL_LIGHT0, GL_SPECULAR, fSpecularLight);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	//glEnable(GL_POLYGON_SMOOTH);

	// Mostly use material tracking
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glMaterialfv(GL_FRONT, GL_SPECULAR, fSpecularLight);
	glMateriali(GL_FRONT, GL_SHININESS, shine*128);

	PlaceItems();

	PrepareFrame(DisplaySpheresMovie_CurrentFrame);
}

void DisplayText(const std::string & text, double x, double y)
{
	glColor3f(0.0, 0.0, 0.0);
	glRasterPos2f(x, y);
	for ( int i = 0; i < text.size(); ++i ) 
	{
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, text[i]);
	}
}

//save the figure as a povray script
void SavePovRay(std::string name="")
{
	auto t = std::time(nullptr);

	struct tm * now = localtime( & t );
	std::stringstream ss;
	if (name == "")
	{
		ss  << SavePovray_Prefix << '_'
			<< (now->tm_year + 1900) << '_'
			<< (now->tm_mon + 1) << '_'
			<< now->tm_mday << '_'
			<< now->tm_hour << '_'
			<< now->tm_min << '_'
			<< now->tm_sec
			<< ".pov";
	}
	else
		ss << name << ".pov";
	std::fstream ofile(ss.str(), std::fstream::out);
	M3DVector3f co, ca;
	ofile<<"#version 3.6;\n";

	if (OrthoProjection)
	{
		ofile << "camera { orthographic \n right 1.00*x ";
		frameCamera.GetOrigin(co);
		ofile << "  location <" << (-1)*(co[0]) << "," << co[1] << "," << co[2] << ">";
		frameCamera.GetForwardVector(ca);
		ofile << "  look_at <" << (-1)*(co[0] + ca[0]) << "," << co[1] + ca[1] << "," << co[2] + ca[2] << ">";
		ofile << "  sky y";
		ofile << " angle 80}\n";
	}
	else
	{
		ofile << "camera { right 1.00*x ";
		frameCamera.GetOrigin(co);
		ofile << "  location <" << (-1)*(co[0]) << "," << co[1] << "," << co[2] << ">";
		frameCamera.GetForwardVector(ca);
		ofile << "  look_at <" << (-1)*(co[0] + ca[0]) << "," << co[1] + ca[1] << "," << co[2] + ca[2] << ">";
		ofile << "  sky y";
		ofile << " angle 29}\n";
	}
	ofile<<"background{rgb 1}\n";
	ofile<<"\n";
	ofile<<"#declare REFLECTION = 0.02;\n#declare SPECULAR = 0.5;\n#declare AMBIENT = 0.48;\n#declare DIFFUSE = 0.8;\n";
	struct rgb
	{
		double r, g, b, t;
		rgb(const sphere & s) : r(s.red), g(s.green), b(s.blue), t(1.0 - s.transparency)
		{
		}
		rgb(const line & s) : r(s.red), g(s.green), b(s.blue), t(1.0 - s.transparency)
		{
		}
		bool operator < (const rgb & right) const
		{
			if (this->r < right.r)
				return true;
			else if (this->r>right.r)
				return false;
			else if (this->g < right.g)
				return true;
			else if (this->g>right.g)
				return false;
			else if (this->b < right.b)
				return true;
			else if (this->b>right.b)
				return false;
			else if (this->t < right.t)
				return true;
			else 
				return false;
		}
	};
	std::map<rgb, int> rgbMap;
	int index = 0;
	int index2 = 0;
	std::map<double, int> radiusMap;
	std::stringstream rgbString, radiusString;
	for (int i = 0; i < Spheres.size(); i++)
	{
		rgb c(Spheres[i]);
		if (rgbMap.find(c) == rgbMap.end())
		{
			index++;
			rgbMap.insert(std::make_pair(c, index));
			rgbString << "#declare r" << index << " = " << c.r << ";\n";
			rgbString << "#declare g" << index << " = " << c.g << ";\n";
			rgbString << "#declare b" << index << " = " << c.b << ";\n";
			rgbString << "#declare t" << index << " = " << c.t << ";\n";
		}
		double radius = Spheres[i].radius;
		if (radiusMap.find(radius) == radiusMap.end())
		{
			index2++;
			radiusMap.insert(std::make_pair(radius, index2));
			radiusString << "#declare Radius" << index2 << " = " << radius << ";\n";
		}
	}
	for (int i = 0; i < Lines.size(); i++)
	{
		rgb c(Lines[i]);
		if (rgbMap.find(c) == rgbMap.end())
		{
			index++;
			rgbMap.insert(std::make_pair(c, index));
			rgbString << "#declare r" << index << " = " << c.r << ";\n";
			rgbString << "#declare g" << index << " = " << c.g << ";\n";
			rgbString << "#declare b" << index << " = " << c.b << ";\n";
			rgbString << "#declare t" << index << " = " << c.t << ";\n";
		}
	}

	//if there are too many different colors involved, stop using variables
	if (rgbMap.size() > 100)
		rgbMap.clear();
	else
		ofile << rgbString.str();
	//same for radius
	if (radiusMap.size() > 100)
		radiusMap.clear();
	else
		ofile << radiusString.str();

	ofile<<"light_source {<0,0,1000>, color rgb 1 parallel point_at <0,0,0>}\n";

	for(int i=0; i<spheres.size(); i++)
	{
		M3DMatrix44f m0, m1, m2;
		AllSpheres.GetMatrix(m2);
		spheres[i].GetMatrix(m1);
		m3dMatrixMultiply44(m0, m2, m1);
		ofile << "sphere {<" << (-1)*(m0[12]) << "," << m0[13] << "," << m0[14] << ">,";
		auto iter = radiusMap.find(Spheres[i].radius);
		if (iter != radiusMap.end())
			ofile << "Radius" << iter->second;
		else
			ofile << Spheres[i].radius;

		ofile<<" pigment { color rgbt <";
		rgb c(Spheres[i]);
		auto iter2 = rgbMap.find(c);
		if (iter2 != rgbMap.end())
			ofile << "r" << iter2->second << ", g" << iter2->second << ", b" << iter2->second << ", t" << iter2->second;
		else
			ofile << Spheres[i].red << "," << Spheres[i].green << "," << Spheres[i].blue << "," << 1.0 - Spheres[i].transparency;
		ofile << "> } finish{reflection REFLECTION specular SPECULAR ambient AMBIENT diffuse DIFFUSE} } \n";
	}
	for(int i=0; i<Lines.size(); i++)
	{
		M3DMatrix44f m0, m1, m2;
		AllSpheres.GetMatrix(m2);
		for(int j=0; j<12; j++) m1[j]=0.0;
		m1[12]=Lines[i].x1;
		m1[13]=Lines[i].y1;
		m1[14]=Lines[i].z1;
		m1[15]=0.0;
		m3dMatrixMultiply44(m0, m2, m1);
		ofile<<"cylinder {<"<<(-1)*(m0[12])<<","<<m0[13]<<","<<m0[14]<<">,";
		double x1 = (-1)*(m0[12]), y1 = m0[13], z1 = m0[14];

		AllSpheres.GetMatrix(m2);
		for(int j=0; j<12; j++) m1[j]=0.0;
		m1[12]=Lines[i].x2;
		m1[13]=Lines[i].y2;
		m1[14]=Lines[i].z2;
		m1[15]=0.0;
		m3dMatrixMultiply44(m0, m2, m1);
		ofile<<"<"<<(-1)*(m0[12])<<","<<m0[13]<<","<<m0[14]<<">,"<<Lines[i].width<<" pigment { color rgb <";
		rgb c(Lines[i]);
		auto iter2 = rgbMap.find(c);
		if (iter2 != rgbMap.end())
			ofile << "r" << iter2->second << ", g" << iter2->second << ", b" << iter2->second << ", t" << iter2->second;
		else
			ofile << Lines[i].red << "," << Lines[i].green << "," << Lines[i].blue << "," << 1.0 - Lines[i].transparency;
		ofile<<"> } finish{reflection REFLECTION specular SPECULAR ambient AMBIENT diffuse DIFFUSE}  } \n";

		if (Lines[i].arrowHead)
		{
			double dx1 = (-1)*(m0[12]) - x1;
			double dy1 = (m0[13]) - y1;
			double dz1 = (m0[14]) - z1;
			double coeff = 8.0*Lines[i].width/std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
			ofile << "cone {<" << (-1)*(m0[12]) << "," << m0[13] << "," << m0[14] << ">," << 3.0*Lines[i].width << ",";
			ofile << "<" << (1.0 + coeff)*(-1)*(m0[12]) - coeff*x1 << "," << (1.0 + coeff)*m0[13] - coeff*y1 << "," << (1.0 + coeff)*m0[14] - coeff*z1 << ">," << 0.0 << " pigment { color rgb <";
			if (iter2 != rgbMap.end())
				ofile << "r" << iter2->second << ", g" << iter2->second << ", b" << iter2->second << ", t" << iter2->second;
			else
				ofile << Lines[i].red << "," << Lines[i].green << "," << Lines[i].blue << "," << 1.0 - Lines[i].transparency;
			ofile << "> } finish{reflection REFLECTION specular SPECULAR ambient AMBIENT diffuse DIFFUSE} } \n";
		}
	}
}
///////////////////////////////////////////////////////////////////////
// Draw spheres

//get depth of spheres[i]
double GetDepth(size_t i)
{
	M3DMatrix44f m0, m1, m2;
	frameCamera.GetCameraOrientation(m0);
	AllSpheres.GetMatrix(m1);
	m3dMatrixMultiply44(m2, m0, m1);
	spheres[i].GetMatrix(m1);
	m3dMatrixMultiply44(m0, m2, m1);
	return m0[14];
}
void DrawInhabitants()
{
	// Draw the spheres
	{
		//sort the spheres according to the depth of their centers
		auto DepthCompareFunc = [&] (const size_t i, const size_t j) -> bool
		{
			return GetDepth(i)<GetDepth(j);
		};
		std::sort(RenderOrder.begin(), RenderOrder.end(), DepthCompareFunc);
	}
	glPushMatrix();
	AllSpheres.ApplyActorTransform();
	glEnable(GL_LIGHTING);
	for(auto iter=RenderOrder.begin(); iter!=RenderOrder.end(); iter++)
	{
		size_t i=*iter;
		glPushMatrix();
		glColor4f(Spheres.at(i).red, Spheres.at(i).green, Spheres.at(i).blue, Spheres.at(i).transparency);
		spheres[i].ApplyActorTransform();
		glutSolidSphere(Spheres.at(i).radius, Precision, Precision);
		glPopMatrix();
	}
	glDisable(GL_LIGHTING);
	for(std::vector<line>::const_iterator iter=Lines.begin(); iter!=Lines.end(); iter++)
	{
		glColor4f(iter->red, iter->green, iter->blue, iter->transparency);
		glLineWidth(iter->width/ (Spheres[0].radius*0.1) );
		glBegin(GL_LINES);
		glVertex3d(iter->x1, iter->y1, iter->z1);
		glVertex3d(iter->x2, iter->y2, iter->z2);
		glEnd();
	}
	glPopMatrix();

}


// Called to draw scene
void RenderScene(void)
{
	// Clear the window with current clearing color
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	frameCamera.ApplyCameraTransform();
	// Position light before any other transformations
	glLightfv(GL_LIGHT0, GL_POSITION, fLightPos);
	DrawInhabitants();
	glPopMatrix();


	//2d stuff
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0, WindowWidth, 0, WindowHeight, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);

	(*pExtra2DStuffFunc)(WindowWidth, WindowHeight);

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	// Do the buffer Swap
	glutSwapBuffers();
	glutPostRedisplay();
}


static double OrthoWindowSize=1.0;
void ChangeToOrtho()
{
	OrthoProjection=true;
	// Reset the coordinate system before modifying
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glOrtho(fAspect*(-1)*OrthoWindowSize, fAspect*OrthoWindowSize, ((-1)*OrthoWindowSize), OrthoWindowSize, 0.1, 500);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();    
}
void ChangeToPerspective()
{
	OrthoProjection=false;
	// Reset the coordinate system before modifying
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(35.0f, fAspect, 0.1f, 500.0f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();    
}

// Respond to key strokes
void SpecialKeys(int key, int x, int y)
{
	if(key == GLUT_KEY_UP)
		if(OrthoProjection)
		{
			OrthoWindowSize*=0.9;
			ChangeToOrtho();
		}
		else
			frameCamera.MoveForward(0.1f);

	if(key == GLUT_KEY_DOWN)
		if(OrthoProjection)
		{
			OrthoWindowSize*=1.1;
			ChangeToOrtho();
		}
		else
			frameCamera.MoveForward(-0.1f);

	if(key == GLUT_KEY_LEFT)
		frameCamera.RotateLocalY(0.1);

	if(key == GLUT_KEY_RIGHT)
		frameCamera.RotateLocalY(-0.1);

	// Refresh the Window
	glutPostRedisplay();
}
void KeyStroke(unsigned char key, int x, int y)
{
	if(key == '+')
	{
		::Precision *= 2;
		if( ::Precision==0)
			::Precision=1;
		glutPostRedisplay();
	}
	if(key == '-')
	{
		::Precision /= 2;
		glutPostRedisplay();
	}
	if(key == 'q')
	{
		glutDestroyWindow( ::window_handle);
		glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);
		glutLeaveMainLoop();
	}
	if(key == 'o')
	{
		ChangeToOrtho();
		glutPostRedisplay();
	}
	if(key == 'p')
	{
		ChangeToPerspective();
		glutPostRedisplay();
	}
	if(key == 'a')
	{
		FrameRate/=1.5;
	}
	if(key == 's')
	{
		FrameRate*=1.5;
	}
	if(key == 'd')
	{
		Playing = !Playing;
		LastFrameTime = 0;
	}
	if(key == 'r')
	{
		Reverse = !Reverse;
	}
	if(key=='v')
		SavePovRay();
}

void IdleFunction(void)
{
	if(Playing)
	{
		long Time = GetPreciseClock();
		if(Time-LastFrameTime>FrameRate || LastFrameTime==0)
		{
			int NextFrame = DisplaySpheresMovie_CurrentFrame;
			size_t FrameChange;
			if (LastFrameTime == 0)
				FrameChange = 1;
			else
				FrameChange = (Time - LastFrameTime) / FrameRate;
			if (Reverse)
				NextFrame -= FrameChange;
			else
				NextFrame += FrameChange;
			if(NextFrame>=0 && NextFrame<MaxFrame)
			{
				PrepareFrame(NextFrame);
				glutPostRedisplay();
			}
			else
			{
				Playing=false;
			}
		}
	}
}
///////////////////////////////////////////////////////////
void ChangeSize(int w, int h)
{

	// Prevent a divide by zero, when window is too short
	// (you cant make a window of zero width).
	if(h == 0)
		h = 1;

	glViewport(0, 0, w, h);

	fAspect = (GLfloat)w / (GLfloat)h;
	WindowWidth=w;
	WindowHeight=h;

	// Set the clipping volume
	if(OrthoProjection)
		ChangeToOrtho();
	else
		ChangeToPerspective();
}

///////////////////////////////////////////////////////////
//this part of code deals with mouse dragging
static int CurrentMouseX, CurrentMouseY;
void PassiveMotion(int x, int y)
{
	::CurrentMouseX=x;
	::CurrentMouseY=y;
}
float Sensitivity = 0.003;
void Motion(int x, int y)
{
	int dx=x-CurrentMouseX;
	int dy=y-CurrentMouseY;
	::CurrentMouseX=x;
	::CurrentMouseY=y;

	AllSpheres.RotateWorld( ::Sensitivity*dy, 1.0, 0.0, 0.0);
	AllSpheres.RotateWorld( ::Sensitivity*dx, 0.0, 1.0, 0.0);
}

void DisplaySpheres_MultipleFrames(std::function<void(std::vector<sphere> & VSpheres, std::vector<line> & VLines, long FrameNumber)> GetFrameFunc, long NumFrames, std::function<void(size_t, size_t)> Extra2DStuff)
{
	if (SavePovrayAndExit)
	{
		::MaxFrame = NumFrames;
		pGetFrameFunc = &GetFrameFunc;
		pExtra2DStuffFunc = &Extra2DStuff;

		PlaceItems();

		for (long i = 0; i < NumFrames; i++)
		{
			std::stringstream ss;
			ss << SavePovray_Prefix << "frame_" << i;
			PrepareFrame(i);
			SavePovRay(ss.str());
		}
	}
	else
	{
		::MaxFrame = NumFrames;
		pGetFrameFunc = &GetFrameFunc;
		pExtra2DStuffFunc = &Extra2DStuff;


		int fake_argc = 1;
		char * fake_argv = "DisplaySpheres";
		glutInit(&fake_argc, &fake_argv);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
		glutInitWindowSize(WindowWidth, WindowHeight);
		window_handle = glutCreateWindow("");
		glutReshapeFunc(ChangeSize);
		glutDisplayFunc(RenderScene);
		glutSpecialFunc(SpecialKeys);
		glutIdleFunc(IdleFunction);
		glutMotionFunc(Motion);
		glutPassiveMotionFunc(PassiveMotion);
		glutKeyboardFunc(KeyStroke);

		SetupRC();

		glutMainLoop();
	}
}

#else

//These functions are disabled if OpenGL is disabled
void DisplayText(const std::string & text, double x, double y)
{
	return;
}
void DisplaySpheres_MultipleFrames(std::function<void(std::vector<sphere> & VSpheres, std::vector<line> & VLines, long FrameNumber)> GetFrameFunc, long NumFrames, std::function<void(size_t, size_t)> Extra2DStuff)
{
}
#endif

void DisplaySpheres_Movie(std::function<void(std::vector<sphere> & VSpheres, std::vector<line> & VLines, long FrameNumber)> GetFrameFunc, long NumFrames)
{
	auto ExtraDisplayFunc = [&] (size_t, size_t) ->void
	{
		std::stringstream ss;
		ss<<"Frame "<<DisplaySpheresMovie_CurrentFrame;
		DisplayText(ss.str());
	};
	DisplaySpheres_MultipleFrames(GetFrameFunc, NumFrames, ExtraDisplayFunc);
}
void DisplaySpheres(const std::vector<sphere> & VSpheres, const std::vector<line> & VLines, std::function<void(size_t, size_t)> Extra2DStuff)
{
	auto GetFrameFunc = [&](std::vector<sphere> & V1, std::vector<line> & V2, long FrameNumber) ->void
	{
		V1=VSpheres;
		V2=VLines;
	};
	DisplaySpheres_MultipleFrames(GetFrameFunc, 1, Extra2DStuff);
}
