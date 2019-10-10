
#include <iostream>
#include <vector>
using namespace std;


// Helper graphic libraries
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/transform.hpp>
#include "graphics.h"
#include "shapes.h"
#include "aStar.cpp"
#include <ctime>
# define M_PI           3.14159265358979323846 

// MAIN FUNCTIONS
void startup();
void updateCamera();
void updateSceneElements();
void renderScene();
void ball_bounce();
void explode();
void crowdSim();
void displayAStar();
void updateActor();
// CALLBACK FUNCTIONS
void onResizeCallback(GLFWwindow* window, int w, int h);
void onKeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
void onMouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
void onMouseMoveCallback(GLFWwindow* window, double x, double y);
void onMouseWheelCallback(GLFWwindow* window, double xoffset, double yoffset);


// VARIABLES

bool        quit = false;
float        deltaTime = 0.0f;    // Keep track of time per frame.
float        lastTime = 0.0f;    // variable to keep overall time.
float		timer = 0.0f;
bool        keyStatus[1024];    // Hold key status.
GLfloat startTime = (GLfloat)glfwGetTime();

//bounce variables
glm::vec3 pos = glm::vec3(0.0f, 1.0f, 0.0f);
glm::vec3 vel = glm::vec3(0.0f, 0.0f, 0.0f);
glm::vec3 acc = glm::vec3(0.0f, -9.81f, 0.0f); 
float yk = -0.0005;
bool bounce = false;
float distance1;

//particle variables
bool explosion = false;
const int numberOfParticles = 180;
glm::vec3 epos[numberOfParticles];
glm::vec3 evel[numberOfParticles];
glm::vec3 eacc[numberOfParticles];

//crowd variables
bool crowd = false;
const int numberOfObjects = 4;
glm::vec3 apos[numberOfObjects];
glm::vec3 avel[numberOfObjects];
glm::vec3 aacc[numberOfObjects];


//A* variables
map m;
point s, e(MAP_SIZE-1, MAP_SIZE-1);
aStar as;
std::list<point> path;
int xPosition = 0;
int zPosition = 0;
int path_number = 0;
bool stop = true;
bool a_star = false;


// MAIN GRAPHICS OBJECT
Graphics    myGraphics;        // Runing all the graphics in this object


// OBJECTS
Cube        particle[numberOfParticles];
Cube		tile[MAP_SIZE][MAP_SIZE];
Sphere		piece;
Sphere        ball;
Arrow        arrowX;
Arrow        arrowY;
Arrow        arrowZ;
Cube        myFloor;
Line        myLine;
Cylinder    myCylinder;




int main()
{

	int errorGraphics = myGraphics.Init();        // Launch window and graphics context
	if (errorGraphics) return 0;                // Close if something went wrong...

	startup();                                    // Setup all necessary information for startup (aka. load texture, shaders, models, etc).



	// MAIN LOOP run until the window is closed
	while (!quit) {

		// Update the camera tranform based on interactive inputs or by following a predifined path.
		updateCamera();

		// Update position, orientations and any other relevant visual state of any dynamic elements in the scene.
		updateSceneElements();

		// Render a still frame into an off-screen frame buffer known as the backbuffer.
		renderScene();

		// Swap the back buffer with the front buffer, making the most recently rendered image visible on-screen.
		glfwSwapBuffers(myGraphics.window);        // swap buffers (avoid flickering and tearing)

	}


	myGraphics.endProgram();            // Close and clean everything up...

   // cout << "\nPress any key to continue...\n";
   // cin.ignore(); cin.get(); // delay closing console to read debugging errors.

	return 0;
}


void startup() {
	// Keep track of the running time
	GLfloat currentTime = (GLfloat)glfwGetTime();    // retrieve timelapse
	deltaTime = currentTime;                        // start delta time
	lastTime = currentTime;                            // Save for next frame calculations.

	// Callback graphics and key update functions - declared in main to avoid scoping complexity.
	// More information here : https://www.glfw.org/docs/latest/input_guide.html
	glfwSetWindowSizeCallback(myGraphics.window, onResizeCallback);            // Set callback for resize
	glfwSetKeyCallback(myGraphics.window, onKeyCallback);                    // Set Callback for keys
	glfwSetMouseButtonCallback(myGraphics.window, onMouseButtonCallback);    // Set callback for mouse click
	glfwSetCursorPosCallback(myGraphics.window, onMouseMoveCallback);        // Set callback for mouse move
	glfwSetScrollCallback(myGraphics.window, onMouseWheelCallback);            // Set callback for mouse wheel.

	// Calculate proj_matrix for the first time.
	myGraphics.aspect = (float)myGraphics.windowWidth / (float)myGraphics.windowHeight;
	myGraphics.proj_matrix = glm::perspective(glm::radians(50.0f), myGraphics.aspect, 0.1f, 1000.0f);

	// Load Geometry.
	for (int x = 0; x < numberOfParticles; x++) {
		particle[x].Load();
	}
	
	
	ball.Load();
	ball.fillColor = glm::vec4(0.0f, 1.0f, 0.0f, 1.0f);    // You can change the shape fill colour, line colour or linewidth

	piece.Load();
	piece.fillColor = glm::vec4(1.0f, 1.0f, 0.0f, 1.0f);    // You can change the shape fill colour, line colour or linewidth

	
	myFloor.Load();
	myFloor.fillColor = glm::vec4(130.0f / 255.0f, 96.0f / 255.0f, 61.0f / 255.0f, 1.0f);    // Sand Colour
	myFloor.lineColor = glm::vec4(130.0f / 255.0f, 96.0f / 255.0f, 61.0f / 255.0f, 1.0f);    // Sand again

	myCylinder.Load();
	myCylinder.fillColor = glm::vec4(0.7f, 0.7f, 0.7f, 1.0f);
	myCylinder.lineColor = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);


	// Optimised Graphics
	myGraphics.SetOptimisations();        // Cull and depth testing

}

void updateCamera() {

	// calculate movement for FPS camera
	GLfloat xoffset = myGraphics.mouseX - myGraphics.cameraLastX;
	GLfloat yoffset = myGraphics.cameraLastY - myGraphics.mouseY;    // Reversed mouse movement
	myGraphics.cameraLastX = (GLfloat)myGraphics.mouseX;
	myGraphics.cameraLastY = (GLfloat)myGraphics.mouseY;

	GLfloat sensitivity = 0.05f;
	xoffset *= sensitivity;
	yoffset *= sensitivity;

	myGraphics.cameraYaw += xoffset;
	myGraphics.cameraPitch += yoffset;

	// check for pitch out of bounds otherwise screen gets flipped
	if (myGraphics.cameraPitch > 89.0f) myGraphics.cameraPitch = 89.0f;
	if (myGraphics.cameraPitch < -89.0f) myGraphics.cameraPitch = -89.0f;

	// Calculating FPS camera movement (See 'Additional Reading: Yaw and Pitch to Vector Calculations' in VISION)
	glm::vec3 front;
	front.x = cos(glm::radians(myGraphics.cameraYaw)) * cos(glm::radians(myGraphics.cameraPitch));
	front.y = sin(glm::radians(myGraphics.cameraPitch));
	front.z = sin(glm::radians(myGraphics.cameraYaw)) * cos(glm::radians(myGraphics.cameraPitch));

	myGraphics.cameraFront = glm::normalize(front);

	// Update movement using the keys
	GLfloat cameraSpeed = 1.5f * deltaTime;
	if (keyStatus[GLFW_KEY_W]) myGraphics.cameraPosition += cameraSpeed * myGraphics.cameraFront;
	if (keyStatus[GLFW_KEY_S]) myGraphics.cameraPosition -= cameraSpeed * myGraphics.cameraFront;
	if (keyStatus[GLFW_KEY_A]) myGraphics.cameraPosition -= glm::normalize(glm::cross(myGraphics.cameraFront, myGraphics.cameraUp)) * cameraSpeed;
	if (keyStatus[GLFW_KEY_D]) myGraphics.cameraPosition += glm::normalize(glm::cross(myGraphics.cameraFront, myGraphics.cameraUp)) * cameraSpeed;

	// IMPORTANT PART
	// Calculate my view matrix using the lookAt helper function
	myGraphics.viewMatrix = glm::lookAt(myGraphics.cameraPosition,        // eye
		myGraphics.cameraPosition + myGraphics.cameraFront,                // centre
		myGraphics.cameraUp);                                            // up

}

void updateSceneElements() {

	glfwPollEvents();                                // poll callbacks

	// Calculate frame time/period -- used for all (physics, animation, logic, etc).
	GLfloat currentTime = (GLfloat)glfwGetTime();    // retrieve timelapse
	deltaTime = currentTime - lastTime;                // Calculate delta time
	lastTime = currentTime;                            // Save for next frame calculations.
	GLfloat totalTime = currentTime - startTime;
	// Do not forget your ( T * R * S ) http://www.opengl-tutorial.org/beginners-tutorials/tutorial-3-matrices/
	timer += deltaTime;

	if (bounce) {
		ball_bounce();
	}
	if (!bounce) {
		pos.y = 0.5f;
	}

	if (explosion) {
		explode();
	}
	if (crowd) {
		crowdSim();
	}
	//update actor position
	if (a_star && !stop) {
		updateActor();
	}


	// Calculate particle position
	for (int x = 0; x < numberOfParticles; x++) {
		glm::mat4 mv_matrix_cube =
			glm::translate(glm::vec3(epos[x])) *
			glm::scale(glm::vec3(0.02f, 0.02f, 0.02f)) *
			glm::mat4(1.0f);
		particle[x].mv_matrix = myGraphics.viewMatrix * mv_matrix_cube;
		particle[x].proj_matrix = myGraphics.proj_matrix;
	}

	//calculate tiles position
	for (int x = 0; x < MAP_SIZE; x++) {
		for (int z = 0; z< MAP_SIZE; z++) {
			glm::mat4 mv_matrix_cube;
			if (m(x, z) == 0) {
				mv_matrix_cube =
					glm::translate(glm::vec3(x, 0.5f, z)) *
					glm::mat4(1.0f);
			}
			else {
				mv_matrix_cube =
					glm::translate(glm::vec3(x, 0.5f, z)) *
					glm::scale(glm::vec3(1.0f, 2.0f, 1.0f)) *
					glm::mat4(1.0f);
			}
			tile[x][z].mv_matrix = myGraphics.viewMatrix * mv_matrix_cube;
			tile[x][z].proj_matrix = myGraphics.proj_matrix;
		}
	}


	// calculate ball movement
	glm::mat4 mv_matrix_sphere;
	if (distance1 > 0) {
		mv_matrix_sphere =
			glm::translate(glm::vec3(pos)) *
			glm::rotate(glm::radians(-distance1 * 100), glm::vec3(0.0f, 0.0f, 1.0f)) *
			glm::mat4(1.0f);
	}
	else {
		mv_matrix_sphere =
			glm::translate(glm::vec3(pos)) *

			glm::mat4(1.0f);
	}
	ball.mv_matrix = myGraphics.viewMatrix * mv_matrix_sphere;
	ball.proj_matrix = myGraphics.proj_matrix;


	//calculte actor position
	mv_matrix_sphere =
		glm::translate(glm::vec3(xPosition, 1.5f, zPosition)) *
		glm::mat4(1.0f);
	piece.mv_matrix = myGraphics.viewMatrix * mv_matrix_sphere;
	piece.proj_matrix = myGraphics.proj_matrix;


	// Calculate floor position and resize
	myFloor.mv_matrix = myGraphics.viewMatrix *
		glm::translate(glm::vec3(0.0f, 0.0f, 0.0f)) *
		glm::scale(glm::vec3(1000.0f, 0.001f, 1000.0f)) *
		glm::mat4(1.0f);
	myFloor.proj_matrix = myGraphics.proj_matrix;

	// Calculate cylinder
	myCylinder.mv_matrix = myGraphics.viewMatrix *
		glm::translate(glm::vec3(-1.0f, 0.5f, 2.0f)) *
		glm::mat4(1.0f);
	myCylinder.proj_matrix = myGraphics.proj_matrix;
	


	if (glfwWindowShouldClose(myGraphics.window) == GL_TRUE) quit = true; // If quit by pressing x on window.

}
//-------------------------------------------------------------------------------------------------
//---------------FUNCTIONS-----------------------------------------
//update ball position
void ball_bounce() {

	if (pos.y < 0.5 && vel.y < 0) {
		vel.y = vel.y*-0.7;
		if (sqrt(vel.y*vel.y) < 0.1) {
			acc.y = 0;
			vel.y = 0;
		}

	}
	float oldPos = pos.x;
	vel += acc * deltaTime;
	if (vel.x < 0) {
		vel.x = 0;
	}
	else {
		acc.x += yk * vel.x*vel.x;
	}

	pos += vel * deltaTime;

	distance1 = pos.x;
	distance1 = sqrt(distance1*distance1);
}

//update the particles positions
void explode() {

	for (int x = 0; x < numberOfParticles; x++) {
		
		if (epos[x].y < 0.02 && evel[x].y < 0) {
			evel[x].y = evel[x].y*-0.1;
			if (sqrt(evel[x].y*evel[x].y) < 0.1) {
				eacc[x].y = 0;
				evel[x].y = 0;
				epos[x].y = 0.01;
			}


		}

		evel[x] += eacc[x] * deltaTime;
		if (evel[x].x < 0.05 && evel[x].x > -0.05) {
			evel[x].x = 0;
			eacc[x].x = 0;
		}
		else {
			if (evel[x].x < 0) {
				eacc[x].x *= -1;
			}

			if (evel[x].x < 0) {
				eacc[x].x *= -1;
			}
		}
		if (evel[x].z < 0.05 && evel[x].z > -0.05) {
			evel[x].z = 0;
			eacc[x].z = 0;
		}
		else {
			if (evel[x].z < 0) {
				eacc[x].z *= -1;
			}

			if (evel[x].z < 0) {
				eacc[x].z *= -1;
			}
		}
		epos[x] += evel[x] * deltaTime;

	}

}
void crowdSim() {

}
//update actor position
void updateActor() {
	if (timer > 0.5) {
		timer = 0;
		std::list<point>::iterator it = path.begin();


		std::advance(it, path_number);

		if (it != path.end()) {
			path_number++;
			xPosition = (*it).x;
			zPosition = (*it).y;
		}
		else {
			stop = true;
		}
	}
}
//display tiles
void displayAStar() {
	map n;
	m = n;
	aStar newStar;
	as = newStar;
	path.clear(); 
	if (as.search(s, e, m)) {
		as.path(path);
		stop = false;
		a_star = true;
		for (int x = 0; x < MAP_SIZE; x++) {
			for (int z = 0; z < MAP_SIZE; z++) {
				tile[x][z].Load();
				if (m(x, z) == 0) {
					tile[x][z].fillColor = glm::vec4(0.0f, 0.0f, 1.0f, 1.0f); tile[x][z].lineColor = glm::vec4(0.0f, 0.0f, 1.0f, 1.0f);
				}
				else {
					tile[x][z].fillColor = glm::vec4(1.0f, 0.0f, 0.0f, 1.0f); tile[x][z].lineColor = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
				}

			}
		}
		for (std::list<point>::iterator i = path.begin(); i != path.end(); i++) {
			tile[(*i).x][(*i).y].fillColor = glm::vec4(0.0f, 1.0f, 0.0f, 1.0f); tile[(*i).x][(*i).y].lineColor = glm::vec4(0.0f, 1.0f, 0.0f, 1.0f);
		}
	}
	else {
		displayAStar();
	}
}
void renderScene() {
	// Clear viewport - start a new frame.
	myGraphics.ClearViewport();

	// Draw objects in screen
	myFloor.Draw();
	for (int x = 0; x < numberOfParticles; x++) {
		particle[x].Draw();
	}
	if (a_star) {
		for (int x = 0; x < MAP_SIZE; x++) {
			for (int z = 0; z < MAP_SIZE; z++) {
				tile[x][z].Draw();
			}
		}
		piece.Draw();
	}
	if (bounce) {
		ball.Draw();
	}
	//arrowX.Draw();
	//arrowY.Draw();
	//arrowZ.Draw();

	//myLine.Draw();
	myCylinder.Draw();
}


// CallBack functions low level functionality.
void onResizeCallback(GLFWwindow* window, int w, int h) {    // call everytime the window is resized
	//myGraphics.windowWidth = w;
	//myGraphics.windowHeight = h;

	glfwGetFramebufferSize(window, &myGraphics.windowWidth, &myGraphics.windowHeight);

	myGraphics.aspect = (float)w / (float)h;
	myGraphics.proj_matrix = glm::perspective(glm::radians(50.0f), myGraphics.aspect, 0.1f, 1000.0f);
}

void onKeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) { // called everytime a key is pressed
	if (action == GLFW_PRESS) keyStatus[key] = true;
	else if (action == GLFW_RELEASE) keyStatus[key] = false;

	// toggle showing mouse.
	if (keyStatus[GLFW_KEY_M]) myGraphics.ToggleMouse();

	//if b key is pressed start bouncing ball.
	if (keyStatus[GLFW_KEY_B]) {
		pos = glm::vec3(0.0f, 4.0f, 0.0f);
		vel = glm::vec3(1.0f, 0.0f, 0.0f);
		acc = glm::vec3(0.0f, -9.81f, 0.0f);
		if (bounce) {
			bounce = false;
		}
		else {
			bounce = true;
		}
	}
	//if X key is pressed start the particle explosion.
	if (keyStatus[GLFW_KEY_X]) {
		float angle = 360 / numberOfParticles;
		for (int x = 0; x < numberOfParticles; x++) {
			float r = 3 * (rand() / (RAND_MAX + 1.0));
			float h = 5 * (rand() / (RAND_MAX + 1.0));
			epos[x] = glm::vec3(0.0f, 0.01f, 0.0f);
			evel[x] = glm::vec3((cos(angle*x)*r), h, (sin(angle*x) * r));
			eacc[x] = glm::vec3(-(cos(angle*x)*r), -9.81f, -(sin(angle*x) * r));
		}
		if (!explosion) {


			explosion = true;
		}
		else {
			explosion = false;
		}
	}
	//if V key is pressed start the a* search simulation
	if (keyStatus[GLFW_KEY_V]) {
		if (a_star) {
			a_star = false;
			stop = true;
			path_number = 0;
			xPosition = 0;
			zPosition = 0;
		}
		else {
			displayAStar();
		}

	}
	if (keyStatus[GLFW_KEY_C]) {
		if (crowd) {

		}
		else {

		}
	}
	// If exit key pressed.
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GLFW_TRUE);
}

void onMouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {

}

void onMouseMoveCallback(GLFWwindow* window, double x, double y) {
	int mouseX = static_cast<int>(x);
	int mouseY = static_cast<int>(y);

	myGraphics.mouseX = mouseX;
	myGraphics.mouseY = mouseY;

	// helper variables for FPS camera
	if (myGraphics.cameraFirstMouse) {
		myGraphics.cameraLastX = (GLfloat)myGraphics.mouseX; myGraphics.cameraLastY = (GLfloat)myGraphics.mouseY; myGraphics.cameraFirstMouse = false;
	}
}

void onMouseWheelCallback(GLFWwindow* window, double xoffset, double yoffset) {
	int yoffsetInt = static_cast<int>(yoffset);
}


////AStar.cpp
#include <list>
#include <algorithm>
#include <iostream>
#include <random>
#define MAP_SIZE  20
#define BLOCK_NUMBERS 100

class point {
public:
	point(int a = 0, int b = 0) { x = a; y = b; }
	bool operator ==(const point& o) { return o.x == x && o.y == y; }
	point operator +(const point& o) { return point(o.x + x, o.y + y); }
	int x, y;
};
//create map for the a* search
class map {
public:
	 map() {
		bool random = true;
		char t[MAP_SIZE][MAP_SIZE];
		for (int x = 0; x < MAP_SIZE; x++) {
			for (int z = 0; z < MAP_SIZE; z++) {
				t[x][z] = 0;
			}
		}
		
		if (random) {
			std::random_device rd;
			std::mt19937 eng(rd()); 
			std::uniform_int_distribution<> distr(0, MAP_SIZE-1);

			for (int x = 0; x < BLOCK_NUMBERS; x++) {
				int r = distr(eng);
				int c = distr(eng);
				if (t[r][c] == 0) {
					if (r == 0 && c == 0) {
						
					}
					if (r == MAP_SIZE-1 && c == MAP_SIZE-1) {
						
					}
					else {
						t[r][c] = 1;
					}
				}
				else {
					x--;
				}
			}
		}
		w = h = MAP_SIZE;
		for (int r = 0; r < h; r++)
			for (int s = 0; s < w; s++)
				m[s][r] = t[r][s];
		
	}
	int operator() (int x, int y) { return m[x][y]; }
	char m[MAP_SIZE][MAP_SIZE];
	int w, h;
};

class node {
public:
	bool operator == (const node& o) { return pos == o.pos; }
	bool operator == (const point& o) { return pos == o; }
	bool operator < (const node& o) { return dist + cost < o.dist + o.cost; }
	point pos, parent;
	int dist, cost;
};

//a* search algorithm
class aStar {
public:
	aStar() {
		neighbours[1] = point(0, -1); neighbours[2] = point(-1, 0);
		neighbours[3] = point(0, 1); neighbours[4] = point(1, 0);
	}
	//calculate manhattan distance
	int calculateDistance(point& p) {
		
		int x = end.x - p.x, y = end.y - p.y;
		return(x * x + y * y);
	}

	bool isValid(point& p) {
		return (p.x > -1 && p.y > -1 && p.x < m.w && p.y < m.h);
	}

	bool existPoint(point& p, int cost) {
		std::list<node>::iterator i;
		i = std::find(closed.begin(), closed.end(), p);
		if (i != closed.end()) {
			if ((*i).cost + (*i).dist < cost) return true;
			else { closed.erase(i); return false; }
		}
		i = std::find(open.begin(), open.end(), p);
		if (i != open.end()) {
			if ((*i).cost + (*i).dist < cost) return true;
			else { open.erase(i); return false; }
		}
		return false;
	}

	bool fillOpen(node& n) {
		int stepCost, nc, dist;
		point neighbour;

		for (int x = 0; x < 8; x++) {
			
			stepCost = x < 4 ? 1 : 1;
			neighbour = n.pos + neighbours[x];
			if (neighbour == end) return true;

			if (isValid(neighbour) && m(neighbour.x, neighbour.y) != 1) {
				nc = stepCost + n.cost;
				dist = calculateDistance(neighbour);
				if (!existPoint(neighbour, nc + dist)) {
					node m;
					m.cost = nc; m.dist = dist;
					m.pos = neighbour;
					m.parent = n.pos;
					open.push_back(m);
				}
			}
		}
		return false;
	}

	bool search(point& s, point& e, map& mp) {
		node n; end = e; start = s; m = mp;
		n.cost = 0; n.pos = s; n.parent = 0; n.dist = calculateDistance(s);
		open.push_back(n);
		while (!open.empty()) {
			node n = open.front();
			open.pop_front();
			closed.push_back(n);
			if (fillOpen(n)) return true;
		}
		return false;
	}
	//calculate the path from beginning to end
	int path(std::list<point>& path) {
		path.push_front(end);
		int cost = 1 + closed.back().cost;
		path.push_front(closed.back().pos);
		point parent = closed.back().parent;

		for (std::list<node>::reverse_iterator i = closed.rbegin(); i != closed.rend(); i++) {
			if ((*i).pos == parent && !((*i).pos == start)) {
				path.push_front((*i).pos);
				parent = (*i).parent;
			}
		}
		path.push_front(start);
		return cost;
	}

	map m; point end, start;
	point neighbours[8];
	std::list<node> open;
	std::list<node> closed;
};
