#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
using namespace std;

#include <glad/glad.h> // be sure to include glad before any library that requires openGL
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Delaunay.h"
#include "Shader.h"

glm::mat4 world2screen;
glm::mat4 screen2world;
glm::vec2 cursorWorldPos;

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
	{
		glfwSetWindowShouldClose(window, true);
	}
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
	// convert (xpos, ypos) to (x,y)
	// get current screen dimension
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	glm::vec4 pos = screen2world * glm::vec4(xpos / width * 2 - 1.0, 1.0 - ypos / height * 2, 0.0, 1.0);
	cursorWorldPos.x = pos.x;
	cursorWorldPos.y = pos.y;
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}

double deltaTime;
double lastFrame;


int main()
{
	// read input from spreadsheet (.csv)
	fstream fin;
	fin.open("7.csv", ios::in);
	string line;
	int line_count = 0;

	vector<dt::Vertex> points;

	while (fin >> line)
	{
		if (line_count == 0)
		{
			cout << "Name of properties:" << line << endl;
		}
		else
		{
			// used for breaking words 
			stringstream s(line);
			string word;
			vector<string> row;
			while (getline(s, word, ',')) {

				// add all the column data 
				// of a row to a vector 
				row.push_back(word);
			}

			double x = stod(row[2]);
			double y = stod(row[3]);
			double height = stod(row[1]);
			cout << x << ' ' << y << ' ' << height << endl;
			points.emplace_back(x, y, height);
		}
		line_count++;
	}
	cout << "processed " << (line_count - 1) << " lines" << endl;

	const auto start = std::chrono::high_resolution_clock::now();
	dt::Delaunay delaunay(points);
	const auto end = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double> diff = end - start;
	unsigned int numOfTriangles = delaunay.getNumOfTriangles();
	std::cout << numOfTriangles << " triangles generated in " << diff.count() << "s\n";

	// get contour line
	fstream fout;
	fout.open("475_new.csv", ios::out);
	fout.precision(17);
	auto contours = delaunay.getContours(475);
	int num = 1;
	for (const auto& contour : contours)
	{
		fout << "contour No." << std::to_string(num++) << endl;
		for (const auto& vertex : contour)
		{
			fout << vertex.x << ',' << vertex.y << endl;
		}
	}
	fout.close();
	

	// initialize GLFW
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	//glfwWindowHint(GLFW_SAMPLES, 2);

	// GLFW window creation
	GLFWwindow* window = glfwCreateWindow(1280, 720, "Delaunay Triangulation", NULL, NULL);
	if (window == NULL)
	{
		cout << "Failed to create GLFW window" << endl;
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

	// register callback functions
	glfwSetKeyCallback(window, key_callback);
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	glfwSetCursorPosCallback(window, mouse_callback);

	// glad: load all OpenGL function pointers
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		cout << "Failed to initialize GLAD" << endl;
		return -1;
	}

	// compile shaders
	Shader meshShader("shaders/mesh.vs", "shaders/mesh.fs");

	// create world2screen and screen2world matrix
	world2screen = glm::mat3(1.0);
	world2screen = glm::scale(world2screen, glm::vec3(
		1.8 / (delaunay.getMaxX() - delaunay.getMinX()),
		1.8 / (delaunay.getMaxY() - delaunay.getMinY()),
		2.0 / (delaunay.getMaxH() - delaunay.getMinH())));
	world2screen = glm::translate(world2screen, glm::vec3(
		-(delaunay.getMaxX() + delaunay.getMinX()) / 2,
		-(delaunay.getMaxY() + delaunay.getMinY()) / 2,
		-(delaunay.getMaxH() + delaunay.getMinH()) / 2));
	screen2world = glm::inverse(world2screen);

	
	// get triangulation data
	const vector<dt::Vertex>& vertices = delaunay.getVertices();
	vector<unsigned int> indices = delaunay.getTriangleVertexIndices();
	
	// steup VAO, VBO and EBO
	unsigned int VAO, VBO, EBO;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);

	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(dt::Vertex), vertices.data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, sizeof(dt::Vertex), (void*)0);
	glBindVertexArray(0);

	/*unsigned int contourVBO, contourVAO;
	glGenVertexArrays(1, &contourVAO);
	glGenBuffers(1, &contourVBO);
	glBindVertexArray(contourVAO);
	glBindBuffer(GL_ARRAY_BUFFER, contourVBO);
	glBufferData(GL_ARRAY_BUFFER, numOfTriangles, NULL, GL_STREAM_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, sizeof(dt::Vertex), (void*)0);
	glBindVertexArray(0);*/

	// enable MSAA
	//glEnable(GL_MULTISAMPLE);

	// the render loop
	while (!glfwWindowShouldClose(window))
	{
		double currentFrame = glfwGetTime();
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		// process input here
		dt::Vertex q(cursorWorldPos.x, cursorWorldPos.y);
		//vector<vector<dt::Vertex>> contour_points;
		if (delaunay.interpolate(q))
		{
			cout << "x:" << q.x << " y:"<< q.y << " h:"<< q.h << endl;
			// get contour line
			/*auto contours = delaunay.getContours(q.h);
			for (const auto& contour : contours)
			{
				vector<dt::Vertex> points(contour.begin(), contour.end());
				contour_points.push_back(std::move(points));
			}*/
		}

		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		// draw stuff here
		// draw triangle meshes using GL_TRIANGLES
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		meshShader.Use();
		meshShader.SetMat4("world2screen", world2screen);
		glBindVertexArray(VAO);
		glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);

		// draw contours using GL_LINE_STRIP or GL_LINE_LOOP (depends on whether the countour is closed or not)
		// update contour VBO data
		/*if (!contour_points.empty())
		{
			meshShader.Use();
			meshShader.SetMat4("world2screen", world2screen);
			for (const auto& points : contour_points)
			{
				glBindVertexArray(contourVAO);
				glBufferSubData(GL_ARRAY_BUFFER, 0, points.size(), points.data());
				glDrawArrays(GL_LINE_STRIP, 0, points.size());
				glBindVertexArray(0);
			}
		}*/

		glfwSwapBuffers(window);
		glfwPollEvents();
	}
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
	//glDeleteVertexArrays(1, &contourVAO);
	//glDeleteBuffers(1, &contourVBO);
	glfwTerminate();
	return 0;
}