#include <GLFW/glfw3.h>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>

void geneCUBE(std::vector<double>& out_vertices, std::vector<unsigned int>& out_triangles) {
	out_vertices = {
		-1.0f,-1.0f,-1.0f,
		-1.0f,-1.0f, 2.0f,
		-1.0f, 1.0f, 2.0f,//1
		 1.0f, 1.0f,-1.0f,
		-1.0f,-1.0f,-1.0f,
		-1.0f, 1.0f,-1.0f,//2
		 1.0f,-1.0f, 2.0f,
		-1.0f,-1.0f,-1.0f,
		 1.0f,-1.0f,-1.0f,//3
		 1.0f, 1.0f,-1.0f,
		 1.0f,-1.0f,-1.0f,
		-1.0f,-1.0f,-1.0f,//4
		-1.0f,-1.0f,-1.0f,
		-1.0f, 1.0f, 2.0f,
		-1.0f, 1.0f,-1.0f,//5
		 1.0f,-1.0f, 2.0f,
		-1.0f,-1.0f, 2.0f,
		-1.0f,-1.0f,-1.0f,//6
		-1.0f, 1.0f, 2.0f,
		-1.0f,-1.0f, 2.0f,
		 1.0f,-1.0f, 2.0f,//7
		 1.0f, 1.0f, 2.0f,
		 1.0f,-1.0f,-1.0f,
		 1.0f, 1.0f,-1.0f,//8
		 1.0f,-1.0f,-1.0f,
		 1.0f, 1.0f, 2.0f,
		 1.0f,-1.0f, 2.0f,//9
		 1.0f, 1.0f, 2.0f,
		 1.0f, 1.0f,-1.0f,
		-1.0f, 1.0f,-1.0f,//10
		 1.0f, 1.0f, 2.0f,
		-1.0f, 1.0f,-1.0f,
		-1.0f, 1.0f, 2.0f,//11
		 1.0f, 1.0f, 2.0f,
		-1.0f, 1.0f, 2.0f,
		 1.0f,-1.0f, 2.0f//12
	};

	for (int i = 0; i < 36; i++)
	{
		out_triangles.push_back(i);
	}
	std::cout << out_triangles.size() << std::endl;
}

bool loadOBJ(
	const char* path,

	std::vector<double>& out_vertices,
	std::vector<unsigned int>& out_triangles
) {
	printf("Loading OBJ file %s...\n", path);


	//std::vector<double> temp_vertices;
	std::vector<unsigned int> vertexIndices;

	FILE* file = fopen(path, "r");
	if (file == NULL) {
		printf("cannot open the file !\n");
		getchar();
		return false;
	}

	while (1) {

		char lineHeader[128];
		// read the first word 
		int res = fscanf(file, "%s", lineHeader);
		if (res == EOF)
			break; // EOF = End Of File Quit 
		// else : parse lineHeader

		if (strcmp(lineHeader, "v") == 0) {
			double x, y, z;
			fscanf(file, "%lf %lf %lf\n", &x, &y, &z);
			std::cout << x << " " << y << " " << z << " " << std::endl;
			out_vertices.push_back(x);
			out_vertices.push_back(y);
			out_vertices.push_back(z);
		}
		else if (strcmp(lineHeader, "f") == 0) {
			unsigned int vertexIndex[3];
			int matches = fscanf(file, "%d %d %d/\n", &vertexIndex[0], &vertexIndex[1], &vertexIndex[2]);
			if (matches != 3) {
				printf("File can't be read");
				fclose(file);
				return false;
			}
			//std::cout << vertexIndex[0] << " " << vertexIndex[1] << " " << vertexIndex[2] << " " << std::endl;
			vertexIndices.push_back(vertexIndex[0] - 1 );
			vertexIndices.push_back(vertexIndex[1] - 1);
			vertexIndices.push_back(vertexIndex[2] - 1);
		}
		else {
			// comment
			char stupidBuffer[1000];
			fgets(stupidBuffer, 1000, file);
		}

	}
	//std::cout << temp_vertices[0] << " " << temp_vertices[1] << " " << temp_vertices[2] << std::endl;
	std::cout << "\n print vertices  index is \n";
	std::cout << vertexIndices.size() << std::endl;

	out_triangles = vertexIndices;
	fclose(file);
	return true;
}

static void error_callback(int error, const char* description)
{
  fputs(description, stderr);
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    glfwSetWindowShouldClose(window, GL_TRUE);
}


int main(void)
{
	std::vector<double> points;
	std::vector<unsigned int> indexs;
	//loadOBJ("bunny.obj", points, indexs);
	geneCUBE(points, indexs);

  GLFWwindow* window;
  glfwSetErrorCallback(error_callback);
  if (!glfwInit())
    exit(EXIT_FAILURE);
  window = glfwCreateWindow(840, 680, "Assignment_3_transform_cuboid", NULL, NULL);
  if (!window)
  {
    glfwTerminate();
    exit(EXIT_FAILURE);
  }

  glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, key_callback);
  int iframe = 0;
  while (!glfwWindowShouldClose(window))
  {
    float ratio;
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    ratio = width / (float) height;

    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-ratio*5.f, ratio * 5.f, -1.f * 5.f, 1.f * 5.f, 1.f * 5.f, -1.f * 5.f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    //glRotatef((float) glfwGetTime() * 50.f, 0.f, 0.f, 1.f);
    glBegin(GL_TRIANGLES);
	glColor3f(1.0f,1.0f,1.0f);
	for (unsigned int i = 0; i < indexs.size() / 3; ++i) {

		unsigned int i0 = indexs[i * 3 + 0];
		unsigned int i1 = indexs[i * 3 + 1];
		unsigned int i2 = indexs[i * 3 + 2];
		glVertex3f(points[i0 * 3 + 0], points[i0 * 3 + 1], points[i0 * 3 + 2]);
		glVertex3f(points[i1 * 3 + 0], points[i1 * 3 + 1], points[i1 * 3 + 2]);
		glVertex3f(points[i2 * 3 + 0], points[i2 * 3 + 1], points[i2 * 3 + 2]);
	}
		glEnd();
    glfwSwapBuffers(window);
    glfwPollEvents();
  }
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}