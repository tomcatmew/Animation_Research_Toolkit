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

std::vector<double> matrix_mutiple(std::vector<std::vector<double>> affine_m, std::vector<double> coord) {
	std::vector<double> result;
	for (int i = 0; i < affine_m.size(); i += 1) {
		double tempt_sum = 0;
		for (int j = 0; j < affine_m[0].size(); j += 1) {
			tempt_sum = tempt_sum + coord[j] * affine_m[i][j];
		}
		result.push_back(tempt_sum);
	}
	return result;
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
  int rotate_degree = -1;
  double move_step = 0.1f;
  bool reach_end = false;

  glEnable(GL_DEPTH_TEST);
  // Accept fragment if it closer to the camera than the former one
  glDepthFunc(GL_LESS);

  double scale_ratio = 5.f;
  while (!glfwWindowShouldClose(window))
  {

	  //get current rotation degree
	  if (rotate_degree < 359)
	  {
		  rotate_degree += 1;
	  }
	  else {
		  rotate_degree = 0;
	  }

	  //get current moving position 
	  if ((move_step < 4.f) & (reach_end == false)){
		  move_step += 0.05f;
	  }
	  else {
		  reach_end = true;
	  }

	  if ((move_step > -5.f) & (reach_end == true)) {
		  move_step -= 0.05f;
	  }
	  else {
		  reach_end = false;
	  }

	  // set up affine matrix and calculate new coordinate 

	  //affine matrix Y rotation
	  std::vector<std::vector<double>> affineY;
	  std::vector<double> row1y = { cos(rotate_degree * 3.14159 / 180),0,sin(rotate_degree * 3.14159 / 180),0 };
	  std::vector<double> row2y = { 0, 1, 0, 0 };
	  std::vector<double> row3y = { -sin(rotate_degree * 3.14159 / 180), 0,cos(rotate_degree * 3.14159 / 180), 0 };
	  std::vector<double> row4y = { 0, 0, 0, 1 };
	  affineY.push_back(row1y);
	  affineY.push_back(row2y);
	  affineY.push_back(row3y);
	  affineY.push_back(row4y);

	  //affine matrix X rotation and let cuboid moves on X axis
	  std::vector<std::vector<double>> affineX;
	  std::vector<double> row1x = { 1,0,0,move_step };
	  std::vector<double> row2x = { 0, cos(rotate_degree * 3.14159 / 180),-sin(rotate_degree * 3.14159 / 180), 0 };
	  std::vector<double> row3x = { 0, sin(rotate_degree * 3.14159 / 180), cos(rotate_degree * 3.14159 / 180), 0 };
	  std::vector<double> row4x = { 0, 0, 0, 1 };
	  affineX.push_back(row1x);
	  affineX.push_back(row2x);
	  affineX.push_back(row3x);
	  affineX.push_back(row4x);

	  //store the vertices list of rotated and transformed vertices
	  std::vector<double> transformed_points;

	  for (int i = 0; i < points.size() / 3; i += 1) {
		  std::vector<double> tempt_coord;
		  std::vector<double> new_coordY;
		  std::vector<double> new_coordX;

		  tempt_coord.push_back(points[i * 3 + 0]);
		  tempt_coord.push_back(points[i * 3 + 1]);
		  tempt_coord.push_back(points[i * 3 + 2]);
		  tempt_coord.push_back(1.f);

		  //first calculating the new coordinate from Y rotation
		  new_coordY = matrix_mutiple(affineY,tempt_coord);
		  //then calculating the new coordinate from X rotation
		  new_coordX = matrix_mutiple(affineX, new_coordY);
		  //std::cout << new_coord.size() << std::endl;

		  //put the rotated coordinate into the vertex list
		  for (int j = 0; j < new_coordX.size() - 1 ; j += 1) {
			  transformed_points.push_back(new_coordX[j]);
		  }
	  }
	  // end calculation of affine matrix 


    float ratio;
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    ratio = width / (float) height;

    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-ratio* scale_ratio, ratio * scale_ratio, -1.f * scale_ratio, 1.f * scale_ratio, 1.f * scale_ratio, -1.f * scale_ratio);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    //glRotatef((float) glfwGetTime() * 50.f, 0.f, 0.f, 1.f);
    glBegin(GL_TRIANGLES);

	for (unsigned int i = 0; i < indexs.size() / 3; ++i) {

		unsigned int i0 = indexs[i * 3 + 0];
		unsigned int i1 = indexs[i * 3 + 1];
		unsigned int i2 = indexs[i * 3 + 2];
		glColor3f(1.0f, 1.0f, 1.0f);
		glVertex3f(transformed_points[i0 * 3 + 0], transformed_points[i0 * 3 + 1], transformed_points[i0 * 3 + 2]);
		glColor3f(0.7f, 1.f, 0.f);
		glVertex3f(transformed_points[i1 * 3 + 0], transformed_points[i1 * 3 + 1], transformed_points[i1 * 3 + 2]);
		glColor3f(0.1f, 0.5f, 1.f);
		glVertex3f(transformed_points[i2 * 3 + 0], transformed_points[i2 * 3 + 1], transformed_points[i2 * 3 + 2]);
	}
		glEnd();
    glfwSwapBuffers(window);
    glfwPollEvents();
  }
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}