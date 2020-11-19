#include <GLFW/glfw3.h>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>

void geneCUBE_length(double x, double y, double z, std::vector<double>& out_vert, std::vector<unsigned int>& out_tri){
	out_vert = {
	-x / 2.0f,-y / 2.0f,-z / 2.0f,  //0
	-x / 2.0f,-y / 2.0f, z / 2.0f,  //1
	-x / 2.0f, y / 2.0f, z / 2.0f,  //2
	 x / 2.0f, y / 2.0f,-z / 2.0f,  //3 
	-x / 2.0f, y / 2.0f,-z / 2.0f,  //4
	 x / 2.0f,-y / 2.0f, z / 2.0f,  //5
	 x / 2.0f,-y / 2.0f,-z / 2.0f,  //6
	 x / 2.0f, y / 2.0f, z / 2.0f   //7
	};

	out_tri = {
		0,1,2,
		3,0,4,
		5,0,6,
		3,6,0,
		0,2,4,
		5,1,0,
		2,1,5,
		7,6,3,
		6,7,5,
		7,3,4,
		7,4,2,
		7,2,5
	};
}

void geneCUBE(std::vector<double>& out_vertices, std::vector<unsigned int>& out_triangles) {
	out_vertices = {
		-2.0f,-1.0f,-1.0f,  //0
		-2.0f,-1.0f, 1.0f,  //1
		-2.0f, 1.0f, 1.0f,  //2
		 2.0f, 1.0f,-1.0f,  //3 
		-2.0f, 1.0f,-1.0f,  //4
		 2.0f,-1.0f, 1.0f,  //5
		 2.0f,-1.0f,-1.0f,  //6
		 2.0f, 1.0f, 1.0f   //7
	};

	out_triangles = {
		0,1,2,
		3,0,4,
		5,0,6,
		3,6,0,
		0,2,4,
		5,1,0,
		2,1,5,
		7,6,3,
		6,7,5,
		7,3,4,
		7,4,2,
		7,2,5
	};

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

std::vector<double> matrix_mutiple(const std::vector<double>& affine_m, const std::vector<double>& coord) {
	std::vector<double> result;
	for (int i = 0; i < 4; i += 1) {
		double tempt_sum = 0;
		for (int j = 0; j < 4; j += 1) {
			tempt_sum = tempt_sum + coord[j] * affine_m[i * 4 +j];
		}
		result.push_back(tempt_sum);
	}
	return result;
}

std::vector<double> matrix_m_matrix(const std::vector<double>& matrix_a, const std::vector<double>& matrix_b) {
	std::vector<double> result;
	for (int i = 0; i < 4; i += 1) {
		for (int j = 0; j < 4; j += 1) {
			double tempt = 0.f;
			for (int k = 0; k < 4; k += 1) {
				tempt = tempt + matrix_a[i * 4 + k] * matrix_b[k * 4 + j];
			}
			result.push_back(tempt);
		}
	}
	return result;
}

std::vector<double> a_rotate_degree(const std::vector<double>& coordinate, int degree, double length) {
	std::vector<double> transformed_coordinate;
	//affine matrix Z rotation 
	std::vector<double> affine_tran_a = { 1,0,0, length/2,
							   0,1,0,0,
							   0,0,1,0,
							   0,0,0,1 };

	std::vector<double> affineZ_a = { cos(degree * 3.14159 / 180),-sin(degree * 3.14159 / 180),0,0,
							   sin(degree * 3.14159 / 180), cos(degree * 3.14159 / 180),0, 0 ,
							  0, 0, 0, 0 ,
							   0, 0, 0, 1 };

	std::vector<double> tempt_affine;
	tempt_affine = matrix_m_matrix(affineZ_a, affine_tran_a);

	transformed_coordinate = matrix_mutiple(tempt_affine, coordinate);
	return transformed_coordinate;
}

std::vector<double> b_rotate_degree(const std::vector<double>& coordinate, int degree_a,int degree_b, double length_a, double length_b) {
	std::vector<double> transformed_coordinate;
	//affine matrix Z rotation 
	std::vector<double> affine_tran_a = { 1,0,0, length_a,
							   0,1,0,0,
							   0,0,1,0,
							   0,0,0,1 };

	std::vector<double> affine_tran_b = { 1,0,0, length_b/2,
						   0,1,0,0,
						   0,0,1,0,
						   0,0,0,1 };

	std::vector<double> affineZ_a = { cos(degree_a * 3.14159 / 180),-sin(degree_a * 3.14159 / 180),0,0,
							   sin(degree_a * 3.14159 / 180), cos(degree_a * 3.14159 / 180),0, 0 ,
							  0, 0, 0, 0 ,
							   0, 0, 0, 1 };

	std::vector<double> affineZ_b = { cos(degree_b * 3.14159 / 180),-sin(degree_b * 3.14159 / 180),0,0,
						   sin(degree_b * 3.14159 / 180), cos(degree_b * 3.14159 / 180),0, 0 ,
						  0, 0, 0, 0 ,
						   0, 0, 0, 1 };

	std::vector<double> tempt_affine;
	std::vector<double> tempt_affine2;
	std::vector<double> tempt_affine3;

	tempt_affine = matrix_m_matrix(affineZ_a, affine_tran_a);  // R2(theta1) * T (l1x)
	tempt_affine2 = matrix_m_matrix(tempt_affine, affineZ_b); // R2(theta1) * T (l1x) * R(theta2)
	tempt_affine3 = matrix_m_matrix(tempt_affine2, affine_tran_b);  // R2(theta1) * T (l1x) * R(theta2) * T(l2x/2)

	transformed_coordinate = matrix_mutiple(tempt_affine3, coordinate);
	return transformed_coordinate;
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

	std::vector<double> points_b;
	std::vector<unsigned int> indexs_b;
	//loadOBJ("bunny.obj", points, indexs);
	//geneCUBE(points, indexs);
	double cube_x = 6.0f;
	double cube_y = 2.0f;
	double cube_z = 2.0f;

	double cube_b_x = 4.0f;
	double cube_b_y = 2.0f;
	double cube_b_z = 2.0f;
	geneCUBE_length(cube_x, cube_y, cube_z,points, indexs);
	geneCUBE_length(cube_b_x, cube_b_y, cube_b_z, points_b, indexs_b);


  GLFWwindow* window;
  glfwSetErrorCallback(error_callback);
  if (!glfwInit())
    exit(EXIT_FAILURE);
  window = glfwCreateWindow(1080, 680, "Assignment_4_transform_arm", NULL, NULL);
  if (!window)
  {
    glfwTerminate();
    exit(EXIT_FAILURE);
  }

  glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, key_callback);
  int iframe = 0;
  int rotate_degree_a = -1;
  int rotate_degree_b = -1;

  bool r_reach_end = false;
  bool reach_end = false;

  glEnable(GL_DEPTH_TEST);
  // Accept fragment if it closer to the camera than the former one
  glDepthFunc(GL_LESS);

  double scale_ratio = 10.f;
  while (!glfwWindowShouldClose(window))
  {

	  //wiggle the arm repeatedly 
	  if ((rotate_degree_a < 60) & (r_reach_end == false)) {
		  rotate_degree_a += 1;
	  }
	  else {
		  r_reach_end = true;
	  }

	  if ((rotate_degree_a > 0) & (r_reach_end == true)) {
		  rotate_degree_a -= 1;
	  }
	  else {
		  r_reach_end = false;
	  }

	  if ((rotate_degree_b < 60) & (reach_end == false)) {
		  rotate_degree_b += 1;
	  }
	  else {
		  reach_end = true;
	  }

	  if ((rotate_degree_b > 0) & (reach_end == true)) {
		  rotate_degree_b -= 1;
	  }
	  else {
		  reach_end = false;
	  }

	  //affine matrix transformation 
	  //store the vertices list of rotated and transformed vertices
	  std::vector<double> transformed_points;
	  std::vector<double> transformed_points_b;
	  // calculation of cube A
	  for (int i = 0; i < points.size() / 3; i += 1) {
		  std::vector<double> tempt_coord;
		  std::vector<double> new_coord_final;

		  tempt_coord.push_back(points[i * 3 + 0]);
		  tempt_coord.push_back(points[i * 3 + 1]);
		  tempt_coord.push_back(points[i * 3 + 2]);
		  tempt_coord.push_back(1.f);

		  new_coord_final = a_rotate_degree(tempt_coord, rotate_degree_a, cube_x);

		  //put the rotated coordinate into the vertex list
		  for (int j = 0; j < new_coord_final.size() - 1 ; j += 1) {
			  transformed_points.push_back(new_coord_final[j]);
		  }
	  }
	  // calculation of cube B
	  for (int i = 0; i < points_b.size() / 3; i += 1) {
		  std::vector<double> tempt_coord_b;
		  std::vector<double> new_coord_b_final;

		  tempt_coord_b.push_back(points_b[i * 3 + 0]);
		  tempt_coord_b.push_back(points_b[i * 3 + 1]);
		  tempt_coord_b.push_back(points_b[i * 3 + 2]);
		  tempt_coord_b.push_back(1.f);

		  new_coord_b_final = b_rotate_degree(tempt_coord_b, rotate_degree_a, rotate_degree_b, cube_x,cube_b_x);

		  //put the rotated coordinate into the vertex list
		  for (int j = 0; j < new_coord_b_final.size() - 1; j += 1) {
			  transformed_points_b.push_back(new_coord_b_final[j]);
		  }
	  }

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

	// draw first cube
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
	// draw second cube 
	for (unsigned int i = 0; i < indexs_b.size() / 3; ++i) {

		unsigned int i0 = indexs_b[i * 3 + 0];
		unsigned int i1 = indexs_b[i * 3 + 1];
		unsigned int i2 = indexs_b[i * 3 + 2];
		glColor3f(1.0f, 0.7f, 1.0f);
		glVertex3f(transformed_points_b[i0 * 3 + 0], transformed_points_b[i0 * 3 + 1], transformed_points_b[i0 * 3 + 2]);
		glColor3f(0.3f, 0.2f, 0.7f);
		glVertex3f(transformed_points_b[i1 * 3 + 0], transformed_points_b[i1 * 3 + 1], transformed_points_b[i1 * 3 + 2]);
		glColor3f(0.1f, 0.2f, 0.4f);
		glVertex3f(transformed_points_b[i2 * 3 + 0], transformed_points_b[i2 * 3 + 1], transformed_points_b[i2 * 3 + 2]);
	}

		glEnd();
    glfwSwapBuffers(window);
    glfwPollEvents();
  }
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}