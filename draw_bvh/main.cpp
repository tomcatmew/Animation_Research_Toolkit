#include <GLFW/glfw3.h>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>

bool loadBVH(std::vector<std::string>& aBone,
	std::vector<std::string>& aChannelRotTransBone,
	int& nframe,
	std::vector<double>& aValueRotTransBone,
	const std::string& path_bvh)
{
	std::ifstream fin;
	fin.open(path_bvh.c_str());
	if (!fin.is_open()) {
		std::cout << "cannot open file" << std::endl;
		return false;
	}
	aBone.clear();
	aChannelRotTransBone.clear();
	//
	std::string line;
	std::vector<int> stackIndBone;
	while (std::getline(fin, line))
	{
		if (line[line.size() - 1] == '\n') line.erase(line.size() - 1); // remove the newline code
		if (line[line.size() - 1] == '\r') line.erase(line.size() - 1); // remove the newline code

		std::vector<std::string> aToken;
		std::istringstream iss(line);
		for (std::string line; iss >> line; )
			aToken.push_back(line);

		//for (int i = 0; i < aToken.size(); i += 1) {
		//	std::cout << aToken[i];
		//}
		//std::cout << std::endl;

		if (aToken[0] == "HIERARCHY") {
			assert(aBone.empty());
		}
		else if (aToken[0] == "ROOT") {
			assert(aBone.size() == 0);
			CRigBone br;
			assert(aToken.size() == 2);
			br.name = aToken[1];
			aBone.push_back(br);
		}
		else if (aToken[0] == "{") {
			stackIndBone.push_back(aBone.size() - 1);
			if (stackIndBone.size() > 1) {
				int ibp = stackIndBone[stackIndBone.size() - 2];
				int ib = aBone.size() - 1;
				aBone[ib].ibone_parent = ibp;
			}
		}
		else if (aToken[0] == "}") {
			stackIndBone.resize(stackIndBone.size() - 1);
		}
		else if (aToken[0] == "OFFSET") {
			assert(aToken.size() == 4);
			int ib = aBone.size() - 1;
			double org_x = rig_v3q::myStod(aToken[1]);
			double org_y = rig_v3q::myStod(aToken[2]);
			double org_z = rig_v3q::myStod(aToken[3]);
			aBone[ib].invBindMat[3] = -org_x;
			aBone[ib].invBindMat[7] = -org_y;
			aBone[ib].invBindMat[11] = -org_z;
			if (stackIndBone.size() > 1) {
				const int ibp = stackIndBone[stackIndBone.size() - 2];
				assert(ibp < (int)aBone.size());
				aBone[ib].invBindMat[3] += aBone[ibp].invBindMat[3];
				aBone[ib].invBindMat[7] += aBone[ibp].invBindMat[7];
				aBone[ib].invBindMat[11] += aBone[ibp].invBindMat[11];
			}
		}
		else if (aToken[0] == "CHANNELS") {
			assert(aToken.size() >= 2);
			int nch = rig_v3q::myStoi(aToken[1]);
			assert((int)aToken.size() == nch + 2);
			assert(!aBone.empty());
			const std::size_t ib = aBone.size() - 1;
			for (int ich = 0; ich < nch; ++ich) {
				const std::string& type_ch = aToken[ich + 2];
				if (type_ch == "Xposition") { aChannelRotTransBone.emplace_back(ib, 0, false); }
				else if (type_ch == "Yposition") { aChannelRotTransBone.emplace_back(ib, 1, false); }
				else if (type_ch == "Zposition") { aChannelRotTransBone.emplace_back(ib, 2, false); }
				else if (type_ch == "Xrotation") { aChannelRotTransBone.emplace_back(ib, 0, true); }
				else if (type_ch == "Yrotation") { aChannelRotTransBone.emplace_back(ib, 1, true); }
				else if (type_ch == "Zrotation") { aChannelRotTransBone.emplace_back(ib, 2, true); }
				else {
					std::cout << "ERROR-->undefiend type" << std::endl;
				}
			}
		}
		else if (aToken[0] == "JOINT") {
			CRigBone br;
			assert(aToken.size() == 2);
			br.name = aToken[1];
			aBone.push_back(br);
		}
		else if (aToken[0] == "End") {
			assert(aToken[1] == "Site");
			CRigBone br;
			assert(aToken.size() == 2);
			br.name = aToken[1];
			aBone.push_back(br);
		}
		else if (aToken[0] == "MOTION") {
			break;
		}
	}
	nframe = 0;
	{
		std::string stmp0;
		{
			std::getline(fin, line);
			std::stringstream ss(line);
			ss >> stmp0 >> nframe;
			//      std::cout << "frame: " << nframe << std::endl;
		}
		std::getline(fin, line);
		//    std::cout << "frametime: " << line << std::endl;
	}
	const int nchannel = aChannelRotTransBone.size();
	aValueRotTransBone.resize(nframe* nchannel);
	for (int iframe = 0; iframe < nframe; ++iframe) {
		std::getline(fin, line);
		line = rig_v3q::MyReplace(line, '\t', ' ');
		if (line[line.size() - 1] == '\n') line.erase(line.size() - 1); // remove the newline code
		if (line[line.size() - 1] == '\r') line.erase(line.size() - 1); // remove the newline code
		std::vector<std::string> aToken = rig_v3q::MySplit(line, ' ');
		//    std::cout << aToken.size() << " " << aChannelRotTransBone.size() << std::endl;
		assert(aToken.size() == aChannelRotTransBone.size());
		for (int ich = 0; ich < nchannel; ++ich) {
			aValueRotTransBone[iframe * nchannel + ich] = rig_v3q::myStod(aToken[ich]);
		}
	}
	// ---------------
	for (std::size_t ibone = 0; ibone < aBone.size(); ++ibone) {
		CRigBone& bone = aBone[ibone];
		bone.scale = 1.0;
		bone.quatRelativeRot[0] = 1.0;
		bone.quatRelativeRot[1] = 0.0;
		bone.quatRelativeRot[2] = 0.0;
		bone.quatRelativeRot[3] = 0.0;
		bone.transRelative[0] = 0.0;
		bone.transRelative[1] = 0.0;
		bone.transRelative[2] = 0.0;
		if (bone.ibone_parent != -1) {
			const CRigBone& bone_p = aBone[bone.ibone_parent];
			bone.transRelative[0] = (-bone.invBindMat[3]) - (-bone_p.invBindMat[3]);
			bone.transRelative[1] = (-bone.invBindMat[7]) - (-bone_p.invBindMat[7]);
			bone.transRelative[2] = (-bone.invBindMat[11]) - (-bone_p.invBindMat[11]);
		}
	}
	for (auto& bone : aBone) {
		for (int i = 0; i < 16; ++i) { bone.affmat3Global[i] = bone.invBindMat[i]; }
		int info; rig_v3q::CalcInvMat(bone.affmat3Global, 4, info);
	}
	return true;
}


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

	std::vector<double> affineX = { 1,0,0,0,
						   0, cos(20 * 3.14159 / 180),-sin(20 * 3.14159 / 180), 0 ,
						  0, sin(20 * 3.14159 / 180), cos(20 * 3.14159 / 180), 0 ,
						   0, 0, 0, 1 };

	std::vector<double> affineY = { cos(45 * 3.14159 / 180),0,sin(45 * 3.14159 / 180),0,
								 0, 1, 0, 0 ,
								 -sin(45 * 3.14159 / 180), 0,cos(45 * 3.14159 / 180), 0 ,
								 0, 0, 0, 1 };

	std::vector<double> tempt_affine;
	std::vector<double> tempt_affine2;

	tempt_affine = matrix_m_matrix(affineZ_a, affine_tran_a);
	//tempt_affine2 = matrix_m_matrix(tempt_affine, affineY);

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

	std::vector<double> affineX = { 1,0,0,0,
							   0, cos(20 * 3.14159 / 180),-sin(20 * 3.14159 / 180), 0 ,
							  0, sin(20 * 3.14159 / 180), cos(20 * 3.14159 / 180), 0 ,
							   0, 0, 0, 1 };

	std::vector<double> affineY = { cos(45 * 3.14159 / 180),0,sin(45 * 3.14159 / 180),0,
							 0, 1, 0, 0 ,
							 -sin(45 * 3.14159 / 180), 0,cos(45 * 3.14159 / 180), 0 ,
							 0, 0, 0, 1 };

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

	//std::vector<double> tempt_affine4;

	tempt_affine = matrix_m_matrix(affineZ_a, affine_tran_a);  // R2(theta1) * T (l1x)
	tempt_affine2 = matrix_m_matrix(tempt_affine, affineZ_b); // R2(theta1) * T (l1x) * R(theta2)
	tempt_affine3 = matrix_m_matrix(tempt_affine2, affine_tran_b);  // R2(theta1) * T (l1x) * R(theta2) * T(l2x/2)

	//tempt_affine4 = matrix_m_matrix(tempt_affine3, affineY);

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

  GLFWwindow* window;
  glfwSetErrorCallback(error_callback);
  if (!glfwInit())
    exit(EXIT_FAILURE);
  window = glfwCreateWindow(1080, 680, "Assignment_5_read_BVH", NULL, NULL);
  if (!window)
  {
    glfwTerminate();
    exit(EXIT_FAILURE);
  }

  glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, key_callback);

  glEnable(GL_DEPTH_TEST);
  // Accept fragment if it closer to the camera than the former one
  glDepthFunc(GL_LESS);

  std::vector<std::string> aBone;
  std::vector<std::string> aChannelRotTransBone;
  int nframe = 0;
  std::vector<double> aValRotTransBone;

  std::string path_bvh = "02_01.bvh";

  loadBVH(aBone, aChannelRotTransBone, nframe, aValRotTransBone,path_bvh);

  double scale_ratio = 10.f;
  while (!glfwWindowShouldClose(window))
  {

    float ratio;
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    ratio = width / (float) height;

    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity(); 
    glOrtho(-ratio* scale_ratio, ratio * scale_ratio, -1.f * scale_ratio, 1.f * scale_ratio, 1.f * scale_ratio, -1.f * scale_ratio);
	//gluPerspective(120,ratio,0.001f,1.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    //glRotatef((float) glfwGetTime() * 50.f, 0.f, 0.f, 1.f);
	glBegin(GL_TRIANGLES);
	glColor3f(0.7f, 1.f, 0.f);
	glVertex3f(-0.6f, -0.4f, 0.f);

	glColor3f(0.7f, 1.f, 0.f);
	glVertex3f(0.6f, -0.4f, 0.f);

	glColor3f(0.1f, 1.f, 1.f);
	glVertex3f(0.0f, 0.8f, 0.f);


	glEnd();
    glfwSwapBuffers(window);
    glfwPollEvents();
  }
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}