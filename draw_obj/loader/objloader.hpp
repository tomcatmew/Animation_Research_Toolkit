#ifndef OBJLOADER_H
#define OBJLOADER_H

bool loadOBJ(
	const char * path, 
	//std::vector<glm::vec3> & out_vertices,
	std::vector<double>& out_vertices,
	std::vector<unsigned int>& out_triangles

);



bool loadAssImp(
	const char * path, 
	std::vector<unsigned short> & indices,
	std::vector<glm::vec3> & vertices
);

#endif