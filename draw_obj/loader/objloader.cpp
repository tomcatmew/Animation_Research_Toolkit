#include <vector>
#include <stdio.h>
#include <string>
#include <cstring>
#include <iostream>

#include <glm/glm.hpp>

#include "objloader.hpp"

//loadOBJ function to load simple obj file 
bool loadOBJ(
	const char* path,
	std::vector<glm::vec3>& out_vertices,
	std::vector<unsigned int>& out_triangles

){
	printf("Loading OBJ file %s...\n", path);

	std::vector<glm::vec3> temp_vertices; 
	std::vector<unsigned int> vertexIndices;

	FILE * file = fopen(path, "r");
	if( file == NULL ){
		printf("cannot open the file !\n");
		getchar();
		return false;
	}

	while( 1 ){

		char lineHeader[128];
		// read the first word 
		int res = fscanf(file, "%s", lineHeader);
		if (res == EOF)
			break; // EOF = End Of File Quit 
		// else : parse lineHeader
		
		if ( strcmp( lineHeader, "v" ) == 0 ){
			glm::vec3 vertex;
			fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z );
			temp_vertices.push_back(vertex);
		}
		else if ( strcmp( lineHeader, "f" ) == 0 ){
			unsigned int vertexIndex[3];
			int matches = fscanf(file, "%d %d %d/\n", &vertexIndex[0], &vertexIndex[1], &vertexIndex[2]);
			if (matches != 3){
				printf("File can't be read");
				fclose(file);
				return false;
			}
			vertexIndices.push_back(vertexIndex[0]);
			vertexIndices.push_back(vertexIndex[1]);
			vertexIndices.push_back(vertexIndex[2]);
		}else{
			// comment
			char stupidBuffer[1000];
			fgets(stupidBuffer, 1000, file);
		}

	}
	std::cout << "print vertex\n";
	std::cout << temp_vertices.size();
	std::cout << "\n";
	for ( int i = 0; i < temp_vertices.size(); i++) {
		for (int j = 0; j < 3; j++) {
			std::cout << temp_vertices[i][j];
			std::cout << "\n";
		}

	}
	std::cout << "print indices size is \n";
	std::cout << vertexIndices.size()/3;
	std::cout << "\n";
	for ( int i = 0; i < vertexIndices.size(); i++) {
		std::cout << vertexIndices[i] << '\n ';
		std::cout << "\n";
	}
	// For each vertex of each triangle
	for( unsigned int i=0; i<vertexIndices.size(); i++ ){
		// Get the indices of its attributes
		// Get the attributes to the index
		glm::vec3 vertex = temp_vertices[vertexIndices[i] -1 ];
		// Put the attributes in buffers
		out_vertices.push_back(vertex);
	}
	out_triangles = vertexIndices;
	fclose(file);
	return true;
}

