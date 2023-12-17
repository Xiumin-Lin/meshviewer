#pragma once
#include "myFace.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <map>
#include <vector>
#include <string>

class myMesh
{
public:
	std::vector<myVertex *> vertices;
	std::vector<myHalfedge *> halfedges;
	std::vector<myFace *> faces;
	std::string name;

	void checkMesh();
	bool readFile(std::string filename);
	void computeNormals();
	void normalize();

	void subdivisionCatmullClark();

	void splitFaceTRIS(myFace*, myPoint3D*);
	void splitFaceTRIS(myFace *, myVertex*); // custom

	void splitEdge(myHalfedge*, myPoint3D*);
	void splitEdge(myHalfedge *, myVertex*);	// custom

	void splitFaceQUADS(myFace*, myPoint3D*);
	void splitFaceQUADS(myFace*, myVertex*, std::map<myHalfedge*, myVertex*>&); // custom

	void triangulate();
	bool triangulate(myFace *);

	void collapse(myHalfedge *);
	void cleanOverlappingMesh(myVertex* v);
	myHalfedge* collapseFace(myHalfedge *);
	myHalfedge* getShortestEdge();

	void clear();

	myMesh(void);
	~myMesh(void);
};

