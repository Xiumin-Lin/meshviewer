#include "myFace.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <GL/glew.h>
#include <iostream>
using namespace std;

static int face_id_cpt = 0;

myFace::myFace(void)
{
	adjacent_halfedge = NULL;
	normal = new myVector3D(1.0, 1.0, 1.0);
	id = face_id_cpt++;
	index = id;
	//cout << "create face : " << id << endl;
}

myFace::~myFace(void)
{
	//cout << "delete face : " << index << endl;
	if (normal != NULL || normal != nullptr) delete normal;
}

void myFace::computeNormal()
{
	/**** TODO ****/
	myPoint3D* p1 = this->adjacent_halfedge->source->point;
	myPoint3D* p2 = this->adjacent_halfedge->next->source->point;
	myPoint3D* p3 = this->adjacent_halfedge->prev->source->point;

	normal->setNormal(p1, p2, p3);
}
