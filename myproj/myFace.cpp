#include "myFace.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <GL/glew.h>

myFace::myFace(void)
{
	adjacent_halfedge = NULL;
	normal = new myVector3D(1.0, 1.0, 1.0);
}

myFace::~myFace(void)
{
	if (normal) delete normal;
}

void myFace::computeNormal()
{
	/**** TODO ****/
	myPoint3D* p1 = this->adjacent_halfedge->source->point;
	myPoint3D* p2 = this->adjacent_halfedge->next->source->point;
	myPoint3D* p3 = this->adjacent_halfedge->prev->source->point;

	myVector3D v1 (p2->X - p1->X, p2->Y - p1->Y, p2->Z - p1->Z);
	myVector3D v2 (p3->X - p1->X, p3->Y - p1->Y, p3->Z - p1->Z);

	if (normal) delete normal;

	normal = &(v1.crossproduct(v2));
}
