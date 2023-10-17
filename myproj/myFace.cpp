#include "myFace.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <GL/glew.h>

static int face_id_cpt = 0;

myFace::myFace(void)
{
	adjacent_halfedge = NULL;
	normal = new myVector3D(1.0, 1.0, 1.0);
	index = face_id_cpt++;
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

	normal->setNormal(p1, p2, p3);
}
