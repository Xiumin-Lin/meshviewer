#include "myVertex.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myFace.h"
#include <iostream>
using namespace std;

myVertex::myVertex(void)
{
	point = NULL;
	originof = NULL;
	normal = new myVector3D(1.0,1.0,1.0);
}

myVertex::~myVertex(void)
{
	if (normal) delete normal;
}

void myVertex::computeNormal()
{
	delete normal;
	normal = new myVector3D(0.0, 0.0, 0.0);
	myHalfedge* steph = originof;

	do {
		*normal += *(steph->adjacent_face->normal);
		steph = steph->twin->next;
	} while (originof != steph);
	normal->normalize();
}
