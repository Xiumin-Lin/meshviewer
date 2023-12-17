#include "myVertex.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myFace.h"
#include <iostream>
using namespace std;

static int vertex_id_cpt = 1;

myVertex::myVertex(void)
{
	id = vertex_id_cpt++;
	index = id;
	//cout << "create vertex : " << id << endl;
	point = NULL;
	originof = NULL;
	normal = new myVector3D(1.0,1.0,1.0);
}

myVertex::~myVertex(void)
{
	if (normal) delete normal;
	if (point) delete point;
}

void myVertex::computeNormal()
{
	delete normal;
	normal = new myVector3D(0.0, 0.0, 0.0);
	myHalfedge* steph = originof;
	//cout << "do steph : " << steph->id << endl;
	do {
		if(steph->adjacent_face != NULL) *normal += *(steph->adjacent_face->normal);
		steph = steph->twin->next;
	} while (originof != steph);
	normal->normalize();
}
