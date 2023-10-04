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
	//cout << index << " : originof = " << endl;
	/**** TODO ****/
	/*delete normal;
	normal = new myVector3D(0.0, 0.0, 0.0);
	myHalfedge* nextVector = originof;
	int count = 0;
	do {
		*normal += *(nextVector->adjacent_face->normal);
		nextVector = nextVector->next->next->twin;
		count++;
	} while (nextVector != originof);

	*normal = *normal / count;*/

	delete normal;
	normal = new myVector3D(0.0, 0.0, 0.0);
	myHalfedge* h = originof;
	myHalfedge* steph = h;
	int count_outer = 0;

	do {
		*normal += *(steph->adjacent_face->normal);
		myHalfedge* steph_inner = steph;

		do {
			steph_inner = steph_inner->next;
		} while (steph_inner != steph);

		steph = steph->twin->next;
		count_outer++;
	} while (h != steph);
	//cout << index << " : originof = " << originof->index << " - normal = " << normal->dX << ", " << normal->dY << ", " << normal->dZ;
	*normal = *normal / count_outer;
	//cout << " / " << count_outer << " = " << normal->dX << ", " << normal->dY << ", " << normal->dZ << endl;
}
