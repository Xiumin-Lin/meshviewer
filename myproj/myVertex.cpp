#include "myVertex.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myFace.h"

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
	/**** TODO ****/
	delete normal;
	normal = new myVector3D(0.0, 0.0, 0.0);
	myHalfedge* nextVector = originof;
	int count = 0;
	do {
		*normal += *(nextVector->adjacent_face->normal);
		nextVector = nextVector->next->next->twin;
		count++;
	} while (nextVector != originof);

	*normal = *normal / count;
}
