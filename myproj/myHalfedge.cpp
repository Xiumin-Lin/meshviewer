#include "myHalfedge.h"
#include <iostream>

static int halfede_id_cpt = 0;

myHalfedge::myHalfedge(void)
{
	source = NULL; 
	adjacent_face = NULL; 
	next = NULL;  
	prev = NULL;  
	twin = NULL;  
	id = halfede_id_cpt++;
	//std::cout << "create halfedge : " << id << std::endl;
}

void myHalfedge::copy(myHalfedge *ie)
{
/**** TODO ****/
}

myHalfedge::~myHalfedge(void)
{
}
