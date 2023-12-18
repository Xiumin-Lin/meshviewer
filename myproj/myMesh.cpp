#include "myMesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <GL/glew.h>
#include <algorithm>
#include "myvector3d.h"

using namespace std;

myMesh::myMesh(void)
{
	/**** TODO ****/
}

myMesh::~myMesh(void)
{
	/**** TODO ****/
	for (myVertex* v : vertices)	if (v) delete v;
	for (myHalfedge* h : halfedges)	if (h) delete h;
	for (myFace* f : faces)			if (f) delete f;
}

void myMesh::clear()
{
	for (unsigned int i = 0; i < vertices.size(); i++)	if (vertices[i] != NULL)	delete vertices[i];
	for (unsigned int i = 0; i < halfedges.size(); i++) if (halfedges[i] != NULL)	delete halfedges[i];
	for (unsigned int i = 0; i < faces.size(); i++)		if (faces[i] != NULL)		delete faces[i];

	vector<myVertex *> empty_vertices;    vertices.swap(empty_vertices);
	vector<myHalfedge *> empty_halfedges; halfedges.swap(empty_halfedges);
	vector<myFace *> empty_faces;         faces.swap(empty_faces);
}

/* Test function to assert every mesh is correct */
void myMesh::checkMesh()
{
	// check if each halfedge has a next, prev and twin
	bool allHasTwin = true;
	bool allHasNext = true;
	bool allHasPrev = true;
	for (myHalfedge* he : halfedges)
	{
		if (allHasNext) {
			if (he->next == NULL || he->next == nullptr)
			{
				std::cout << "Error! Not all halfedge have their next!" << std::endl;
				allHasNext = false;
			}
			else if (he->next->prev != he) {
				std::cout << "Error! halfedge->next->prev should be the halfedge itself" << std::endl;
				allHasNext = false;
			}
		}
		if (allHasPrev) {
			if (he->prev == NULL || he->prev == nullptr)
			{
				std::cout << "Error! Not all halfedge have their prev!" << std::endl;
				allHasPrev = false;
			}
			else if (he->prev->next != he) {
				std::cout << "Error! halfedge->prev->next should be the halfedge itself!" << std::endl;
				allHasPrev = false;
			}
		}
		if (allHasTwin) {
			if (he->twin == NULL || he->twin == nullptr)
			{
				std::cout << "Error! Not all halfedge have their twins!" << std::endl;
				allHasTwin = false;
			}
			else if (he->twin->twin != he) {
				std::cout << "Error! Twin of halfedge is not the halfedge itself!" << std::endl;
				allHasTwin = false;
			}
			else if (he->next == NULL || he->next == nullptr || he->twin->source != he->next->source) {
				std::cout << "Error! Twin of halfedge has not the same source as the next of halfedge!" << std::endl;
				allHasTwin = false;
			}
		}
	}
	if (allHasNext)	std::cout << "Each halfedge has a next!" << std::endl;
	if (allHasPrev) std::cout << "Each halfedge has a prev!" << std::endl;
	if (allHasTwin) std::cout << "Each halfedge has a twin!" << std::endl;

	// check if each face can move around
	bool allFaceCanMoveAround = true;
	for (myFace* fa : faces)
	{
		map<int, bool> visited_vertex;
		myHalfedge* step_h = fa->adjacent_halfedge;
		int cpt_v = 0;
		do {
			if (visited_vertex.find(step_h->id) != visited_vertex.end()) {
				std::cout << "Error! Face loop indefinitely for : " << fa->id << std::endl;
				allFaceCanMoveAround = false;
				break;
			}
			cpt_v++;
			visited_vertex[step_h->id] = true;
			step_h = step_h->next;
		} while (step_h != fa->adjacent_halfedge);

		if (cpt_v < 3) {
			std::cout << "Error! Face have less than 3 vertices : " << fa->id << std::endl;
		}
	}
	if (allFaceCanMoveAround) std::cout << "Each face can move around!" << std::endl;

	// check if each vertex can move around
	bool allVertexCanMoveAround = true;
	for (myVertex* v : vertices)
	{
		map<int, bool> visited_faces;
		myHalfedge* step_h = v->originof;
		int cpt_v = 0;
		do {
			if (step_h->adjacent_face != NULL) {
				if (visited_faces.find(step_h->adjacent_face->id) != visited_faces.end()) {
					std::cout << "Error! Vertex loop indefinitely for " << v->id << std::endl;
					allVertexCanMoveAround = false;
					break;
				}
				visited_faces[step_h->adjacent_face->id] = true;
			}
			cpt_v++;
			step_h = step_h->twin->next;
		} while (step_h != v->originof);

		if (cpt_v < 2) {
			std::cout << "Error! Vertex only have one edge for " << v->id << std::endl;
		}

		if (allVertexCanMoveAround) break;
	}
	if (allVertexCanMoveAround) std::cout << "Each vertex can move around!" << std::endl;

	std::cout << std::endl;
}

bool myMesh::readFile(std::string filename)
{
	string s, t, u;
	//vector<int> faceids;
	//myHalfedge **hedges;

	ifstream fin(filename);
	if (!fin.is_open()) {
		std::cout << "Unable to open file!\n";
		return false;
	}
	name = filename;

	map<pair<int, int>, myHalfedge *> twin_map;
	map<pair<int, int>, myHalfedge *>::iterator it;

	while (getline(fin, s))
	{
		stringstream myline(s);
		myline >> t;
		if (t == "g") {}
		else if (t == "v")
		{
			float x, y, z;
			myline >> x >> y >> z;
			//std::cout << "v " << x << " " << y << " " << z << endl;

			myVertex* vertex = new myVertex();
			vertex->point = new myPoint3D(x, y, z);
			vertices.push_back(vertex);
		}
		else if (t == "mtllib") {}
		else if (t == "usemtl") {}
		else if (t == "s") {}
		else if (t == "f")
		{
			vector<int> list;
			int vSize = vertices.size();

			//std::cout << "f";
			while (myline >> u)
			{
				int vertexIdx = atoi((u.substr(0, u.find("/"))).c_str());
				vertexIdx = vertexIdx > 0 ? vertexIdx - 1 : (vertexIdx < 0 ? (vertexIdx % vSize) + vSize : 0);
				//std::cout << " " << vertexIdx;
				list.push_back(vertexIdx);
			}
			//std::cout << endl;

			int nbVertex = list.size();
			myFace* face = new myFace();
			myHalfedge* prevHalfedge = nullptr;
			for (int i = 0; i < nbVertex; i++)
			{
				int vertex_a = list[i];
				int vertex_b = list[(i + 1) % nbVertex];
				myHalfedge* e = new myHalfedge();

				/* neighbours setting */
				e->source = this->vertices.at(vertex_a);
				e->source->originof = e;
				e->adjacent_face = face;
				e->prev = prevHalfedge;

				if (i == 0) face->adjacent_halfedge = e;	// first halfedge
				else prevHalfedge->next = e;				// between first and last halfedge
				if (i == nbVertex - 1) {					// last halfedge, don't move up this line
					e->next = face->adjacent_halfedge;
					face->adjacent_halfedge->prev = e;
				}
				prevHalfedge = e;							// don't move up this line

				/* twin setting */
				it = twin_map.find(make_pair(vertex_b, vertex_a));
				if (it == twin_map.end()) { 
					twin_map[make_pair(vertex_a, vertex_b)] = e;
				}
				else { 
					e->twin = it->second;
					e->twin->twin = e;
					twin_map.erase(it);	// erase the twin from the map so we can create a twin for every halfedge still in the map at the end
				}

				halfedges.push_back(e);
			}
			faces.push_back(face);
		}
	}

	// create twin for every halfedge still in the map, it's considerate as a boundary halfedge
	for (auto& kv: twin_map)
	{
		pair<int, int> key = kv.first;
		int vertex_a = key.first;
		int vertex_b = key.second;
		myHalfedge* twin = kv.second;	// halfedge a->b
		myHalfedge* e = new myHalfedge();
		/* neighbours setting */
		e->source = this->vertices.at(vertex_b);
		e->source->originof = e;
		e->adjacent_face = NULL;	// no face for boundary halfedge
		e->twin = twin;
		twin->twin = e;
		/* setup prev */
		myHalfedge* steph = twin->next;
		while (steph->twin != NULL && steph->twin->adjacent_face != NULL)
		{
			steph = steph->twin->next;
		}
		if (steph->twin != NULL)
		{
			e->prev = steph->twin;
			e->prev->next = e;
		}
		/* setup prev */
		steph = twin->prev;
		while (steph->twin != NULL && steph->twin->adjacent_face != NULL)
		{
			steph = steph->twin->prev;
		}
		if (steph->twin != NULL)
		{
			e->next = steph->twin;
			e->next->prev = e;
		}
		
		halfedges.push_back(e);
	}

	checkMesh();
	normalize();
	return true;
}

void myMesh::computeNormals()
{
	/**** TODO ****/
	//std::cout << "-----------------start faces computeNormals" << endl;
	for (size_t i = 0; i < faces.size(); i++)
	{
		//std::cout << i << " : " << faces[i]->id << endl;
		faces[i]->computeNormal();
	}
	//std::cout << "-----------------start vertices computeNormals" << endl;
	for (size_t i = 0; i < vertices.size(); i++)
	{
		//std::cout << i << " : " << vertices[i]->id << endl;
		vertices[i]->computeNormal();
	}
	//std::cout << "-----------------end vertices computeNormals" << endl;
}

void myMesh::normalize()
{
	if (vertices.size() < 1) return;

	int tmpxmin = 0, tmpymin = 0, tmpzmin = 0, tmpxmax = 0, tmpymax = 0, tmpzmax = 0;

	for (unsigned int i = 0; i < vertices.size(); i++) {
		if (vertices[i]->point->X < vertices[tmpxmin]->point->X) tmpxmin = i;
		if (vertices[i]->point->X > vertices[tmpxmax]->point->X) tmpxmax = i;

		if (vertices[i]->point->Y < vertices[tmpymin]->point->Y) tmpymin = i;
		if (vertices[i]->point->Y > vertices[tmpymax]->point->Y) tmpymax = i;

		if (vertices[i]->point->Z < vertices[tmpzmin]->point->Z) tmpzmin = i;
		if (vertices[i]->point->Z > vertices[tmpzmax]->point->Z) tmpzmax = i;
	}

	double xmin = vertices[tmpxmin]->point->X, xmax = vertices[tmpxmax]->point->X,
		ymin = vertices[tmpymin]->point->Y, ymax = vertices[tmpymax]->point->Y,
		zmin = vertices[tmpzmin]->point->Z, zmax = vertices[tmpzmax]->point->Z;

	double scale = (xmax - xmin) > (ymax - ymin) ? (xmax - xmin) : (ymax - ymin);
	scale = scale > (zmax - zmin) ? scale : (zmax - zmin);

	for (unsigned int i = 0; i < vertices.size(); i++) {
		vertices[i]->point->X -= (xmax + xmin) / 2;
		vertices[i]->point->Y -= (ymax + ymin) / 2;
		vertices[i]->point->Z -= (zmax + zmin) / 2;

		vertices[i]->point->X /= scale;
		vertices[i]->point->Y /= scale;
		vertices[i]->point->Z /= scale;
	}
}

void myMesh::splitFaceTRIS(myFace* f, myPoint3D* p) { /**** TODO ****/ }
void myMesh::splitFaceTRIS(myFace* f, myVertex* v) { // custom
	/**** TODO ****/
}

void myMesh::splitEdge(myHalfedge* e1, myPoint3D* p) { /**** TODO ****/ }
void myMesh::splitEdge(myHalfedge* e1, myVertex* v) { // custom
	/**** TODO ****/
	myHalfedge* etwin = e1->twin;

	// create new halfedges
	myHalfedge* new_he = new myHalfedge();
	myHalfedge* new_etwin = new myHalfedge();
	new_he->source = v;
	new_etwin->source = v;
	new_he->source->originof = new_he;

	// link twins
	new_he->twin = e1->twin;
	new_he->twin->twin = new_he;
	new_etwin->twin = e1;
	new_etwin->twin->twin = new_etwin;

	new_he->adjacent_face = e1->adjacent_face;
	new_etwin->adjacent_face = etwin->adjacent_face;

	// link halfedges
	new_he->prev = e1;
	new_he->next = e1->next;
	new_he->prev->next = new_he;
	new_he->next->prev = new_he;
	
	new_etwin->prev = etwin;
	new_etwin->next = etwin->next;
	new_etwin->prev->next = new_etwin;
	new_etwin->next->prev = new_etwin;

	halfedges.push_back(new_he);
	halfedges.push_back(new_etwin);
}

void myMesh::splitFaceQUADS(myFace* f, myPoint3D* p) { /**** TODO ****/ }
void myMesh::splitFaceQUADS(myFace* f, myVertex* center_v, map<myHalfedge*, myVertex*>& halfedgesPointsMap) { // custom
	/**** TODO ****/
	vector<myHalfedge*> faceHalfedgeList;
	vector<myVertex*> middleVertexList;
	myHalfedge* step_h = f->adjacent_halfedge;
	// get halfedges and middle vertex of the face
	do
	{
		if (halfedgesPointsMap.find(step_h) != halfedgesPointsMap.end()) {
			faceHalfedgeList.push_back(step_h);
			middleVertexList.push_back(halfedgesPointsMap[step_h]);
		}
		else if (halfedgesPointsMap.find(step_h->next->twin) != halfedgesPointsMap.end())  {
			faceHalfedgeList.push_back(step_h);	// do not change step_h
			middleVertexList.push_back(halfedgesPointsMap[step_h->next->twin]);
		}
		step_h = step_h->next;
	} while (f->adjacent_halfedge != step_h);

	// for each middle vertex, create new halfedge, twins and link them
	size_t size = faceHalfedgeList.size();
	map<pair<int, int>, myHalfedge*> face_edge_map;
	for (size_t i = 0; i < size; i++)
	{
		myHalfedge* edge = faceHalfedgeList[i];
		myVertex* middle_v = middleVertexList[i];
		
		// create new halfedges
		myHalfedge* new_e = new myHalfedge();		// center to edge
		myHalfedge* new_twin = new myHalfedge();	// edge to center
		halfedges.push_back(new_e);
		halfedges.push_back(new_twin);

		new_e->source = center_v;
		new_twin->source = middle_v;
		new_e->source->originof = new_e;

		// link twins
		new_e->twin = new_twin;
		new_twin->twin = new_e;

		// link halfedges
		new_e->next = edge->next;
		new_e->next->prev = new_e;

		new_twin->prev = edge;
		new_twin->prev->next = new_twin;
		
		// map to retrieve halfedge from center_v id and middle_v id
		int prevMiddleId = (i - 1 < 0) ? middleVertexList[(i - 1 + size) % size]->id : middleVertexList[(i - 1) % size]->id;
		auto prev_it = face_edge_map.find(make_pair(center_v->id, prevMiddleId));
		if (prev_it == face_edge_map.end()) {
			face_edge_map[make_pair(middle_v->id, center_v->id)] = new_twin;
		}
		else {
			new_twin->next = prev_it->second;
			new_twin->next->prev = new_twin;
		}

		int nextMiddleId = middleVertexList[(i + 1) % size]->id;
		auto next_it = face_edge_map.find(make_pair(nextMiddleId, center_v->id));
		if (next_it == face_edge_map.end()) {
			face_edge_map[make_pair(center_v->id, middle_v->id)] = new_e;
		}
		else {
			new_e->prev = next_it->second;
			new_e->prev->next = new_e;
		}
	}
	// create new faces and link them
	myFace* new_little_face = f;
	for (size_t i = 0; i < faceHalfedgeList.size(); i++)
	{
		myHalfedge* edge = faceHalfedgeList[i];
		// create face (use the initial face first before create new one)
		if (i > 0) {
			new_little_face = new myFace();
			faces.push_back(new_little_face);
		}
		new_little_face->adjacent_halfedge = edge;

		// setup adjacent_halfedge for each halfedge of the little face
		myHalfedge* step_h = edge;
		do
		{
			step_h->adjacent_face = new_little_face;
			step_h = step_h->next;
		} while (step_h != edge);
	}
}

myVertex* computeFaceCenterPoint(myFace* face) {
	myPoint3D* center_p = new myPoint3D();

	double nb_vertex = 0;
	myHalfedge* step_h = face->adjacent_halfedge;
	do
	{
		*center_p += *(step_h->source->point);
		nb_vertex++;
		step_h = step_h->next;
	} while (face->adjacent_halfedge != step_h);
	*center_p /= nb_vertex;

	myVertex* center_vertex = new myVertex();
	center_vertex->point = center_p;
	return center_vertex;
}

myVertex* computeEdgeMidPoint(myHalfedge* e) {
	myVertex* middle_e = new myVertex();
	middle_e->point = new myPoint3D();
	*(middle_e->point) = (*(e->source->point) + *(e->twin->source->point)) / 2;
	return middle_e;
}

// https://www.youtube.com/watch?v=mfp1Z1mBClc
void myMesh::subdivisionCatmullClark()
{
	vector<myVertex*>			newVertices;
	map<myFace*, myVertex*>		facePointsMap;
	map<myHalfedge*, myVertex*>	halfedgesPointsMap;
	
	// Step 1: Create new vertices as face center point
	for (myFace* face : faces) {
		myVertex* center_v = computeFaceCenterPoint(face);
		facePointsMap[face] = center_v;
		newVertices.push_back(center_v);
	}

	// Step 2: Create new vertices as edge middle point
	for (myHalfedge* edge : halfedges) {
		// skip if already computed
		auto already_computed_it = halfedgesPointsMap.find(edge->twin);
		if (already_computed_it != halfedgesPointsMap.end()) continue;
		// compute edge middle point
		myVertex* middle_v = computeEdgeMidPoint(edge);
		halfedgesPointsMap[edge] = middle_v;
		newVertices.push_back(middle_v);
	}

	// Step 3: Move newEdgepoints
	map<myVertex*, myPoint3D*>	vertexNewPositionMap;
	for (auto& kv : halfedgesPointsMap) {
		myHalfedge* edge = kv.first;
		myVertex* middle_v = kv.second;

		// get points
		myPoint3D originP_1 = *edge->source->point;
		myPoint3D originP_2 = *edge->twin->source->point;
		myPoint3D newFaceCenterP_1;
		auto it = facePointsMap.find(edge->adjacent_face);
		if (it != facePointsMap.end()) newFaceCenterP_1 += *(it->second->point);

		myPoint3D newFaceCenterP_2;
		it = facePointsMap.find(edge->twin->adjacent_face);
		if (it != facePointsMap.end()) newFaceCenterP_2 += *(it->second->point);
		
		// compute new position
		myPoint3D* e = new myPoint3D();
		*e = (originP_1 + originP_2 + newFaceCenterP_1 + newFaceCenterP_2) / 4;

		vertexNewPositionMap[middle_v] = e;
	}

	// Step 4: Move origins point
	for (myVertex* origin_v : vertices) {
		myPoint3D Q;
		myPoint3D R;
		myPoint3D S;
		S += *origin_v->point;
		int n = 0;
		myHalfedge* step_h = origin_v->originof;
		do
		{
			auto face_it = facePointsMap.find(step_h->adjacent_face);
			if (face_it != facePointsMap.end()) {
				Q += *(face_it->second->point);
			}

			auto hedge_it = halfedgesPointsMap.find(step_h);
			if (hedge_it != halfedgesPointsMap.end()) {
				R += *(hedge_it->second->point);
			}
			else {
				hedge_it = halfedgesPointsMap.find(step_h->twin);
				if (hedge_it != halfedgesPointsMap.end()) 
				{
					R += *(hedge_it->second->point);
				}
			}
			n++;
			step_h = step_h->twin->next;
		} while (step_h != origin_v->originof);
		Q /= n;
		R /= n;
		S *= (n - 3);
		S /= n;

		myPoint3D* newOriginP = new myPoint3D();
		*newOriginP = (Q + R * 2 + S) / n;
		vertexNewPositionMap[origin_v] = newOriginP;
	}

	// Step 5: Connection
	// Connect each middle point to the associated original edges
	for (auto& kv : halfedgesPointsMap) {
		myVertex* middle_v = kv.second;
		myHalfedge* edge = kv.first;
		splitEdge(edge, middle_v);
	}

	// Connect each new face point to the new edge points of all original edges defining the original face
	for (auto& kv : facePointsMap) {
		myVertex* middle_v = kv.second;
		myFace* face = kv.first;
		splitFaceQUADS(face, middle_v, halfedgesPointsMap);
	}

	// Step 6: Update vertices position
	for (auto& kv : vertexNewPositionMap) {
		myVertex* v = kv.first;
		v->point = kv.second;
	}

	// add new vertices to the mesh
	for (myVertex* v : newVertices) {
		vertices.push_back(v);
	}

	checkMesh();
}

/*void myMesh::subdivisionCatmullClark()
{
	map<myFace*, myVertex*> facePointsMap;
	map<myHalfedge*, myVertex*> halfedgesPointsMap;
	map<pair<int, int>, myHalfedge*> twin_map;

	vector<myVertex*> newVertices;
	vector<myFace*> newFaces;
	vector<myHalfedge*> newHalfedges;

	int sizeV = vertices.size();
	int sizeH = halfedges.size();
	int sizeF = faces.size();

	for (size_t idx = 0; idx < sizeF; idx++) {
		myPoint3D* centerFacePoint = new myPoint3D();
		myHalfedge* currentEdge = faces[idx]->adjacent_halfedge;
		int cpt = 0;
		do
		{
			*centerFacePoint += *(currentEdge->source->point);
			currentEdge = currentEdge->next;
			cpt++;
		} while (currentEdge != faces[idx]->adjacent_halfedge);

		*centerFacePoint /= cpt;
		myVertex* newVertex = new myVertex();
		newVertex->point = centerFacePoint;
		newVertices.push_back(newVertex);
		facePointsMap[faces[idx]] = newVertex;
	}

	vector<myHalfedge*> listHalfedges;
	for (size_t idx = 0; idx < sizeH; idx++)
	{
		auto it = std::find(listHalfedges.begin(), listHalfedges.end(), halfedges[idx]->twin);
		if (it == listHalfedges.end()) {
			listHalfedges.push_back(halfedges[idx]);
		}
		else {
			halfedgesPointsMap[halfedges[idx]] = halfedgesPointsMap[halfedges[idx]->twin];// peut-être une erreur ici
			continue;
		}

		myPoint3D* facePoint1 = halfedges[idx]->source->point;
		myPoint3D* facePoint2 = halfedges[idx]->twin->source->point;
		myPoint3D* facePointMap1 = facePointsMap[halfedges[idx]->adjacent_face]->point;
		myPoint3D* facePointMap2 = facePointsMap[halfedges[idx]->twin->adjacent_face]->point;

		myPoint3D* edgePoint = new myPoint3D();
		*edgePoint += *facePoint1 + *facePoint2 + *facePointMap1 + *facePointMap2;
		*edgePoint /= 4;

		myVertex* newVertex = new myVertex();
		newVertex->point = edgePoint;

		newVertices.push_back(newVertex);
		halfedgesPointsMap[halfedges[idx]] = newVertex;
	}
	for (size_t idx = 0; idx < sizeV; idx++)
	{
		myPoint3D* avgFacePoint = new myPoint3D();
		myPoint3D* avgEdgePoint = new myPoint3D();

		myHalfedge* currentEdge = vertices[idx]->originof;
		myHalfedge* startEdge = currentEdge;

		int cpt = 0;

		do
		{
			*avgFacePoint += *facePointsMap[currentEdge->adjacent_face]->point;
			*avgEdgePoint += *halfedgesPointsMap[currentEdge]->point;
			currentEdge = currentEdge->twin->next;
			cpt++;
		} while (currentEdge != startEdge);

		*avgFacePoint /= cpt;
		*avgEdgePoint /= cpt;

		myPoint3D* newPoint = new myPoint3D();
		// ERREUR ICI SUREMENT
		*newPoint += (*avgFacePoint + *avgEdgePoint * 2.0 + *vertices[idx]->point * (cpt - 2)) / (cpt + 1); // Marche mais pas la formule
		//*newPoint += (*avgFacePoint + *avgEdgePoint * 2.0 + *vertices[idx]->point * (cpt - 3)) / cpt; // Marche mais n'est pas la formule


		vertices[idx]->point = newPoint;

	}
	for (size_t idx = 0; idx < sizeF; idx++) {
		myFace* oldFace = faces[idx];
		myHalfedge* currentEdge = oldFace->adjacent_halfedge;

		int numVertices = 0;

		do {
			numVertices++;
			currentEdge = currentEdge->next;
		} while (currentEdge != oldFace->adjacent_halfedge);
		currentEdge = oldFace->adjacent_halfedge;

		for (size_t i = 0; i < numVertices; i++) {
			myVertex* v1 = currentEdge->source;
			myVertex* v2 = halfedgesPointsMap[currentEdge];
			myVertex* v3 = facePointsMap[oldFace];
			myVertex* v4 = halfedgesPointsMap[currentEdge->prev];

			myFace* newFace = new myFace();
			newFaces.push_back(newFace);

			myHalfedge* e1 = new myHalfedge();
			myHalfedge* e2 = new myHalfedge();
			myHalfedge* e3 = new myHalfedge();
			myHalfedge* e4 = new myHalfedge();

			////////////// DEBUG - A SUPPRIMER //////////////
			e1->twin = e3;
			e2->twin = e4;
			e3->twin = e2;
			e4->twin = e1;
			/////////////////////////////////////////////////

			e1->next = e2;
			e2->next = e3;
			e3->next = e4;
			e4->next = e1;

			e1->prev = e4;
			e2->prev = e1;
			e3->prev = e2;
			e4->prev = e3;

			e1->source = v1;
			e2->source = v2;
			e3->source = v3;
			e4->source = v4;

			v1->originof = e1;
			v2->originof = e2;
			v3->originof = e3;
			v4->originof = e4;

			newFace->adjacent_halfedge = e1;

			e1->adjacent_face = newFace;
			e2->adjacent_face = newFace;
			e3->adjacent_face = newFace;
			e4->adjacent_face = newFace;

			faces.push_back(newFace);
			newHalfedges.push_back(e1);
			newHalfedges.push_back(e2);
			newHalfedges.push_back(e3);
			newHalfedges.push_back(e4);

			
			// NE MARCHE PAS
			associateTwin(twin_map, e1->source->index, e4->source->index, e2);
			associateTwin(twin_map, e2->source->index, e1->source->index, e3);
			associateTwin(twin_map, e3->source->index, e2->source->index, e4);
			associateTwin(twin_map, e4->source->index, e3->source->index, e1);
			
			currentEdge = currentEdge->prev;
		}
	}

	for (size_t idx = 0; idx < sizeF; idx++) {
		delete faces[idx];
	}
	faces.clear();
	faces = newFaces;

	int sizeNV = newVertices.size();
	for (size_t i = 0; i < sizeNV; i++)
	{
		vertices.push_back(newVertices[i]);
	}

	for (size_t idx = 0; idx < sizeH; idx++) {
		delete halfedges[idx];
	}
	halfedges.clear();
	halfedges = newHalfedges;

	checkMesh();
}
*/

void associateTwin(map<pair<int, int>, myHalfedge*>& twin_map, int idx_vertex_a, int idx_vertex_b, myHalfedge* e)
{
	map<pair<int, int>, myHalfedge*>::iterator it = twin_map.find(make_pair(idx_vertex_b, idx_vertex_a));
	if (it == twin_map.end()) {
		twin_map[make_pair(idx_vertex_a, idx_vertex_b)] = e;
	}
	else {
		e->twin = it->second;
		e->twin->twin = e;
	}
}

void myMesh::triangulate()
{
	int cpt = 0;
	int delete_cpt = 0;
	/**** TODO ****/
	for (size_t idx = 0; idx < faces.size(); idx++)
	{
		myFace* face = faces[idx];
		if (triangulate(face)) continue;
		myPoint3D* center_p = new myPoint3D();
		
		double vertex_cpt = 0;
		myHalfedge* step_halfedge = face->adjacent_halfedge;
		do
		{
			*center_p += *(step_halfedge->source->point);
			vertex_cpt++;
			step_halfedge = step_halfedge->next;
		} while (face->adjacent_halfedge != step_halfedge);

		if (vertex_cpt >= 4) {
			map<pair<int, int>, myHalfedge*> twin_map;

			// create center vertex
			*center_p /= vertex_cpt;
			myVertex* center_vertex = new myVertex();
			center_vertex->point = center_p;
			vertices.push_back(center_vertex);

			myFace* new_face = face;
			myHalfedge* originHalfedge = face->adjacent_halfedge;
			myHalfedge* stepHalfedge = face->adjacent_halfedge;
			int face_cpt = 0;
			do
			{
				// create face (use the initial face first before create new one)
				if (face_cpt > 0) {
					new_face = new myFace();
					faces.push_back(new_face);
				}
				face_cpt++;
				new_face->adjacent_halfedge = stepHalfedge;
				stepHalfedge->adjacent_face = new_face;

				// create halfedge
				myHalfedge* nextHalfedge = stepHalfedge->next;

				myHalfedge* nextHalf = new myHalfedge();
				nextHalf->source = nextHalfedge->source;
				nextHalf->source->originof = nextHalf;
				nextHalf->adjacent_face = new_face;

				myHalfedge* prevHalf = new myHalfedge();
				prevHalf->source = center_vertex;
				prevHalf->source->originof = prevHalf;
				prevHalf->adjacent_face = new_face;

				// associate halfedge
				nextHalf->prev = stepHalfedge;
				nextHalf->next = prevHalf;

				prevHalf->prev = nextHalf;
				prevHalf->next = stepHalfedge;

				stepHalfedge->prev = prevHalf;
				stepHalfedge->next = nextHalf;

				// create twin
				associateTwin(twin_map, nextHalf->source->id, prevHalf->source->id, nextHalf);
				associateTwin(twin_map, prevHalf->source->id, stepHalfedge->source->id, prevHalf);

				// save created halfedge
				halfedges.push_back(nextHalf);
				halfedges.push_back(prevHalf);

				// pass to next globale halfedge
				stepHalfedge = nextHalfedge;
			} while (originHalfedge != stepHalfedge);
		}
		else
		{
			std::cout << "Error when tryng triangulate face id " << face->id << " : have only " << vertex_cpt << " halfedges !" << endl;
		}
	}
	checkMesh();
}

//return false if already triangle, true othewise.
bool myMesh::triangulate(myFace *f)
{
	/**** TODO ****/
	return f->adjacent_halfedge->next->next->next == f->adjacent_halfedge ? true : false;
}

void myMesh::cleanOverlappingMesh(myVertex* v) {
	// iterate around face of v
	myHalfedge* step_h = v->originof;
	map<int, bool> visited_faces;
	map<string, myFace*> face_map;
	vector<myFace*> toDeteledFaces;
	// iterate around face of v
	do {
		myHalfedge* innerstep = step_h;
		if (step_h->adjacent_face == NULL)
		{
			//std::cout << "cleanOverlappingMesh >> boundary halfedge for " << step_h->id << endl;
			if (step_h->twin->next == step_h) break; // if the halfedge is alone, break
		}
		else {
			//std::cout << "cleanOverlappingMesh >> face " << step_h->adjacent_face->id << endl
			if (visited_faces.find(step_h->adjacent_face->id) != visited_faces.end()) {
				std::cout << "Error! Vertex loop indefinitely for " << v->id << std::endl;
				break;
			}
			visited_faces[step_h->adjacent_face->id] = true;

			// get the list of vertex of the face
			vector<int> vertex_list;
			do {
				vertex_list.push_back(innerstep->source->id);
				innerstep = innerstep->next;
			} while (innerstep->source->id != step_h->source->id);

			// convert list to string like "1,2,3,4"
			sort(vertex_list.begin(), vertex_list.end());
			stringstream ss;
			for (size_t i = 0; i < vertex_list.size(); i++)
			{
				ss << vertex_list[i];
				if (i != vertex_list.size() - 1) ss << ",";
			}
			string key = ss.str();
			//cout << "face key = " << key << endl;

			// if the face already exist, detele it, else add it 
			if (face_map.find(key) == face_map.end()) face_map[key] = step_h->adjacent_face;
			else toDeteledFaces.push_back(step_h->adjacent_face);
		}
		// go to next face
		step_h = step_h->twin->next;
	} while (step_h != v->originof);

	// delete faces
	for (myFace* f : toDeteledFaces) {
		std::cout << "delete face " << f->id << endl;
		myHalfedge* step = f->adjacent_halfedge;
		do {
			step->adjacent_face = NULL;
			step = step->next;
		} while (step != f->adjacent_halfedge);

		for (int i = static_cast<int>(faces.size()) - 1; i >= 0; i--)
		{
			if (faces[i] == f) {
				faces.erase(faces.begin() + i);
			}
		}
	}
}

void myMesh::collapse(myHalfedge* e) {
	if (e == NULL) return;
	//cout << "collapse edge " << e->id << endl;
	// calculate the middle point of the edge
	myVertex* v1 = e->source;
	myVertex* v2 = e->twin->source;
	//std::cout << "collapse v1 = " << v1->id << " & v2 = " << v2->id << endl;
	*v1->point += *v2->point;
	*v1->point /= 2; // v1 is the new_middle_vertex

	// change all halfedge of v2 to v1
	myHalfedge* step = v2->originof;
	do {
		step = step->twin->next;
		step->source = v1;
	} while (step != v2->originof);
	v1->originof = e->prev->twin;
	
	// collapse face
	myFace* toDelete_f1 = e->adjacent_face;
	myHalfedge* toDelete_half_1 = collapseFace(e);
	//cout << "toDelete_half_1 = " << toDelete_half_1 << endl;
	myHalfedge* toDelete_half_1_twin = (toDelete_half_1 != nullptr) ? toDelete_half_1->twin : nullptr;

	// collapse twin face
	myFace* toDelete_f2 = e->twin->adjacent_face;
	myHalfedge* toDelete_half_2 = collapseFace(e->twin);
	//cout << "toDelete_half_2 = " << toDelete_half_2 << endl;
	myHalfedge* toDelete_half_2_twin = (toDelete_half_2 != nullptr) ? toDelete_half_2->twin : nullptr;

	// delete vertice v2
	for (int i = static_cast<int>(vertices.size()) - 1; i >= 0; i--)
	{
		if (vertices[i] == v2) {
			//std::cout << "Delete vertice num " << vertices[i]->id << endl;
			vertices.erase(vertices.begin() + i);
		}
	}
	// delete faces
	for (int i = static_cast<int>(faces.size()) - 1; i >= 0; i--)
	{
		if ((faces[i] == toDelete_f1 && toDelete_half_1 != nullptr) || (faces[i] == toDelete_f2 && toDelete_half_2 != nullptr)) {
			//std::cout << "Delete face num " << faces[i]->id << endl;
			faces.erase(faces.begin() + i);
		}
	}
	// delete halfedge
	myHalfedge* e_twin = e->twin;
	for (int i = static_cast<int>(halfedges.size()) - 1; i >= 0; i--)
	{
		if (halfedges[i] == toDelete_half_1 || halfedges[i] == toDelete_half_1_twin || halfedges[i] == toDelete_half_2 || halfedges[i] == toDelete_half_2_twin || halfedges[i] == e || halfedges[i] == e_twin) {
			//std::cout << "Delete halfedge num " << halfedges[i]->id << endl;
			halfedges.erase(halfedges.begin() + i);
		}
	}
	cleanOverlappingMesh(v1);
	checkMesh();
}

myHalfedge* myMesh::collapseFace(myHalfedge* e) {
	if (e->adjacent_face == NULL)
	{
		//cout << "Boundary halfedge !" << endl;
		e->next->prev = e->prev;
		e->prev->next = e->next;
		e->next->source->originof = e->next;	// don't move this line, it's the result of 3h of debug 
		return nullptr;
	}

	myHalfedge* step = e;
	int cpt = 0;
	do
	{
		step = step->next;
		cpt++;
	} while (step != e);
	//std::cout << "mesh of " << cpt << " vertices for face :" << e->adjacent_face->id << endl;

	if (cpt == 3) {
		myHalfedge* halfedge_a = e->prev;
		myHalfedge* halfedge_b = e->next->twin;
		// change source orignof if it's b
		if (halfedge_b->source->originof == halfedge_b) {
			halfedge_b->source->originof = halfedge_a;
		}

		// replace halfedge b by a
		halfedge_a->next = halfedge_b->next;
		halfedge_a->prev = halfedge_b->prev;
		halfedge_a->adjacent_face = halfedge_b->adjacent_face;
		if (halfedge_a->adjacent_face != NULL) halfedge_a->adjacent_face->adjacent_halfedge = halfedge_a;

		halfedge_a->next->prev = halfedge_a;
		halfedge_a->prev->next = halfedge_a;

		return halfedge_b; // return halfedge to delete
	}
	else if (cpt > 3) {
		myHalfedge* next = e->next;
		myHalfedge* prev = e->prev;

		next->prev = prev;
		prev->next = next;

		e->adjacent_face->adjacent_halfedge = next;
		return nullptr;
	}
	else {
		std::cout << "Can't collapse face, not enought vertices !" << endl;
		return nullptr;
	}
}

bool isBoundary(myHalfedge* e) {
	return e->adjacent_face == NULL || e->twin->adjacent_face == NULL;
}

myHalfedge* myMesh::getShortestEdge() {
	map<int, int> twin_map;
	double min_dist = DBL_MAX;
	myHalfedge* min_edge = NULL;
	for (size_t i = 0; i < halfedges.size(); i++)
	{
		myHalfedge* h = halfedges[i];
		auto twin_it = twin_map.find(h->twin->id);
		if (twin_it == twin_map.end()) {
			twin_map[h->id] = h->twin->id;
		}
		else continue;

		// ====== check if collapse is possible ======
#pragma region Check if vertice is boundary
		bool h_vertex_is_boundary = false;
		bool twin_vertex_is_boundary = false;
		// check if h->source is a boundary vertex
		myHalfedge* steph = h;
		do {
			if (steph->adjacent_face == NULL) {
				h_vertex_is_boundary = true;
				break;
			}
			steph = steph->twin->next;
		} while (steph != h);

		// check if twin->source is a boundary vertex
		steph = h->twin;
		do {
			if (steph->adjacent_face == NULL) {
				twin_vertex_is_boundary = true;
				break;
			}
			steph = steph->twin->next;
		} while (steph != h->twin);
#pragma endregion

		if (h_vertex_is_boundary && twin_vertex_is_boundary && !isBoundary(h)) continue;

#pragma region Get vertex neighbor vertices
		vector<int> h_neighbor_vertex_id;
		vector<int> twin_neighbor_vertex_id;
		// get all neighbor vertex around the h source vertex
		myHalfedge* neighbor_v_steph = h->twin;
		do {
			h_neighbor_vertex_id.push_back(neighbor_v_steph->source->id);
			neighbor_v_steph = neighbor_v_steph->next->twin;
		} while (neighbor_v_steph != h->twin);

		// get all neighbor vertex around the twin h source vertex
		neighbor_v_steph = h->next->twin;
		do {
			twin_neighbor_vertex_id.push_back(neighbor_v_steph->source->id);
			neighbor_v_steph = neighbor_v_steph->next->twin;
		} while (neighbor_v_steph != h->next->twin);
#pragma endregion

		// intersection of the 2 vectors
		vector<int> intersection;
		sort(h_neighbor_vertex_id.begin(), h_neighbor_vertex_id.end());
		sort(twin_neighbor_vertex_id.begin(), twin_neighbor_vertex_id.end());
		set_intersection(
			h_neighbor_vertex_id.begin(),
			h_neighbor_vertex_id.end(),
			twin_neighbor_vertex_id.begin(),
			twin_neighbor_vertex_id.end(),
			back_inserter(intersection)
		);
		int nb_union = h_neighbor_vertex_id.size() + twin_neighbor_vertex_id.size() - intersection.size();
		//cout << "is Boundary" << isBoundary(h) << ", nb_union = " << nb_union << ", intersection.size() = " << intersection.size() << endl;
		// if it's a boundary halfedge and it's a triangle ==> no collapse
		if (isBoundary(h) && nb_union == 3 && intersection.size() == 1)
		{
			std::cout << "getShortestEdge >> ignore collapse solo triangle" << endl;
			continue;
		}
		myPoint3D* p1 = h->source->point;
		myPoint3D* p2 = h->twin->source->point;

		double d = p1->dist(*p2);
		if (d < min_dist) {
			min_dist = d;
			min_edge = h;
		}
	}
	return min_edge;
}