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
	for (unsigned int i = 0; i < vertices.size(); i++) if (vertices[i]) delete vertices[i];
	for (unsigned int i = 0; i < halfedges.size(); i++) if (halfedges[i]) delete halfedges[i];
	for (unsigned int i = 0; i < faces.size(); i++) if (faces[i]) delete faces[i];

	vector<myVertex *> empty_vertices;    vertices.swap(empty_vertices);
	vector<myHalfedge *> empty_halfedges; halfedges.swap(empty_halfedges);
	vector<myFace *> empty_faces;         faces.swap(empty_faces);
}

void myMesh::checkMesh()
{
	bool allHasTwin = true;
	bool allHasNext = true;
	bool allHasPrev = true;
	for (myHalfedge* he : halfedges)
	{
		if (allHasNext) {
			if (he->next == NULL || he->next == nullptr)
			{
				std::cout << "Error! Not all halfedge have their next!\n";
				allHasNext = false;
			}
			else if (he->next->prev != he) {
				std::cout << "Error! halfedge->next->prev should be the halfedge itself\n";
				allHasNext = false;
			}
		}
		if (allHasPrev) {
			if (he->prev == NULL || he->prev == nullptr)
			{
				std::cout << "Error! Not all halfedge have their prev!\n";
				allHasPrev = false;
			}
			else if (he->prev->next != he) {
				std::cout << "Error! halfedge->prev->next should be the halfedge itself!\n";
				allHasPrev = false;
			}
		}
		if (allHasTwin) {
			if (he->twin == NULL || he->twin == nullptr)
			{
				std::cout << "Error! Not all halfedge have their twins!\n";
				allHasTwin = false;
			}
			else if (he->twin->twin != he) {
				std::cout << "Error! Twin of halfedge is not the halfedge itself!\n";
				allHasTwin = false;
			}
			else if (he->next == NULL || he->next == nullptr || he->twin->source != he->next->source) {
				std::cout << "Error! Twin of halfedge has not the same source as the next of halfedge!\n";
				allHasTwin = false;
			}
		}
	}
	if (allHasNext)	std::cout << "Each halfedge has a next!\n";
	if (allHasPrev) std::cout << "Each halfedge has a prev!\n";
	if (allHasTwin) std::cout << "Each halfedge has a twin!\n";
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
			for (size_t i = 0; i < nbVertex; i++)
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


void myMesh::splitFaceTRIS(myFace *f, myPoint3D *p)
{
	/**** TODO ****/
}

void myMesh::splitEdge(myHalfedge *e1, myPoint3D *p)
{

	/**** TODO ****/
}

void myMesh::splitFaceQUADS(myFace *f, myPoint3D *p)
{
	/**** TODO ****/
}

myVertex* computeFaceCenterPoint(myFace* face) {
	myPoint3D* center_p = new myPoint3D();

	double vertex_cpt = 0;
	myHalfedge* step_halfedge = face->adjacent_halfedge;
	do
	{
		*center_p += *(step_halfedge->source->point);
		vertex_cpt++;
		step_halfedge = step_halfedge->next;
	} while (face->adjacent_halfedge != step_halfedge);
	*center_p / vertex_cpt;

	myVertex* center_vertex = new myVertex();
	center_vertex->point = center_p;
	return center_vertex;
}

myVertex* computeEdgeMidPoint(myHalfedge* e) {
	myVertex* middle_e = new myVertex();
	*(middle_e->point) = (*(e->source->point) + *(e->twin->source->point)) / 2;
	return middle_e;
}

// https://www.youtube.com/watch?v=mfp1Z1mBClc
void myMesh::subdivisionCatmullClark()
{
	/**** TODO ****/
	//vector<myVertex*> newFacePoints;
	//vector<myVertex*> newEdgepoints;
	//vector<myFace*> newFaces;

	//// Step 1: Create new vertices as face center point
	//for (auto& face : faces) {
	//	newFacePoints.push_back(computeFaceCenterPoint(face));
	//}

	//// Step 2: Create new vertices as edge middle point
	//for (auto& edge : halfedges) {
	//	newEdgepoints.push_back(computeEdgeMidPoint(edge));
	//}

	// Step 3: Connection
	// addNewVertexToEdge
		// create 2 new halfedge
		// create 2 new twin
		// link new halfedge to edge vertices
		// delete old halfedge and twin

	// addVertexToFaceCenter
		// link center_v to n middle edge point of the face
		// create n new faces
		// create 2 new halfedges for each new face
		// link halfedges
		// link twin if possible
		// delete old face

	// Step 4: Move newEdgepoints
		// e = (originP_1 + originP_2 + newFaceCenterP_1 + newFaceCenterP_2) / 4

	// Step 5: Move origins point
		// new_OriginP = Q / n + 2R / n + (n - 3)S/n => n is the number of faces the originP touch
		// Q = sum(all newFaceCenterP of face that touch the originP) / n
		// R = sum(all newMiddleEdgeP of edge that touch the originP) / n
		// S = originP

	map<myFace*, myVertex*> facePointsMap;
	map<myHalfedge*, myVertex*> halfedgesPointsMap;

	vector<myVertex*> newVertices;
	vector<myFace*> newFaces;

	for (size_t idx = 0; idx < faces.size(); idx++) {
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

	for (size_t idx = 0; idx < halfedges.size(); idx++)
	{
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

	for (size_t idx = 0; idx < vertices.size(); idx++)
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
		*newPoint += *avgFacePoint + *avgEdgePoint * 2 + *vertices[idx]->point * (cpt - 3);

		vertices[idx]->point = newPoint;
	}
}

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

			myHalfedge* originHalfedge = face->adjacent_halfedge;

			do
			{
				// create face
				myFace* new_face = new myFace();
				new_face->adjacent_halfedge = originHalfedge;
				originHalfedge->adjacent_face = new_face;

				// create halfedge
				myHalfedge* originHalfedgeNext = originHalfedge->next;

				myHalfedge* nextHalf = new myHalfedge();
				nextHalf->source = originHalfedgeNext->source;
				nextHalf->source->originof = nextHalf;
				nextHalf->adjacent_face = new_face;

				myHalfedge* prevHalf = new myHalfedge();
				prevHalf->source = center_vertex;
				prevHalf->source->originof = prevHalf;
				prevHalf->adjacent_face = new_face;

				// associate halfedge
				nextHalf->prev = originHalfedge;
				nextHalf->next = prevHalf;

				prevHalf->prev = nextHalf;
				prevHalf->next = originHalfedge;

				originHalfedge->prev = prevHalf;
				originHalfedge->next = nextHalf;

				// create twin
				associateTwin(twin_map, nextHalf->source->id, prevHalf->source->id, nextHalf);
				associateTwin(twin_map, prevHalf->source->id, originHalfedge->source->id, prevHalf);

				// save created face and halfedge
				faces.push_back(new_face);

				halfedges.push_back(nextHalf);
				halfedges.push_back(prevHalf);

				// pass to next globale halfedge
				originHalfedge = originHalfedgeNext;
			} while (face->adjacent_halfedge != originHalfedge);

			// delete old face
			delete this->faces[idx];
			this->faces[idx] = nullptr;
		}
		else
		{
			std::cout << "Error when tryng triangulate face id " << face->id << " : have only " << vertex_cpt << " halfedges !" << endl;
		}
	}

	for (int i = static_cast<int>(faces.size()) - 1; i >= 0; i--)
	{
		if (faces[i] == nullptr) {
			faces.erase(faces.begin() + i);
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

void myMesh::collapse(myHalfedge* e) {
	// calculate the center point
	myVertex* v1 = e->source;
	myVertex* v2 = e->twin->source;
	//std::cout << "collapse v1 = " << v1->id << " & v2 = " << v2->id << endl;
	*v1->point += *v2->point;
	*v1->point /= 2; // v1 is the new_center_vertex

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
	myHalfedge* toDelete_half_1_twin = (toDelete_half_1 != nullptr) ? toDelete_half_1->twin : nullptr;

	// collapse twin face
	myFace* toDelete_f2 = e->twin->adjacent_face;
	myHalfedge* toDelete_half_2 = collapseFace(e->twin);
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
	checkMesh();
}

myHalfedge* myMesh::collapseFace(myHalfedge* e) {
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
			halfedge_b->source->originof = halfedge_b->twin->next;
		}

		// replace halfedge b by a
		halfedge_a->next = halfedge_b->next;
		halfedge_a->prev = halfedge_b->prev;
		halfedge_a->adjacent_face = halfedge_b->adjacent_face;
		halfedge_a->adjacent_face->adjacent_halfedge = halfedge_a;

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

myHalfedge* myMesh::getShortestEdge() {
	map<int, int> twin_map;
	vector<pair<double, myHalfedge*>> distances;
	for (size_t i = 0; i < halfedges.size(); i++)
	{
		myHalfedge* h = halfedges[i];
		auto twin_it = twin_map.find(h->twin->id);
		if (twin_it == twin_map.end()) {
			twin_map[h->id] = h->twin->id;
		}
		else continue;

		myPoint3D* p1 = h->source->point;
		myPoint3D* p2 = h->twin->source->point;

		distances.push_back(make_pair(p1->dist(*p2), h));
	}

	sort(distances.begin(), distances.end()); // TODO to optimize
	return distances[0].second;
}