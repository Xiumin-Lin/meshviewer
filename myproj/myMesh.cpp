#include "myMesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <GL/glew.h>
#include "myvector3d.h"

using namespace std;

myMesh::myMesh(void)
{
	/**** TODO ****/
}


myMesh::~myMesh(void)
{
	/**** TODO ****/
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
	vector<myHalfedge *>::iterator it;
	for (it = halfedges.begin(); it != halfedges.end(); it++)
	{
		if ((*it)->twin == NULL)
			break;
	}
	if (it != halfedges.end())
		cout << "Error! Not all edges have their twins!\n";
	else cout << "Each edge has a twin!\n";
}


bool myMesh::readFile(std::string filename)
{
	string s, t, u;
	//vector<int> faceids;
	//myHalfedge **hedges;

	ifstream fin(filename);
	if (!fin.is_open()) {
		cout << "Unable to open file!\n";
		return false;
	}
	name = filename;

	map<pair<int, int>, myHalfedge *> twin_map;
	map<pair<int, int>, myHalfedge *>::iterator it;

	int face_id_cpt = 0;
	int halfedge_id_cpt = 0;

	while (getline(fin, s))
	{
		stringstream myline(s);
		myline >> t;
		if (t == "g") {}
		else if (t == "v")
		{
			float x, y, z;
			myline >> x >> y >> z;
			//cout << "v " << x << " " << y << " " << z << endl;

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
			int verticesSize = vertices.size();

			//cout << "f"; 
			while (myline >> u)
			{
				int vertexIdx = atoi((u.substr(0, u.find("/"))).c_str());
				//cout << " " << vertexIdx;
				vertexIdx = vertexIdx >= 0 ? vertexIdx - 1 : vertexIdx % verticesSize;
				list.push_back(vertexIdx);
			}
			//cout << endl;

			int listSize = list.size();
			myFace* face = new myFace();
			myHalfedge* originHalfedge = new myHalfedge();

			myHalfedge* prevHalfedge = nullptr;
			for (size_t i = 0; i < listSize; i++)
			{
				myHalfedge* e = new myHalfedge();
				e->index = halfedge_id_cpt++;	// debug with indice
				e->source = this->vertices.at(list[i]);
				e->source->originof = e;
				e->adjacent_face = face;
				e->prev = prevHalfedge;

				if (i == 0) 
				{
					originHalfedge = e;
					prevHalfedge = e;
				}
				else prevHalfedge->next = e;

				int vertex_a = list[i];
				int vertex_b = list[(i + 1) % listSize];
				it = twin_map.find(make_pair(vertex_b, vertex_a));
				if (it == twin_map.end()) { 
					twin_map[make_pair(vertex_a, vertex_b)] = e;
				}
				else { 
					e->twin = it->second;
					e->twin->twin = e;
				}

				if (i == listSize - 1) {
					e->next = originHalfedge;
					originHalfedge->prev = e;
				}

				halfedges.push_back(e);
				prevHalfedge = e;
			}

			face->adjacent_halfedge = originHalfedge;
			face->index = face_id_cpt++;
			faces.push_back(face);
		}
	}

	checkMesh();
	normalize();

	return true;
}


void myMesh::computeNormals()
{
	/**** TODO ****/
	for (size_t i = 0; i < faces.size(); i++)
	{
		faces[i]->computeNormal();
	}

	for (size_t i = 0; i < vertices.size(); i++)
	{
		vertices[i]->computeNormal();
	}
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


void myMesh::subdivisionCatmullClark()
{
	/**** TODO ****/
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
				associateTwin(twin_map, nextHalf->source->index, prevHalf->source->index, nextHalf);
				associateTwin(twin_map, prevHalf->source->index, originHalfedge->source->index, prevHalf);

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
			cout << "Error when tryng triangulate face id " << face->index << " : have only " << vertex_cpt << " halfedges !" << endl;
		}
	}

	for (vector<myFace*>::iterator it = faces.begin(); it != faces.end(); ) {
		if (*it == nullptr) {
			faces.erase(it);
		}
		else {
			it++;
		}
	}
}

//return false if already triangle, true othewise.
bool myMesh::triangulate(myFace *f)
{
	/**** TODO ****/
	return f->adjacent_halfedge->next->next->next == f->adjacent_halfedge ? true : false;
}
