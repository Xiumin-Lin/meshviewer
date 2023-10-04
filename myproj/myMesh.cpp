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
	int face_id_cpt = 0;
	int halfedge_id_cpt = 0;
	int vertex_id_cpt = 1;

	ifstream fin(filename);
	if (!fin.is_open()) {
		cout << "Unable to open file!\n";
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
			cout << "v " << x << " " << y << " " << z << endl;

			myVertex* vertex = new myVertex();
			vertex->index = vertex_id_cpt++;	// debug with indice start at 1
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

			cout << "f"; 
			while (myline >> u)
			{
				int vertexIdx = atoi((u.substr(0, u.find("/"))).c_str());
				cout << " " << vertexIdx;
				vertexIdx = vertexIdx >= 0 ? vertexIdx - 1 : vertexIdx % verticesSize;
				list.push_back(vertexIdx);
			}
			cout << endl;

			myFace* face = new myFace();
			myHalfedge* originHalfedge = new myHalfedge();
			originHalfedge->index = halfedge_id_cpt++;		// debug with indice
			face->adjacent_halfedge = originHalfedge;
			
			int listSize = list.size();
			originHalfedge->source = this->vertices.at(list[0]);
			originHalfedge->source->originof = originHalfedge;
			originHalfedge->adjacent_face = face;

			halfedges.push_back(originHalfedge);

			myHalfedge* prevHalfedge = originHalfedge;
			for (size_t i = 1; i < listSize; i++)
			{
				myHalfedge* e = new myHalfedge();
				e->index = halfedge_id_cpt++;	// debug with indice
				e->source = this->vertices.at(list[i]);
				e->source->originof = e;
				e->adjacent_face = face;
				e->prev = prevHalfedge;

				prevHalfedge->next = e;

				it = twin_map.find(make_pair(list[(i + 1) % listSize], list[i]));
				if (it == twin_map.end()) { 
					twin_map[make_pair(list[i], list[(i + 1) % listSize])] = e;
				}
				else { 
					e->twin = it->second;
					e->twin->twin = e;
				}

				prevHalfedge = e;
				if (i == listSize - 1) {
					e->next = originHalfedge;
					originHalfedge->prev = e;

					it = twin_map.find(make_pair(list[1], list[0]));
					if (it == twin_map.end()) {
						twin_map[make_pair(list[0], list[1])] = originHalfedge;
					}
					else {
						originHalfedge->twin = it->second;
						originHalfedge->twin->twin = originHalfedge;
					}
				}

				halfedges.push_back(e);
			}
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
	/*for (size_t i = 0; i < faces.size(); i++)
	{
		faces[i]->computeNormal();
	}*/

	/*for (size_t i = 0; i < vertices.size(); i++)
	{
		vertices[i]->computeNormal();
	}*/
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


void myMesh::triangulate()
{
	/**** TODO ****/
}

//return false if already triangle, true othewise.
bool myMesh::triangulate(myFace *f)
{
	/**** TODO ****/
	return false;
}

