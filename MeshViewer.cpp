#include<iostream>
#include<map>
#include<math.h>
#include<fstream>
#include<string>
#include<vector>
#include<sstream>
#include<windows.h>
#include<stdio.h>
#include<glut.h>

using namespace std;

struct HE_vert;
struct HE_face;
struct HE_edge;

int  press_x, press_y;
float x_angle = 0.0, y_angle = 0.0, r, scal;
static int obj_mode = 0, mouse_mode = 0;
static float scale_size = 1;
float xcenter, ycenter, zcenter, xmin, xmax, ymin, ymax, zmin, zmax;
char modelname;

#define OBJ_POINTS	  1
#define OBJ_WIREFRAME 2
#define OBJ_FLAT	  3 
#define OBJ_SMOOTH	  4 

#define TRANSFORM_NONE    0 
#define TRANSFORM_ROTATE  1
#define TRANSFORM_SCALE   2 
#define TRANSFORM_TRANSLATE 3

struct HE_edge {
	HE_vert* vert;
	HE_edge* pair;
	HE_face* face;
	HE_edge* prev;
	HE_edge* next;
};

struct HE_vert {
	int vid;
	double x, y, z;
	double vnx = 0;
	double vny = 0;
	double vnz = 0;
	HE_edge* edge = NULL;
	double FArea;
	std::vector<double> Areaofadjacentindividualface;
	std::vector<int> adjacentfaceid;
};

struct HE_face {
	int fid;
	int v1, v2, v3;
	HE_edge* edge;
	double fnx, fny, fnz;
	double Area;
};

vector<HE_vert*> vVertex;
vector<HE_face*> vFaces;
vector<HE_edge*> vEdges;

void Initialize()
{
	glEnable(GL_DEPTH_TEST);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, 1, .1, 100);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(2, 2, 2, 0, 0, 0, 0, 0, 1);

	//lighting
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_COLOR_MATERIAL);
	GLfloat light_position[] = { 1.8, 1.8, 1.8, 0.0 };
	GLfloat Alight[] = { 0.2, 0.2, 0.2, 0.2 };
	GLfloat Dlight[] = { 1.0, 1.0, 1.0, 1.0 };
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_AMBIENT, Alight);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, Dlight);

	glRotatef(x_angle, 0, 1, 0);
	glRotatef(y_angle, 1, 0, 0);

	glScalef(scale_size, scale_size, scale_size);
}

void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, 1, .1, 100);
	glMatrixMode(GL_MODELVIEW);
}

void mouseButton(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN)
	{
		press_x = x; press_y = y;
		if (button == GLUT_LEFT_BUTTON)
			mouse_mode = TRANSFORM_ROTATE;
		else if (button == GLUT_MIDDLE_BUTTON)
			mouse_mode = TRANSFORM_SCALE;
		else if (button == GLUT_RIGHT_BUTTON)
			mouse_mode = TRANSFORM_TRANSLATE;

	}
	else if (state == GLUT_UP)
	{
		mouse_mode = TRANSFORM_NONE;
	}
}

void mouseMove(int x, int y)
{
	if (mouse_mode == TRANSFORM_ROTATE)
	{
		x_angle += (x - press_x) / 5.0;

		if (x_angle > 180)
			x_angle -= 360;
		else if (x_angle < -180)
			x_angle += 360;

		press_x = x;

		y_angle += (y - press_y) / 5.0;

		if (y_angle > 180)
			y_angle -= 360;
		else if (y_angle < -180)
			y_angle += 360;

		press_y = y;
	}
	else if (mouse_mode == TRANSFORM_SCALE)
	{
		float old_size = scale_size;

		scale_size *= (1 + (y - press_y) / 60.0);

		if (scale_size < 0)
			scale_size = old_size;
		press_y = y;
	}

	glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'P':
	case 'p':
		obj_mode = OBJ_POINTS;
		break;

	case 'W':
	case 'w':
		obj_mode = OBJ_WIREFRAME;
		break;

	case 'F':
	case 'f':
		obj_mode = OBJ_FLAT;
		break;

	case 'S':
	case 's':
		obj_mode = OBJ_SMOOTH;
		break;
	}
	glutPostRedisplay();
}

void Gridlines() {
	//Drawing grid lines

	int xmn = -100, xmx = 100;
	int ymn = -100, ymx = 100;
	glColor3f(1, 1, 1);
	glBegin(GL_LINES);
	for (int x = xmn; x <= 100; x++)
	{
		glVertex3f(x, ymn, 0.0);
		glVertex3f(x, ymx, 0.0);
	}
	for (int y = ymn; y <= 100; y++)
	{
		glVertex3f(xmn, y, 0.0);
		glVertex3f(xmx, y, 0.0);
	}
	glEnd();
}

void DrawXAxis() {
	glPushMatrix();
	glColor3f(0, 0, 1);
	glRotated(90, 0, 1, 0);
	GLUquadricObj* xobj = gluNewQuadric();
	gluCylinder(xobj, 0.1, 0.1, 1.5, 3, 4);
	glTranslated(0, 0, 1.5);
	glutSolidCone(0.2, 0.5, 3, 4);
	glPopMatrix();
}

void DrawYAxis() {
	glPushMatrix();
	glColor3f(0, 1, 0);
	glRotated(-90, 1, 0, 0);
	GLUquadricObj* yobj = gluNewQuadric();
	gluCylinder(yobj, 0.1, 0.1, 1.5, 3, 4);
	glTranslated(0, 0, 1.5);
	glutSolidCone(0.2, 0.5, 3, 4);
	glPopMatrix();
}

void DrawZAxis() {
	glPushMatrix();
	glColor3f(1, 0, 0);
	GLUquadricObj* zobj = gluNewQuadric();
	gluCylinder(zobj, 0.1, 0.1, 1.5, 3, 4);
	glTranslated(0, 0, 1.5);
	glutSolidCone(0.2, 0.5, 3, 4);
	glPopMatrix();
}

void Axes() {
	DrawXAxis();
	DrawYAxis();
	DrawZAxis();
}

void Boundingbox() {

	xmin = vVertex[0]->x, ymin = vVertex[0]->y, zmin = vVertex[0]->z;
	xmax = vVertex[0]->x, ymax = vVertex[0]->y, zmax = vVertex[0]->z;

	for (int i = 1; i != vVertex.size(); i++)
	{
		if (vVertex[i]->x < xmin) xmin = vVertex[i]->x;
		if (vVertex[i]->y < ymin) ymin = vVertex[i]->y;
		if (vVertex[i]->z < zmin) zmin = vVertex[i]->z;
		if (vVertex[i]->x > xmax) xmax = vVertex[i]->x;
		if (vVertex[i]->y > ymax) ymax = vVertex[i]->y;
		if (vVertex[i]->z > zmax) zmax = vVertex[i]->z;
	}

	xcenter = ((xmin + xmax) / 2);
	ycenter = ((ymin + ymax) / 2);
	zcenter = ((zmin + zmax) / 2);

	r = sqrt((xmax - xcenter) * (xmax - xcenter) + (ymax - ycenter) * (ymax - ycenter) + (zmax - zcenter) * (zmax - zcenter));
	scal = 1 / r;


	glColor3f(0.99, 0.50, 1.00);
	glScalef(scal, scal, scal);
	glTranslatef(-xcenter, -ycenter, -zcenter);

	glBegin(GL_LINES);
	{
		glVertex3f(xmin, ymin, zmin); glVertex3f(xmax, ymin, zmin);
		glVertex3f(xmax, ymin, zmin); glVertex3f(xmax, ymax, zmin);
		glVertex3f(xmax, ymax, zmin); glVertex3f(xmin, ymax, zmin);
		glVertex3f(xmin, ymax, zmin); glVertex3f(xmin, ymin, zmin);
		glVertex3f(xmin, ymin, zmax); glVertex3f(xmax, ymin, zmax);
		glVertex3f(xmax, ymin, zmax); glVertex3f(xmax, ymax, zmax);
		glVertex3f(xmax, ymax, zmax); glVertex3f(xmin, ymax, zmax);
		glVertex3f(xmin, ymax, zmax); glVertex3f(xmin, ymin, zmax);
		glVertex3f(xmin, ymin, zmin); glVertex3f(xmin, ymin, zmax);
		glVertex3f(xmax, ymin, zmin); glVertex3f(xmax, ymin, zmax);
		glVertex3f(xmax, ymax, zmin); glVertex3f(xmax, ymax, zmax);
		glVertex3f(xmin, ymax, zmin); glVertex3f(xmin, ymax, zmax);
	}
	glEnd();
}

void DrawModelPoints() {
	for (int a = 0; a != vVertex.size(); a++) {
		glBegin(GL_POINTS);
		glVertex3f(vVertex[a]->x, vVertex[a]->y, vVertex[a]->z);
		glEnd();
	}
}

void DrawModelWireframe() {
	for (int a = 0; a != vFaces.size(); a++) {
		glBegin(GL_LINE_LOOP);
		glVertex3f(vVertex[(vFaces[a]->v1) - 1]->x, vVertex[(vFaces[a]->v1) - 1]->y, vVertex[(vFaces[a]->v1) - 1]->z);
		glVertex3f(vVertex[(vFaces[a]->v2) - 1]->x, vVertex[(vFaces[a]->v2) - 1]->y, vVertex[(vFaces[a]->v2) - 1]->z);
		glVertex3f(vVertex[(vFaces[a]->v3) - 1]->x, vVertex[(vFaces[a]->v3) - 1]->y, vVertex[(vFaces[a]->v3) - 1]->z);
		glEnd();
	}
}

void DrawFlatShadingModel() {
	glShadeModel(GL_FLAT);
	for (int a = 0; a != vFaces.size(); a++) {
		glBegin(GL_TRIANGLES);
		glNormal3f(vVertex[(vFaces[a]->v1) - 1]->vnx, vVertex[(vFaces[a]->v1) - 1]->vny, vVertex[(vFaces[a]->v1) - 1]->vnz);
		glVertex3f(vVertex[(vFaces[a]->v1) - 1]->x, vVertex[(vFaces[a]->v1) - 1]->y, vVertex[(vFaces[a]->v1) - 1]->z);
		glNormal3f(vVertex[(vFaces[a]->v2) - 1]->vnx, vVertex[(vFaces[a]->v2) - 1]->vny, vVertex[(vFaces[a]->v2) - 1]->vnz);
		glVertex3f(vVertex[(vFaces[a]->v2) - 1]->x, vVertex[(vFaces[a]->v2) - 1]->y, vVertex[(vFaces[a]->v2) - 1]->z);
		glNormal3f(vVertex[(vFaces[a]->v3) - 1]->vnx, vVertex[(vFaces[a]->v3) - 1]->vny, vVertex[(vFaces[a]->v3) - 1]->vnz);
		glVertex3f(vVertex[(vFaces[a]->v3) - 1]->x, vVertex[(vFaces[a]->v3) - 1]->y, vVertex[(vFaces[a]->v3) - 1]->z);
		glEnd();
	}
}

void DrawShadeModel() {
	glShadeModel(GL_SMOOTH);
	for (int a = 0; a != vFaces.size(); a++) {
		glBegin(GL_TRIANGLES);
		glNormal3f(vVertex[(vFaces[a]->v1) - 1]->vnx, vVertex[(vFaces[a]->v1) - 1]->vny, vVertex[(vFaces[a]->v1) - 1]->vnz);
		glVertex3f(vVertex[(vFaces[a]->v1) - 1]->x, vVertex[(vFaces[a]->v1) - 1]->y, vVertex[(vFaces[a]->v1) - 1]->z);
		glNormal3f(vVertex[(vFaces[a]->v2) - 1]->vnx, vVertex[(vFaces[a]->v2) - 1]->vny, vVertex[(vFaces[a]->v2) - 1]->vnz);
		glVertex3f(vVertex[(vFaces[a]->v2) - 1]->x, vVertex[(vFaces[a]->v2) - 1]->y, vVertex[(vFaces[a]->v2) - 1]->z);
		glNormal3f(vVertex[(vFaces[a]->v3) - 1]->vnx, vVertex[(vFaces[a]->v3) - 1]->vny, vVertex[(vFaces[a]->v3) - 1]->vnz);
		glVertex3f(vVertex[(vFaces[a]->v3) - 1]->x, vVertex[(vFaces[a]->v3) - 1]->y, vVertex[(vFaces[a]->v3) - 1]->z);
		glEnd();
	}
}

void Draw()
{
	Initialize();
	Gridlines();
	Axes();
	Boundingbox();

	if (obj_mode == OBJ_POINTS) {
		DrawModelPoints();
	}
	else if (obj_mode == OBJ_WIREFRAME) {
		DrawModelWireframe();
	}
	else if (obj_mode == OBJ_FLAT) {
		DrawFlatShadingModel();
	}
	else if (obj_mode == OBJ_SMOOTH) {
		DrawShadeModel();
	}

	glutSwapBuffers();
}

vector<string> ReadSelectedModel(char modelname) {

	string fileLine;
	vector<string> fileLines;

	ifstream infile;

	switch (modelname) {
	case 'B':case 'b':
		infile.open("bunny.m");
		if (!infile.is_open()) {
			cerr << "Error opening file b.m" << endl;
		}
		while (getline(infile, fileLine))
		{
			fileLines.push_back(fileLine);
		}
		infile.close();
		break;
	case 'C':case 'c':
		infile.open("cap.m");
		while (!infile.eof())
		{
			getline(infile, fileLine);
			fileLines.push_back(fileLine);
		}

		infile.close();
		break;
	case 'E':case 'e':
		infile.open("eight.m");
		while (!infile.eof())
		{
			getline(infile, fileLine);
			fileLines.push_back(fileLine);
		}

		infile.close();
		break;
	case 'G':case 'g':
		infile.open("gargoyle.m");
		while (!infile.eof())
		{
			getline(infile, fileLine);
			fileLines.push_back(fileLine);
		}

		infile.close();
		break;
	case 'K':case 'k':infile.open("knot.m");
		while (!infile.eof())
		{
			getline(infile, fileLine);
			fileLines.push_back(fileLine);
		}

		infile.close();
		break;
	}

	return fileLines;
}

void BuildVerticesandFaces(vector<string> fileLines) {

	vector<std::string>::iterator tokensiterator;

	for (int i = 0; i != fileLines.size(); i++)
	{
		string word;
		vector<string> tokens;
		stringstream ss(fileLines[i]);
		while (ss >> word) {
			tokens.push_back(word);
		}

		for (tokensiterator = tokens.begin(); tokensiterator != tokens.end(); tokensiterator++)
		{
			if ((*tokensiterator).compare("Vertex") == 0)
			{
				HE_vert* V = new HE_vert();

				tokensiterator++;
				V->vid = stoi((*tokensiterator));
				tokensiterator++;
				V->x = stod((*tokensiterator));
				tokensiterator++;
				V->y = stod((*tokensiterator));
				tokensiterator++;
				V->z = stod((*tokensiterator));
				vVertex.push_back(V);
			}

			if ((*tokensiterator).compare("Face") == 0)
			{
				HE_face* F = new HE_face();

				tokensiterator++;
				F->fid = stoi((*tokensiterator));
				tokensiterator++;
				F->v1 = stoi((*tokensiterator));
				tokensiterator++;
				F->v2 = stoi((*tokensiterator));
				tokensiterator++;
				F->v3 = stoi((*tokensiterator));
				vFaces.push_back(F);
			}
		}
	}
}

char Welcome() {
	cout << "Choose which model M file need to be opened" << endl;
	cout << "\n" << endl;
	cout << "B or b for bunny model" << endl;
	cout << "C or c for cap model" << endl;
	cout << "E or e for eight model" << endl;
	cout << "G or g for gargoyle model" << endl;
	cout << "K or k for knot model" << endl;
	cout << "\n" << endl;
	cout << "Please use below alphabet keys for rendering model in different form" << endl;
	cout << "\n" << endl;
	cout << "P or p for Points" << endl;
	cout << "W or w for wireframe" << endl;
	cout << "F or f for FlatShading" << endl;
	cout << "S or s for SmoothShading" << endl;
	cin >> modelname;

	return modelname;
}

void  BuildHalfEdgeDataStructure(map<pair<int, int>, HE_edge*>& edges,
	pair<int, int> edgekey,
	pair<int, int> nextEdgeKey,
	pair<int, int> prevEdgeKey,
	int emantingVertexIndex,
	int endingVertexIndex,
	int faceIndex) {

	edges[edgekey]->face = vFaces[faceIndex];
	edges[edgekey]->vert = vVertex[(endingVertexIndex)-1];
	edges[edgekey]->next = edges[nextEdgeKey];
	edges[edgekey]->prev = edges[prevEdgeKey];

	vFaces[faceIndex]->edge = edges[prevEdgeKey];

	if (vVertex[(emantingVertexIndex)-1]->edge == NULL)
	{
		vVertex[(emantingVertexIndex)-1]->edge = edges[edgekey];
	}

	vEdges.push_back(edges[edgekey]);
}

map<pair<int, int>, HE_edge*> BuildFacesHalfEdgeDataStructure() {

	map<pair<int, int>, HE_edge*> edges;

	for (int i = 0; i != vFaces.size(); i++)
	{
		int faceVertex1Index = vFaces[i]->v1;
		int faceVertex2Index = vFaces[i]->v2;
		int faceVertex3Index = vFaces[i]->v3;

		pair<int, int> edgekey = make_pair(faceVertex1Index, faceVertex2Index);
		pair<int, int> nextEdgeKey = make_pair(faceVertex2Index, faceVertex3Index);
		pair<int, int> prevEdgeKey = make_pair(faceVertex3Index, faceVertex1Index);
		int emantingVertexIndex = faceVertex1Index;
		int endingVertexIndex = faceVertex2Index;

		edges[edgekey] = new HE_edge();
		edges[nextEdgeKey] = new HE_edge();
		edges[prevEdgeKey] = new HE_edge();

		BuildHalfEdgeDataStructure(edges, edgekey, nextEdgeKey, prevEdgeKey, emantingVertexIndex, endingVertexIndex, i);

		edgekey = make_pair(faceVertex2Index, faceVertex3Index);
		nextEdgeKey = make_pair(faceVertex3Index, faceVertex1Index);
		prevEdgeKey = make_pair(faceVertex1Index, faceVertex2Index);
		emantingVertexIndex = faceVertex2Index;
		endingVertexIndex = faceVertex3Index;

		BuildHalfEdgeDataStructure(edges, edgekey, nextEdgeKey, prevEdgeKey, emantingVertexIndex, endingVertexIndex, i);

		edgekey = make_pair(faceVertex3Index, faceVertex1Index);
		nextEdgeKey = make_pair(faceVertex1Index, faceVertex2Index);
		prevEdgeKey = make_pair(faceVertex2Index, faceVertex3Index);
		emantingVertexIndex = faceVertex3Index;
		endingVertexIndex = faceVertex1Index;

		BuildHalfEdgeDataStructure(edges, edgekey, nextEdgeKey, prevEdgeKey, emantingVertexIndex, endingVertexIndex, i);
	}

	return edges;
}

void BuildTwinEdgeDataStructure(map<pair<int, int>, HE_edge*>& edges) {
	for (int i = 0; i != vFaces.size(); i++)
	{
		int faceVertex1Index = vFaces[i]->v1;
		int faceVertex2Index = vFaces[i]->v2;
		int faceVertex3Index = vFaces[i]->v3;

		pair<int, int> edgekey = make_pair(faceVertex1Index, faceVertex2Index);
		pair<int, int> twinEdgekey = make_pair(faceVertex2Index, faceVertex1Index);

		if (edges.count(edgekey) > 0 && edges.count(twinEdgekey) > 0)
		{
			edges[edgekey]->pair = edges[twinEdgekey];
			edges[twinEdgekey]->pair = edges[edgekey];
		}

		edgekey = make_pair(faceVertex2Index, faceVertex3Index);
		twinEdgekey = make_pair(faceVertex3Index, faceVertex2Index);

		if (edges.count(edgekey) > 0 && edges.count(twinEdgekey) > 0)
		{
			edges[edgekey]->pair = edges[twinEdgekey];
			edges[twinEdgekey]->pair = edges[edgekey];
		}

		edgekey = make_pair(faceVertex3Index, faceVertex1Index);
		twinEdgekey = make_pair(faceVertex1Index, faceVertex3Index);

		if (edges.count(edgekey) > 0 && edges.count(twinEdgekey) > 0)
		{
			edges[edgekey]->pair = edges[twinEdgekey];
			edges[twinEdgekey]->pair = edges[edgekey];
		}
	}
}

void CalculateFaceNormal() {
	for (int i = 0; i != vFaces.size(); i++)
	{
		double edge[3];
		double adjacentEdge[3];
		double Nx, Ny, Nz, length;

		int faceVertex1Index = vFaces[i]->v1;
		int faceVertex2Index = vFaces[i]->v2;
		int faceVertex3Index = vFaces[i]->v3;

		edge[0] = ((vVertex[(faceVertex2Index)-1]->x) - (vVertex[(faceVertex1Index)-1]->x));
		edge[1] = ((vVertex[(faceVertex2Index)-1]->y) - (vVertex[(faceVertex1Index)-1]->y));
		edge[2] = ((vVertex[(faceVertex2Index)-1]->z) - (vVertex[(faceVertex1Index)-1]->z));

		adjacentEdge[0] = ((vVertex[(faceVertex3Index)-1]->x) - (vVertex[(faceVertex1Index)-1]->x));
		adjacentEdge[1] = ((vVertex[(faceVertex3Index)-1]->y) - (vVertex[(faceVertex1Index)-1]->y));
		adjacentEdge[2] = ((vVertex[(faceVertex3Index)-1]->z) - (vVertex[(faceVertex1Index)-1]->z));

		Nx = (((edge[1]) * (adjacentEdge[2])) - ((edge[2]) * (adjacentEdge[1])));
		Ny = (((edge[2]) * (adjacentEdge[0])) - ((edge[0]) * (adjacentEdge[2])));
		Nz = (((edge[0]) * (adjacentEdge[1])) - ((edge[1]) * (adjacentEdge[0])));

		length = sqrt(((Nx) * (Nx)) + ((Ny) * (Ny)) + ((Nz) * (Nz)));

		vFaces[i]->Area = ((length) / 2);

		vFaces[i]->fnx = ((Nx) / (length));
		vFaces[i]->fny = ((Ny) / (length));
		vFaces[i]->fnz = ((Nz) / (length));
	}
}

void BuidVertexAdjacentFaceDetails() {
	for (int i = 0; i != vVertex.size(); i++)
	{
		HE_edge* outgoing_he = vVertex[i]->edge;
		HE_edge* curr = outgoing_he;

		vVertex[i]->Areaofadjacentindividualface.push_back(vFaces[(curr->face->fid) - 1]->Area);
		vVertex[i]->adjacentfaceid.push_back((curr->face->fid));

		bool found = 0;
		while (curr != NULL && curr->pair != NULL && curr->pair->next != NULL && curr->pair->next != outgoing_he)
		{
			found = 1;
			curr = curr->pair->next;
			vVertex[i]->adjacentfaceid.push_back((curr->face->fid));
			vVertex[i]->Areaofadjacentindividualface.push_back(vFaces[(curr->face->fid) - 1]->Area);
		}

		if (found == 0)
		{
			while (curr != NULL && curr->prev != NULL && curr->prev->pair != NULL && curr->prev->pair != outgoing_he)
			{
				curr = curr->prev->pair;
				vVertex[i]->adjacentfaceid.push_back((curr->face->fid));
				vVertex[i]->Areaofadjacentindividualface.push_back(vFaces[(curr->face->fid) - 1]->Area);
			}
		}

		double temparea = 0;

		for (int j = 0; j != vVertex[i]->Areaofadjacentindividualface.size(); j++)
		{
			temparea = temparea + vVertex[i]->Areaofadjacentindividualface[j];
		}

		vVertex[i]->FArea = temparea;
	}
}

void calculateVertexNormal() {
	for (int i = 0; i != vVertex.size(); i++)
	{
		vector<double> vertexAdjacentFaceAreas = vVertex[i]->Areaofadjacentindividualface;
		double vertexFaceArea = vVertex[i]->FArea;
		for (int j = 0; j != vertexAdjacentFaceAreas.size(); j++)
		{
			vVertex[i]->vnx = vVertex[i]->vnx + (((vertexAdjacentFaceAreas[j]) / (vertexFaceArea)) * (vFaces[(vVertex[i]->adjacentfaceid[j]) - 1]->fnx));
			vVertex[i]->vny = vVertex[i]->vny + (((vertexAdjacentFaceAreas[j]) / (vertexFaceArea)) * (vFaces[(vVertex[i]->adjacentfaceid[j]) - 1]->fny));
			vVertex[i]->vnz = vVertex[i]->vnz + (((vertexAdjacentFaceAreas[j]) / (vertexFaceArea)) * (vFaces[(vVertex[i]->adjacentfaceid[j]) - 1]->fnz));
		}
	}
}

void InitializeGlut(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Mesh Viewer");
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glutDisplayFunc(Draw);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouseButton);
	glutMotionFunc(mouseMove);

	glutMainLoop();
	system("pause");
}

int main(int argc, char** argv)
{
	char modelname = Welcome();
	vector<string> fileLines = ReadSelectedModel(modelname);
	BuildVerticesandFaces(fileLines);
	map<pair<int, int>, HE_edge*> edges = BuildFacesHalfEdgeDataStructure();
	BuildTwinEdgeDataStructure(edges);
	CalculateFaceNormal();
	BuidVertexAdjacentFaceDetails();
	calculateVertexNormal();
	InitializeGlut(argc, argv);

	return 0;
}
