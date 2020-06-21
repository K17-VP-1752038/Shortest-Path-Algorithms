#include "Graph.h"
#include <iostream>
#define MAX 999999
using namespace std;

Graph::Graph()
{
	nVertices = nEdges = 0;
	directed = 0;
	nDegs = NULL;
	adjMatrix = NULL;
}

Graph::Graph(const Graph & g)
{
	nEdges = g.nEdges;
	nVertices = g.nVertices;
	autoMemoryAllocated();
	for (int i = 0; i < nVertices; i++)
		for (int j = 0; j < nVertices; j++)
			adjMatrix[i][j] = g.adjMatrix[i][j];

	nDegs = new int*[2];
	for (int i = 0; i < 2; i++)
		nDegs[i] = new int[nVertices];
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < nVertices; j++)
			nDegs[i][i] = g.nDegs[i][j];
}

Graph & Graph::operator=(Graph g)
{
	nEdges = g.nEdges;
	nVertices = g.nVertices;
	autoMemoryAllocated();
	for (int i = 0; i < nVertices; i++)
		for (int j = 0; j < nVertices; j++)
			adjMatrix[i][j] = g.adjMatrix[i][j];

	nDegs = new int*[2];
	for (int i = 0; i < 2; i++)
		nDegs[i] = new int[nVertices];
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < nVertices; j++)
			nDegs[i][i] = g.nDegs[i][j];

	return *this;
}

int Graph::loadGraphFromAdjMatrixFile(const char* fname)
{
	ifstream file;
	file.open(fname);
	if (file.fail())
		return 0;		// khong mo duoc file

	// Doc so dinh cua do thi
	file >> nVertices;

	// Cap phat dong ma tran ke
	autoMemoryAllocated();

	// Doc ma tran ke
	for (int i = 0; i < nVertices; i++)
		for (int j = 0; j < nVertices; j++)
			file >> adjMatrix[i][j];

	file.close();

	// Kiem tra do thi co huong hay vo huong
	checkDirectedGraph();

	// Dem so canh
	nEdges = 0;
	if (directed)
	{
		for (int i = 0; i < nVertices; i++)
			for (int j = 0; j < nVertices; j++)
				if (adjMatrix[i][j])
					nEdges++;
	}
	else
	{
		for (int i = 0; i < nVertices; i++)
			for (int j = i + 1; j < nVertices; j++)
				if (adjMatrix[i][j])
					nEdges++;
	}

	// Dem so bac cua dinh
	insertDegs();

	return 1;
}

// Ham doc danh sach ke tu file
int Graph::loadGraphFromAdjListFile(const char * fname)
{
	ifstream file;
	file.open(fname);
	if (file.fail())
		return 0;

	// Doc so dinh cua do thi
	file >> nVertices;
	
	// Cap phat dong ma tran ke
	autoMemoryAllocated();

	// Set ma tran ke = 0
	for (int i = 0; i < nVertices; i++)
		for (int j = 0; j < nVertices; j++)
			adjMatrix[i][j] = 0;

	// Doc danh sach ke vao ma tran ke
	int size;
	int data;
	for (int i = 0; i < nVertices; i++)
	{
		file >> size;
		for (int j = 0; j < size; j++)
		{
			file >> data;
			adjMatrix[i][data] = 1;
		}
	}
	file.close();
	
	// Kiem tra do thi co huong hay vo huong
	checkDirectedGraph();

	// Dem so canh
	nEdges = 0;
	if (directed)
	{
		for (int i = 0; i < nVertices; i++)
			for (int j = 0; j < nVertices; j++)
				if (adjMatrix[i][j])
					nEdges++;
	}
	else
	{
		for (int i = 0; i < nVertices; i++)
			for (int j = i + 1; j < nVertices; j++)
				if (adjMatrix[i][j])
					nEdges++;
	}

	// Dem so bac cua dinh
	insertDegs();

	return 1;
}

// Ham kiem tra do thi co huong/vo huong
void Graph::checkDirectedGraph()
{
	directed = false;
	for (int i = 0; i < nVertices; i++)
		for (int j = i + 1; j < nVertices; j++)
			if (adjMatrix[i][j] != adjMatrix[j][i])
			{
				directed = true;
				break;
			}
}

// Ham cap phat bo nho dong
void Graph::autoMemoryAllocated()
{
	adjMatrix = new int*[nVertices];
	for (int i = 0; i < nVertices; i++)
		adjMatrix[i] = new int[nVertices];
}

// Ham doc bac cua dinh
void Graph::insertDegs()
{
	// Cap phat dong mang bac cua dinh
	nDegs = new int*[2];
	for (int i = 0; i < 2; i++)
		nDegs[i] = new int[nVertices];
	
	// Luu lai so bac cua dinh trong do thi vo huong/co huong
	if (directed)
	{
		for (int i = 0; i < nVertices; i++)
		{
			int inDeg = 0, outDeg = 0;
			for (int j = 0; j < nVertices; j++)
			{
				if (adjMatrix[i][j])
					outDeg++;
				if (adjMatrix[j][i])
					inDeg++;
			}
			nDegs[0][i] = outDeg;
			nDegs[1][i] = inDeg;
		}
	}
	else
	{
		for (int i = 0; i < nVertices; i++)
		{
			int Deg = 0;
			for (int j = 0; j < nVertices; j++)
				if (adjMatrix[i][j])
					Deg++;
			nDegs[0][i] = Deg;
		}
	}
}

// ham kiem tra do thi phan doi du
bool Graph::completeBipartiteGraph()
{
	// tim ra 2 loai bac cua dinh
	int degX = nDegs[0][0], degY;
	for (int i = 1; i < nVertices; i++)
		if (nDegs[0][i] != degX)
		{
			degY = nDegs[0][i];
			break;
		}
	
	// ton tai dinh co bac khac X, Y?
	for (int i = 0; i < nVertices; i++)
		if (nDegs[0][i] != degX && nDegs[0][i] != degY)
			return false;
	// phan loai 2 tap dinh X, Y
	int *X = new int[degY];
	int *Y = new int[degX];
	for (int i = 0; i < nVertices; i++)
	{
		if (nDegs[0][i] == degX)
		{
			*X = i;
			X++;
		}
		if (nDegs[0][i] == degX)
		{
			*Y = i;
			Y++;
		}
	}
	// kiem tra tap X
	for (int i = 0; i < degY; i++)
		for (int j = i; j < degY; j++)
			if (adjMatrix[X[i]][X[j]])
				return false;
	// kiem tra tap Y
	for (int i = 0; i < degY; i++)
		for (int j = i; j < degY; j++)
			if (adjMatrix[X[i]][X[j]])
				return false;

	delete[] X; delete[] Y;
	return true;
}

// ham kiem tra do thi day du
bool Graph::completeGraph()
{
	for (int i = 0; i < nVertices; i++)
		if (adjMatrix[i][i])
			return false;
	if (nEdges == (nVertices*(nVertices - 1)) / 2)
		return true;
	return false;
}

// tao do thi vo huong nen cho do thi
Graph Graph::createUnderUndirGraph()
{
	Graph tmp(*this);
	if (directed)
	{
		for (int i = 0; i < tmp.nVertices; i++)
			for (int j = 0; j < tmp.nVertices; j++)
			{
				if(tmp.adjMatrix[i][j])
					tmp.adjMatrix[j][i] = tmp.adjMatrix[i][j];
				if (tmp.adjMatrix[j][i])
					tmp.adjMatrix[i][j] = tmp.adjMatrix[j][i];
			}
	}
	return tmp;
}

// tao do thi bu cho do thi vo huong
void Graph::createCompGraph()
{
	if (!directed)
	{
		Graph tmp;
		// Cap phat dong ma tran ke
		tmp.nVertices = nVertices;
		tmp.autoMemoryAllocated();
		
		// Sao chep ma tran ke
		tmp.nVertices = nVertices;
		for (int i = 0; i < nVertices; i++)
			for (int j = 0; j < nVertices; j++)
				tmp.adjMatrix[i][j] = adjMatrix[i][j];

		for (int i = 0; i < nVertices; i++)
			for (int j = 0; j < nVertices; j++)
			{
				if (adjMatrix[i][j])
					adjMatrix[i][j] = 0;
				else
					adjMatrix[i][j] = 1;
			}
		cout << "Ma tran ke cua do thi bu:" << endl;
		tmp.printMatrix();
	}
}

// tao do thi nguoc cho do thi co huong
void Graph::createConvGraph()
{
	if (directed)
	{
		Graph tmp;
		// Cap phat dong ma tran ke
		tmp.nVertices = nVertices;
		tmp.autoMemoryAllocated();

		// sao chep nguoc ma tran ke
		tmp.nVertices = nVertices;
		for (int i = 0; i < nVertices; i++)
			for (int j = 0; j < nVertices; j++)
				tmp.adjMatrix[j][i] = adjMatrix[i][j];

		cout << "Ma tran ke cua do thi nguoc:" << endl;
		tmp.printMatrix();
	}
}

// ham in ma tran
void Graph::printMatrix()
{
	for (int i = 0; i < nVertices; i++)
	{
		for (int j = 0; j < nVertices; j++)
			cout << adjMatrix[i][j] << " ";
		cout << endl;
	}
}

// Ham xuat mot vai thong tin cua do thi
void Graph::printGraph()
{
	if (directed)
		cout << "Do thi co huong." << endl;
	else
		cout << "Do thi vo huong." << endl;
	
	cout << "So dinh: " << nVertices << endl;
	cout << "So canh: " << nEdges << endl;

	printMatrix();

	if (directed)	//DO THI CO HUONG
	{
		// So bac cua dinh
		cout << "\t Bac vao     Bac ra" << endl;
		for (int j = 0; j < nVertices; j++)
			cout << "Dinh " << j << ":     " << nDegs[0][j] << "\t\t" << nDegs[1][j] << endl;

		// Danh sach dinh co lap
		cout << "--Danh sach dinh co lap--" << endl;
		bool nonVertex = false;
		for (int i = 0; i < nVertices; i++)
			if (nDegs[0][i] + nDegs[1][i] == 0)
			{
				nonVertex = true;
				cout << "Dinh " << i << endl;
			}
		if (!nonVertex)
			cout << "NULL" << endl;

		// Danh sach dinh treo
		cout << "--Danh sach dinh treo--" << endl;
		bool externalVertex = false;
		for (int i = 0; i < nVertices; i++)
			if (nDegs[0][i] + nDegs[1][i] == 1)
			{
				externalVertex = true;
				cout << "Dinh " << i << endl;
			}
		if (!externalVertex)
			cout << "NULL" << endl;
	}
	else  //DO THI VO HUONG
	{
		// So bac cua dinh
		cout << "\t Bac" << endl;
		for (int i = 0; i < nVertices; i++)
			cout << "Dinh " << i << ":    " << nDegs[0][i] << endl;

		// Danh sach dinh co lap
		cout << "--Danh sach dinh co lap--" << endl;
		bool nonVertex = false;
		for (int i = 0; i < nVertices; i++)
			if (nDegs[0][i] == 0)
			{
				nonVertex = true;
				cout << "Dinh " << i << endl;
			}
		if (!nonVertex)
			cout << "NULL" << endl;

		// Danh sach dinh treo
		cout << "--Danh sach dinh treo--" << endl;
		bool externalVertex = false;
		for (int i = 0; i < nVertices; i++)
			if (nDegs[0][i] == 1)
			{
				externalVertex = true;
				cout << "Dinh " << i << endl;
			}
		if (!externalVertex)
			cout << "NULL" << endl;

		// Kiem tra do thi day du
		if (completeGraph())
			cout << "=> Do thi day du" << endl;

		// Kiem tra do thi vong
		for (int i = 0; i < nVertices; i++)
			if (nDegs[0][i] != 2)
				goto next1;
		cout << "=> Do thi vong" << endl;
	next1:
		// Kiem tra do thi banh xe
		int degn = 0, deg3 = 0;
		for (int i = 0; i < nVertices; i++)
		{
			if (nDegs[0][i] == nVertices - 1)
				degn++;
			if (nDegs[0][i] == 3)
				deg3++;
		}
		if (degn == 1 && deg3 == nVertices - 1)
			cout << "=> Do thi banh xe" << endl;
		
		// Kiem tra do thi phan doi du
		if (completeBipartiteGraph())
			cout << "=> Do thi phan doi du" << endl;
	}
}

// Ham duyet sau do thi
void Graph::DFS(int v)
{
	stack<int> S;
	bool *visited = new bool[nVertices];
	for (int i = 0; i < nVertices; i++)
		visited[i] = 0;
	// Danh dau dinh da tham
	visited[v] = 1;
	// In ra dinh
	cout << v << " ";
	// Them dinh vao ngan xep
	S.push(v);

	while (!S.empty())
	{
		int unvisitedVertex = getUnvisitedVertex(S.top(), visited);
		// Khong tim thay dinh can tham
		if (unvisitedVertex == -1)
			S.pop();
		else
		{
			visited[unvisitedVertex] = 1;
			cout << unvisitedVertex << " ";
			S.push(unvisitedVertex);
		}
	}
	delete []visited;
}

// Ham lay dinh lien thong chua tham
int Graph::getUnvisitedVertex(int vertex, bool *&visited)
{
	Graph g = createUnderUndirGraph();
	for (int i = 0; i < g.nVertices; i++)
		if (g.adjMatrix[vertex][i] && !visited[i])
			return i;
	return -1;
}

bool Graph::hasNegativeEdge()
{
	for (int i = 0; i < nVertices; i++)
		for (int j = 0; j < nVertices; j++)
			if (adjMatrix[i][j] < 0)
				return true;
	return false;
}

// Ham duyet rong do thi
void Graph::BFS(int v)
{
	queue<int> Q;
	bool *visited = new bool[nVertices];
	for (int i = 0; i < nVertices; i++)
		visited[i] = 0;
	// Danh dau dinh
	visited[v] = 1;
	// In ra dinh
	cout << v << " ";
	// Them dinh vao hang doi
	Q.push(v);
	
	while (!Q.empty())
	{
		int visitedVertex = Q.front();
		Q.pop();
		while(getUnvisitedVertex(visitedVertex, visited) != -1)
		{
			int i = getUnvisitedVertex(visitedVertex, visited);
			visited[i] = 1;
			cout << i << " ";
			Q.push(i);
		}
	}
	delete []visited;
}

// Ham kiem tra duong di tu 1 dinh toi 1 dinh
bool Graph::reachable(int v1, int v2)
{
	bool flag = false;
	stack<int> S;
	bool *visited = new bool[nVertices];
	for (int i = 0; i < nVertices; i++)
		visited[i] = 0;
	visited[v1] = 1;
	S.push(v1);
	while (!S.empty())
	{
		int unvisitedVertex = getUnvisitedVertex(S.top(), visited);
		if (unvisitedVertex == -1)
			S.pop();
		else if (unvisitedVertex == v2)
		{
			flag = true;
			break;
		}
		else
		{
			visited[unvisitedVertex] = 1;
			S.push(unvisitedVertex);
		}
	}
	delete[]visited;
	return flag;
}

// Ham kiem tra do thi lien thong
bool Graph::connected()
{
	if (countConnectedCompo() == 1)
		return true;
	return false;
}

// Ham kiem tra thanh phan lien thong
void Graph::connectedComponents()
{
	stack<int> S;
	bool *visited = new bool[nVertices];
	for (int i = 0; i < nVertices; i++)
		visited[i] = 0;
	for (int i = 0; i < nVertices;i++)
		if (!visited[i])
		{
			visited[i] = 1;
			cout << i << " ";
			S.push(i);

			while (!S.empty()) {
				int unvisitedVertex = getUnvisitedVertex(S.top(), visited);
				if (unvisitedVertex == -1)
					S.pop();
				else {
					visited[unvisitedVertex] = 1;
					cout << unvisitedVertex << " ";
					S.push(unvisitedVertex);
				}
			}
			cout << endl;
		}
	delete[]visited;
}

// Ham dem so thanh phan lien thong cua do thi
int Graph::countConnectedCompo()
{
	int count = 0;
	stack<int> S;
	bool *visited = new bool[nVertices];
	for (int i = 0; i < nVertices; i++)
		visited[i] = 0;
	for (int i = 0; i < nVertices;i++)
		if (!visited[i]) 
		{
			count++;
			visited[i] = 1;
			S.push(i);
			while (!S.empty())
			{
				int unvisitedVertex = getUnvisitedVertex(S.top(), visited);
				if (unvisitedVertex == -1)
					S.pop();
				else {
					visited[unvisitedVertex] = 1;
					S.push(unvisitedVertex);
				}
			}
		}
	delete[]visited;
	return count;
}

// ham xac dinh canh cau
void Graph::cutEdges()
{
	int ccc = countConnectedCompo();
	for (int i = 0; i < nVertices; i++)
	{
		for (int j = i + 1; j < nVertices; j++) 
		{
			if (adjMatrix[i][j]) 
			{
				Graph G = createUnderUndirGraph();
				G.adjMatrix[i][j] = 0;
				if(G.countConnectedCompo() > ccc)
					cout << "[" << i << "-" << j << "]" << " ";
			}
		}
	}
}

// Ham xuat cay khung theo dfs
void Graph::DFSSpanTree(int v)
{
	stack<int> S;
	bool *visited = new bool[nVertices];
	for (int i = 0; i < nVertices; i++)
		visited[i] = 0;
	// Danh dau dinh da tham
	visited[v] = 1;
	// Them dinh vao ngan xep
	S.push(v);

	while (!S.empty())
	{
		int unvisitedVertex = getUnvisitedVertex(S.top(), visited);
		if (unvisitedVertex == -1)
			S.pop();
		else
		{
			visited[unvisitedVertex] = 1;
			cout << v << " - " << unvisitedVertex << endl;
			v = unvisitedVertex;
			S.push(unvisitedVertex);
		}
	}
	delete[]visited;
}

// Ham xuat cay khung theo bfs
void Graph::BFSSpanTree(int v)
{
	queue<int> Q;
	bool *visited = new bool[nVertices];
	for (int i = 0; i < nVertices; i++)
		visited[i] = 0;
	visited[v] = 1;
	Q.push(v);

	while (!Q.empty())
	{
		int visitedVertex = Q.front();
		Q.pop();
		while (getUnvisitedVertex(visitedVertex, visited) != -1)
		{
			int i = getUnvisitedVertex(visitedVertex, visited);
			visited[i] = 1;
			cout << visitedVertex << " - " << i << endl;
			Q.push(i);
		}
	}
	delete[]visited;
}

void Graph::PrimMinSpanTree(int v)
{
	int count = 0;
	int weight = 0;
	bool *visited = new bool[nVertices];
	for (int i = 0; i < nVertices; i++)
		visited[i] = 0;
	visited[v] = 1;

	int x = 0, y = 0;
	while(count < nVertices-1)
	{
		int min = MAX;
		for (int i = 0; i < nVertices; i++)
			if (visited[i])
				for(int j = 0; j < nVertices; j++)
					if(!visited[j])
						if (min > adjMatrix[i][j] && adjMatrix[i][j])
						{
							x = i;
							y = j;
							min = adjMatrix[i][j];
						}
		count++;
		visited[y] = 1;
		cout << x << " - " << y  << " : " << min << endl;
		weight += min;
	}
	cout << "Tong trong so = " << weight << endl;
	delete[] visited;
}

// Ham tim cay khung toi thieu theo Kruskal
void Graph::KruskalMinSpanTree()
{
	int count = 0, weight = 0;
	int *check = new int[nVertices];
	for (int i = 0; i < nVertices; i++)
		check[i] = i;
	
	while (count < nVertices - 1)
	{
		int min = MAX;
		int x = 0, y = 0;
		for (int i = 0; i < nVertices; i++)
			for (int j = i + 1; j < nVertices; j++)
				if (check[i] != check[j])
					if (min > adjMatrix[i][j] && adjMatrix[i][j])
					{
						x = i;
						y = j;
						min = adjMatrix[i][j];
					}
		cout << x << " - " << y << " : "<< min << endl;
		weight += min;
		count++;
		join(x, y, check);
	}
	cout << "Tong trong so: " << weight << endl;
	delete[]check;
}

// Ham kiem tra 2 dinh cua thuat toan kruskal
void Graph::join(int x, int y, int *&check)
{
	int small = (check[x] < check[y]) ? check[x] : check[y];
	int big = (check[x] > check[y]) ? check[x] : check[y];
	for (int i = 0; i < nVertices; i++)
		if (check[i] == small)
			check[i] = big;
}

void Graph::DijkstraShortestPaths(int v1, int v2)
{
	if (hasNegativeEdge())
		return;
	Vertex *A = new Vertex[nVertices];
	bool *close = new bool[nVertices];
	close[v1] = 1;
	for (int i = 0; i < nVertices; i++)
	{
		close[i] = 0;
		A[i].dist = MAX;
		A[i].lastV = -1;
	}
	A[v1].dist = 0;

	for (int i = 0; i < nVertices; i++)
	{
		int ver = findVertex(v1, A, close);
		close[ver] = 1;
		for(int j = 0; j < nVertices; j++)
			if(adjMatrix[ver][j] && !close[j])
			{
				int distand = A[ver].dist + adjMatrix[ver][j];
				if (A[j].dist == MAX || A[j].dist > distand)
				{
					A[j].dist = distand;
					A[j].lastV = ver;
				}
			}
	}
	PrintPath(A, v1, v2);
	cout << " : " << A[v2].dist;
	delete[]A;
}

void Graph::PrintPath(Vertex *A, int vertex1, int vertex2) {
	if(vertex1 == vertex2) 
		cout << vertex1;
	else 
	{
		int temp = A[vertex2].lastV;
		PrintPath(A, vertex1, temp);
		cout << " -> " << vertex2;
	}
}

int Graph::findVertex(int v, Vertex *&A, bool *close)
{
	int min = MAX, tmp;
	for (int i = 0; i < nVertices; i++)
	{
		if (min > A[i].dist && !close[i]) {
			min = A[i].dist;
			tmp = i;
		}
	}
	return tmp;
}

bool Graph::compareArray(Vertex *a, Vertex *b) 
{
	for (int i = 0; i < nVertices; i++)
		if (a[i].dist != b[i].dist)
			return false;
	return true;
}

void Graph::BellmanShortestPaths(int vertex1, int vertex2)
{
	Vertex *A = new Vertex[nVertices];
	Vertex *newA = new Vertex[nVertices];

	for (int i = 0; i < nVertices; i++) 
	{
		A[i].dist = MAX;
		A[i].lastV = -1;
	}

	A[vertex1].dist = 0;
	int k = 0;
	bool flag = 0;

	while (k <= nVertices) 
	{
		for (int i = 0; i < nVertices; i++) 
		{
			newA[i] = A[i];
			for(int j = 0; j < nVertices; j++)
				if (adjMatrix[j][i] != 0) 
				{
					Vertex temp;
					temp.dist = A[j].dist + adjMatrix[j][i];
					temp.lastV = j;
					if (temp.dist < newA[i].dist) 
					{
						newA[i] = temp;
						if (compareArray(newA, A)) {
							flag = 1;
							break;
						}
						A[i] = newA[i];
					}
				}
		}
		if (flag == 1)
			break;
		k++;
	}
	if (k == nVertices)
		cout << "Chu trinh am (-)!" << endl;
	else {
		PrintPath(newA, vertex1, vertex2);
		cout << " : " << newA[vertex2].dist << endl;
	}
	delete[]A;
	delete[]newA;
}

// Huy cap phat bo nho
Graph::~Graph()
{
	for (int i = 0; i < nVertices; i++)
		delete[] adjMatrix[i];
	delete[] adjMatrix;
	delete[] nDegs[0]; delete[] nDegs[1];
	delete[] nDegs;
}
