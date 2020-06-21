#pragma once
#include <fstream>
#include <stack>
#include <queue>
using namespace std;

struct Vertex
{
	int dist;
	int lastV;
};

class Graph
{
private:
	int nVertices;		// so dinh
	int nEdges;			// so canh
	bool directed;		// do thi co huong/vo huong
	int **nDegs;		// so bac cua dinh
	int **adjMatrix;	// ma tran ke (cap phat dong)
	
	void checkDirectedGraph();
	void autoMemoryAllocated();	// cap phat dong ma tran ke
	void insertDegs();			// luu so bac cua dinh
	int getUnvisitedVertex(int vertex, bool *&visited);
	bool hasNegativeEdge();		// co canh trong so am
public:
	Graph();
	Graph(const Graph &g);
	Graph& operator =(Graph g);

	int loadGraphFromAdjMatrixFile(const char * fname);
	int loadGraphFromAdjListFile(const char * fname);
	
	bool completeBipartiteGraph();
	bool completeGraph();			// kiem tra do thi day du
	Graph createUnderUndirGraph();	// do thi vo huong nen
	void createCompGraph();			// do thi bu (vo huong)
	void createConvGraph();			// do thi nguoc (co huong)
	void printMatrix();
	void printGraph();

	void DFS(int v);				// ham duyet sau
	void BFS(int v);				// ham duyet rong
	bool reachable(int v1, int v2); // ham kiem tra duong di
	bool connected();				// ham kiem tra do thi lien thong
	void connectedComponents();		// ham in cac thanh phan lien thong
	int countConnectedCompo();		// ham dem so thanh phan lien thong
	void cutEdges();				// ham xac dinh canh cầu
	
	void DFSSpanTree(int v);		// Ham tim cay khung theo dfs
	void BFSSpanTree(int v);		// Ham tim cay khung theo bfs
	void PrimMinSpanTree(int v);	// Ham tim cay khung toi thieu theo Prim
	void KruskalMinSpanTree();		// Ham tim cay khung toi thieu theo Kruskal
	void join(int a, int b, int *&check);

	int findVertex(int v, Vertex *&A, bool *close);		// ham tim dinh co duong di den ngan nhat
	void DijkstraShortestPaths(int v1, int v2);			// ham tim duong di ngan nhat theo Djikstra
	void PrintPath(Vertex *A, int vertex1, int vertex2);// ham in duong di
	void BellmanShortestPaths(int vertex1, int vertex2);// ham tim duong di ngan nhat theo Bellman
	bool compareArray(Vertex *a, Vertex *b);			// ham so sanh hai mang

	~Graph();
};
