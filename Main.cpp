#include <conio.h>
#include <iostream>
#include "Graph.h"

int main()
{
	Graph G;
	if (G.loadGraphFromAdjMatrixFile("Graph3.txt"))
	{
		cout << "Find distand from vertex 0 to others:" << endl;
		G.DijkstraShortestPaths(0, 3);
		cout << endl;
		G.BellmanShortestPaths(0, 3);
	}
	else
		cout << "Khong tao duoc do thi!" << endl;
	
	/*
	if (G.loadGraphFromAdjListFile("Graph2.txt"))
	{
		G.printGraph();
		G.createUnderUndirGraph();
	}
	else
		cout << "Khong tao duoc do thi!" << endl;
	*/

	_getch();
	return 0;
}