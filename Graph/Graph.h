#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <queue>
#include <stack>


using namespace std;

const int int_max = numeric_limits<int>::max();
const int inf = int_max >> 1;

typedef void* AdjacencyMatrixTag;
typedef void** IncidenceMatrixTag;
typedef void*** AdjacencyListTag;


class Graph
{
	vector <vector <int>> AdjacencyMatrix;
	int Vertices;
	int Edges;
	bool Directed;

	bool BFSForMaximalFlow(int s, int t, vector <vector <int>>& f, vector <int>& d);

	int DFSForMaximalFlow(int u, int minC, int s, int t, vector<int>& p, vector <vector <int>>& f, vector <int>& d);

	vector <int> BFS(int source, int destination, bool all);

	vector <int> top_sort();

public:
	Graph();

	Graph(AdjacencyMatrixTag);

	Graph(IncidenceMatrixTag);

	void PrintAdjacencyMatrix();

	void PrintIncidenceMatrix();

	void PrintAdjacencyList();

	void PrintDegreeSequence();

	void PrintHangingAndIsolatedVertices();

	void PrintSourcesAndDrains();

	void AddVertex();

	void RemoveVertex(int n);

	void AddEdge(int i, int j);

	void RemoveEdge(int i, int j);

	Graph Adjunction();

	static Graph Conjunction(Graph G1, Graph G2);

	static Graph Disjunction(Graph G1, Graph G2);

	static Graph DisjunctDisjunction(Graph G1, Graph G2);

	static Graph Addition(Graph G1, Graph G2);

	static Graph Multiplication(Graph G1, Graph G2);

	void VertexDoubling(int n);

	void VertexReproduction(int n);

	void EdgeSubdivision(int i, int j);

	void VerticesIdentification(int n, int m);

	void PrintSpanningTreeUsingDFS();

	void PrintSpanningTreeUsingBFS();

	void PrintMinimumWeightedSpanningTree();

	void GetAllDistancesUsingFloyd();

	void GetDistancesInContourlessNetwork(int v1);

	void FindMaxFlow(int s, int t);

	int GetDistance(int n, int m);

	int GetEccentricity(int n);

	int GetRadius();

	int GetDiameter();

	void PrintPeripheralVertices();

	void PrintCentralVertices();

	int GetDegree(int n);

	int GetInSemidegree(int n);

	int GetOutSemidegree(int n);

	int GetNumberOfVertices();

	int GetNumberOfEdges();

	bool GetOrientation();
};