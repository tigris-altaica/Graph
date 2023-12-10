#include "Graph.h"

int main()
{
	void* type = 0;

	Graph* G = nullptr;
	string command;

	while (true)
	{
		getline(cin, command);
		if (command == "gr am")
		{
			G = new Graph((AdjacencyMatrixTag)type);
		}
		else if (command == "gr im")
		{
			G = new Graph((IncidenceMatrixTag)type);
		}
		else if (command == "get v")
		{
			cout << "Vertices: " << G->GetNumberOfVertices() << endl;
		}
		else if (command == "get e")
		{
			cout << "Edges: " << G->GetNumberOfEdges() << endl;
		}
		else if (command == "get deg")
		{
			int n;
			cout << "Vertex N ";
			cin >> n;
			if (!G->GetOrientation())
				cout << "Degree: " << G->GetDegree(n) << endl;
			else
			{
				cout << "In Degree: " << G->GetInSemidegree(n) << endl;
				cout << "Out Degree: " << G->GetOutSemidegree(n) << endl;
			}
		}
		else if (command == "deg seq")
		{
			G->PrintDegreeSequence();
		}
		else if (command == "pr am")
		{
			G->PrintAdjacencyMatrix();
		}
		else if (command == "pr im")
		{
			G->PrintIncidenceMatrix();
		}
		else if (command == "pr al")
		{
			G->PrintAdjacencyList();
		}
		else if (command == "h i")
		{
			G->PrintHangingAndIsolatedVertices();
		}
		else if (command == "s d")
		{
			G->PrintSourcesAndDrains();
		}
		else if (command == "dist")
		{
			int n, m, d;
			cout << "From: ";
			cin >> n;
			cout << "To: ";
			cin >> m;
			d = G->GetDistance(n, m);
			if (d != int_max)
				cout << "Distance: " << d << endl;
			else
				cout << "INFINITY" << endl;
		}
		else if (command == "ec")
		{
			int n, d;
			cout << "Vertex N ";
			cin >> n;
			d = G->GetEccentricity(n);
			if (d != int_max)
				cout << "Eccentricity: " << d << endl;
			else
				cout << "INFINITY" << endl;
		}
		else if (command == "rad")
		{
			int d = G->GetRadius();
			if (d != int_max)
				cout << "Radius: " << d << endl;
			else
				cout << "INFINITY" << endl;
		}
		else if (command == "dia")
		{
			int d = G->GetDiameter();
			if (d != int_max)
				cout << "Diameter: " << d << endl;
			else
				cout << "INFINITY" << endl;
		}
		else if (command == "per")
		{
			G->PrintPeripheralVertices();
		}
		else if (command == "cen")
		{
			G->PrintCentralVertices();
		}
		else if (command == "add v")
		{
			G->AddVertex();
			G->PrintAdjacencyMatrix();
		}
		else if (command == "rem v")
		{
			int n;
			cout << "Vertex N ";
			cin >> n;
			G->RemoveVertex(n);
			G->PrintAdjacencyMatrix();
		}
		else if (command == "add e")
		{
			int i, j;
			cout << "From: ";
			cin >> i;
			cout << "To: ";
			cin >> j;
			G->AddEdge(i, j);
			G->PrintAdjacencyMatrix();
		}
		else if (command == "rem e")
		{
			int i, j;
			cout << "From: ";
			cin >> i;
			cout << "To: ";
			cin >> j;
			G->RemoveEdge(i, j);
			G->PrintAdjacencyMatrix();
		}
		else if (command == "adj")
		{
			Graph H = G->Adjunction();
			H.PrintAdjacencyMatrix();
		}
		else if (command == "e sub")
		{
			int i, j;
			cout << "From: ";
			cin >> i;
			cout << "To: ";
			cin >> j;
			G->EdgeSubdivision(i, j);
			G->PrintAdjacencyMatrix();
		}
		else if (command == "v id")
		{
			int n, m;
			cout << "Vertex 1 N ";
			cin >> n;
			cout << "Vertex 2 N ";
			cin >> m;
			G->VerticesIdentification(n, m);
			G->PrintAdjacencyMatrix();
		}
		else if (command == "v dou")
		{
			int n;
			cout << "Vertex N ";
			cin >> n;
			G->VertexDoubling(n);
			G->PrintAdjacencyMatrix();
		}
		else if (command == "v rep")
		{
			int n;
			cout << "Vertex N ";
			cin >> n;
			G->VertexReproduction(n);
			G->PrintAdjacencyMatrix();
		}
		else if (command == "spt dfs")
		{
			G->PrintSpanningTreeUsingDFS();
		}
		else if (command == "spt bfs")
		{
			G->PrintSpanningTreeUsingBFS();
		}
		else if (command == "min spt")
		{
			G->PrintMinimumWeightedSpanningTree();
		}
		else if (command == "all dist")
		{
			G->GetAllDistancesUsingFloyd();
		}
		else if (command == "from dist")
		{
			int n;
			cout << "Vertex N ";
			cin >> n;
			G->GetDistancesInContourlessNetwork(n - 1);
		}
		else if (command == "max flow")
		{
			G->PrintSourcesAndDrains();
			int s, t;
			cout << "S ";
			cin >> s;
			cout << "T ";
			cin >> t;
			G->FindMaxFlow(s - 1, t - 1);
		}
		else if (command == "ex")
		{
			if (G != nullptr) {
				delete G;
			}
			break;
		}
	}
}