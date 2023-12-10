#include "Graph.h"


bool Graph::BFSForMaximalFlow(int s, int t, vector <vector <int>>& f, vector <int>& d)
{
	d.assign(Vertices, inf);
	d[s] = 0;
	queue <int> Q;
	Q.push(s);
	while (!Q.empty())
	{
		int u = Q.front();
		Q.pop();
		for (int v = 0; v < Vertices; v++)
		{
			if (AdjacencyMatrix[u][v] <= 0)
				continue;
			if (f[u][v] < AdjacencyMatrix[u][v] && d[v] == inf)
			{
				d[v] = d[u] + 1;
				Q.push(v);
			}
		}
	}
	return d[t] != inf;
}

int Graph::DFSForMaximalFlow(int u, int minC, int s, int t, vector<int>& p, vector <vector <int>>& f, vector <int>& d)
{
	if (u == t || minC == 0)
		return minC;
	for (int v = p[u]; v <= Vertices - 1; v++)
	{
		if (d[v] == d[u] + 1)
		{
			int delta = this->DFSForMaximalFlow(v, min(minC, AdjacencyMatrix[u][v] - f[u][v]), s, t, p, f, d);
			if (delta != 0)
			{
				f[u][v] += delta;
				f[v][u] -= delta;
				return delta;
			}
		}
		p[u]++;
	}
	return 0;
}

vector <int> Graph::BFS(int source, int destination, bool all)
{
	vector <int> d(Vertices, inf);
	d[source] = 0;
	queue <int> Q;
	Q.push(source);
	while (!Q.empty())
	{
		int u = Q.front();
		Q.pop();
		for (int v = 0; v < Vertices; v++)
		{
			if (AdjacencyMatrix[u][v] <= 0)
				continue;
			if (d[v] == inf)
			{
				d[v] = d[u] + 1;
				if (!all && v == destination)
					return d;
				Q.push(v);
			}
		}
	}
	return d;
}

vector <int> Graph::top_sort()
{
	vector <int> DegOut(Vertices, 0), Index(Vertices);

	for (int v = 0; v < Vertices; v++)
	{
		for (int w = 0; w < Vertices; w++)
		{
			if (AdjacencyMatrix[w][v] == inf || w == v)
				continue;
			DegOut[w]++;
		}
	}

	queue <int> Q;
	int number = Vertices - 1;
	for (int v = 0; v < Vertices; v++)
	{
		if (DegOut[v] == 0)
			Q.push(v);
	}

	while (!Q.empty())
	{
		int v = Q.front();
		Q.pop();
		Index[number] = v;
		number--;

		for (int w = 0; w < Vertices; w++)
		{
			if (AdjacencyMatrix[w][v] == inf || w == v)
				continue;
			DegOut[w]--;
			if (DegOut[w] == 0)
				Q.push(w);
		}
	}

	return Index;
}

Graph::Graph() {}

Graph::Graph(AdjacencyMatrixTag)
{
	ifstream AdjM("AM.txt");
	string line;
	getline(AdjM, line);

	Vertices = count(line.begin(), line.end(), ' ') + 1;
	Edges = 0;
	Directed = false;

	AdjM.seekg(0);

	for (int i = 0; i < Vertices; i++)
	{
		vector <int> row;
		for (int j = 0; j < Vertices; j++)
		{
			int num;
			string str;
			AdjM >> str;
			if (str == "inf")
				num = inf;
			else
				num = stoi(str);
			row.push_back(num);
			Edges += num;
			if (!Directed && i > j && num != AdjacencyMatrix[j][i])
				Directed = true;
		}
		AdjacencyMatrix.push_back(row);
	}

	if (!Directed)
		Edges /= 2;

	AdjM.close();
}

Graph::Graph(IncidenceMatrixTag)
{
	ifstream IncM("IM.txt");
	string line;
	getline(IncM, line);

	Vertices = 0;
	Edges = count(line.begin(), line.end(), ' ') + 1;
	Directed = false;

	vector <vector <int>> IncidenceMatrix;
	IncM.seekg(0);

	while (!IncM.eof())
	{
		vector <int> v;
		int n;
		for (int j = 0; j < Edges; j++)
		{
			IncM >> n;
			if (!Directed && n < 0)
				Directed = true;
			v.push_back(n);
		}
		IncidenceMatrix.push_back(v);
		Vertices++;
	}

	IncM.close();

	vector <int> row(Vertices, 0);
	AdjacencyMatrix.assign(Vertices, row);

	for (int j = 0; j < Edges; j++)
	{
		int fromto[2], count = 0;
		for (int i = 0; i < Vertices; i++)
		{
			if (IncidenceMatrix[i][j] == 2)
			{
				fromto[0] = i;
				fromto[1] = i;
			}
			if (Directed)
			{
				if (IncidenceMatrix[i][j] == 1)
					fromto[0] = i;
				if (IncidenceMatrix[i][j] == -1)
					fromto[1] = i;
			}
			else
			{
				if (IncidenceMatrix[i][j] == 1)
					fromto[count++] = i;
			}
		}
		AdjacencyMatrix[fromto[0]][fromto[1]]++;
		if (!Directed)
			AdjacencyMatrix[fromto[1]][fromto[0]]++;
	}
}

void Graph::PrintAdjacencyMatrix()
{
	cout << "Adjacency Matrix:" << endl;
	for (int i = 0; i < Vertices; i++)
	{
		for (int j = 0; j < Vertices; j++)
		{
			cout << AdjacencyMatrix[i][j] << " ";
		}
		cout << endl;
	}
}

void Graph::PrintIncidenceMatrix()
{
	vector <vector <int>> IncidenceMatrix;
	vector <int> row(Edges, 0);
	IncidenceMatrix.assign(Vertices, row);
	int column = 0, i, j, k;
	cout << "Incidence Matrix:" << endl;
	for (i = 0; i < Vertices; i++)
	{
		for (j = 0; j < Vertices; j++)
		{
			if (!Directed && j > i)
				break;
			for (k = 0; k < AdjacencyMatrix[i][j]; k++)
			{
				IncidenceMatrix[i][column]++;
				if (!Directed || i == j)
					IncidenceMatrix[j][column]++;
				else if (Directed)
					IncidenceMatrix[j][column]--;
				column++;
				if (i == j)
					break;
			}
		}
	}
	for (i = 0; i < Vertices; i++)
	{
		for (j = 0; j < Edges; j++)
		{
			cout << IncidenceMatrix[i][j] << " ";
		}
		cout << endl;
	}
}

void Graph::PrintAdjacencyList()
{
	cout << "Adjacency List:" << endl;
	for (int i = 0; i < Vertices; i++)
	{
		cout << "{ " << i + 1 << ": ";
		for (int j = 0; j < Vertices; j++)
		{
			for (int k = 0; k < AdjacencyMatrix[i][j] / (i == j ? 2 : 1); k++)
			{
				cout << j + 1 << " ";
			}
		}
		cout << "} ";
	}
	cout << endl;
}

void Graph::PrintDegreeSequence()
{
	if (!Directed)
	{
		vector<int> Degrees(Vertices);
		for (int i = 0; i < Vertices; i++)
		{
			Degrees[i] = GetDegree(i + 1);
		}
		sort(Degrees.begin(), Degrees.end());
		cout << "Degree Sequence:" << endl;
		for (int i = 0; i < Vertices; i++)
		{
			cout << Degrees[i] << " ";
		}
		cout << endl;
	}
	else
	{
		vector<int> InDegrees(Vertices);
		vector<int> OutDegrees(Vertices);
		for (int i = 0; i < Vertices; i++)
		{
			InDegrees[i] = GetInSemidegree(i + 1);
			OutDegrees[i] = GetOutSemidegree(i + 1);
		}
		sort(InDegrees.begin(), InDegrees.end());
		sort(OutDegrees.begin(), OutDegrees.end());
		cout << "In Degrees:" << endl;
		for (int i = 0; i < Vertices; i++)
		{
			cout << InDegrees[i] << " ";
		}
		cout << endl << "Out Degrees:" << endl;
		for (int i = 0; i < Vertices; i++)
		{
			cout << OutDegrees[i] << " ";
		}
		cout << endl;
	}
}

void Graph::PrintHangingAndIsolatedVertices()
{
	vector <int> adj(Vertices, 0);

	for (int i = 0; i < Vertices; i++)
	{
		for (int j = 0; j < Vertices; j++)
		{
			if (AdjacencyMatrix[i][j] != 0 && i != j)
			{
				adj[i]++;
				if (adj[i] > 1)
					break;
			}
		}
	}
	cout << "Hanging:" << endl << "{ ";
	for (int i = 0; i < Vertices; i++)
	{
		if (adj[i] == 1)
			cout << i + 1 << " ";
	}
	cout << "}" << endl << "Isolated:" << endl << "{ ";
	for (int i = 0; i < Vertices; i++)
	{
		if (adj[i] == 0)
			cout << i + 1 << " ";
	}
	cout << "}" << endl;
}

void Graph::PrintSourcesAndDrains()
{
	if (!Directed)
		return;
	cout << "Sources:" << endl << "{ ";
	for (int j = 0; j < Vertices; j++)
	{
		for (int i = 0; i < Vertices; i++)
		{
			if (AdjacencyMatrix[i][j] != 0)
				break;
			if (i == Vertices - 1)
				cout << j + 1 << " ";
		}
	}
	cout << "}" << endl << "Drains:" << endl << "{ ";
	for (int i = 0; i < Vertices; i++)
	{
		for (int j = 0; j < Vertices; j++)
		{
			if (AdjacencyMatrix[i][j] != 0)
				break;
			if (j == Vertices - 1)
				cout << i + 1 << " ";
		}
	}
	cout << "}" << endl;
}

void Graph::AddVertex()
{
	vector <int> row(Vertices + 1, 0);
	AdjacencyMatrix.push_back(row);
	for (int i = 0; i < Vertices; i++)
	{
		AdjacencyMatrix[i].push_back(0);
	}
	Vertices++;
}

void Graph::RemoveVertex(int n)
{
	if (n > Vertices)
		return;
	if (!Directed)
		Edges -= GetDegree(n);
	else
		Edges -= GetInSemidegree(n) + GetOutSemidegree(n);
	for (int i = 0; i < Vertices; i++)
	{
		AdjacencyMatrix[i].erase(AdjacencyMatrix[i].begin() + n - 1);
	}
	AdjacencyMatrix.erase(AdjacencyMatrix.begin() + n - 1);
	Vertices--;
}

void Graph::AddEdge(int i, int j)
{
	if (i > Vertices || j > Vertices)
		return;
	AdjacencyMatrix[i - 1][j - 1]++;
	if (!Directed)
		AdjacencyMatrix[j - 1][i - 1]++;
	Edges++;
}

void Graph::RemoveEdge(int i, int j)
{
	if (i > Vertices || j > Vertices)
		return;
	AdjacencyMatrix[i - 1][j - 1]--;
	if (!Directed)
		AdjacencyMatrix[j - 1][i - 1]--;
	Edges--;
}

Graph Graph::Adjunction()
{
	Graph G;
	G.Vertices = Vertices;
	G.Edges = Edges;
	G.Directed = false;

	for (int i = 0; i < G.Vertices; i++)
	{
		vector <int> row;
		for (int j = 0; j < G.Vertices; j++)
		{
			if (i == j)
				row.push_back(0);
			else
				row.push_back(!AdjacencyMatrix[i][j]);
		}
		G.AdjacencyMatrix.push_back(row);
	}
	return G;
}

Graph Graph::Conjunction(Graph G1, Graph G2)
{
	Graph G3;
	G3.Vertices = min(G1.Vertices, G2.Vertices);
	G3.Edges = 0;
	G3.Directed = G1.Directed;

	for (int i = 0; i < G3.Vertices; i++)
	{
		vector <int>row;
		for (int j = 0; j < G3.Vertices; j++)
		{
			int x = min(G1.AdjacencyMatrix[i][j], G2.AdjacencyMatrix[i][j]);
			row.push_back(x);
			G3.Edges += x;
		}
		G3.AdjacencyMatrix.push_back(row);
	}

	if (!G3.Directed)
		G3.Edges /= 2;

	return G3;
}

Graph Graph::Disjunction(Graph G1, Graph G2)
{
	Graph G3;
	G3.Vertices = max(G1.Vertices, G2.Vertices);
	G3.Edges = 0;
	G3.Directed = G1.Directed;

	Graph* Less = (G1.Vertices <= G2.Vertices ? &G1 : &G2);
	Graph* Greater = (G1.Vertices >= G2.Vertices ? &G1 : &G2);
	int difference = Greater->Vertices - Less->Vertices;
	for (int i = 0; i < difference; i++)
		Less->AddVertex();

	for (int i = 0; i < G3.Vertices; i++)
	{
		vector <int> row;
		for (int j = 0; j < G3.Vertices; j++)
		{
			int x = max(G1.AdjacencyMatrix[i][j], G2.AdjacencyMatrix[i][j]);
			row.push_back(x);
			G3.Edges += x;
		}
		G3.AdjacencyMatrix.push_back(row);
	}

	if (!G3.Directed)
		G3.Edges /= 2;

	return G3;
}

Graph Graph::DisjunctDisjunction(Graph G1, Graph G2)
{
	Graph G3;
	int n1 = G1.Vertices, n2 = G2.Vertices;
	for (int i = 0; i < n1; i++)
	{
		G3.AdjacencyMatrix.push_back(G1.AdjacencyMatrix[i]);
		vector <int> row(n2, 0);
		G3.AdjacencyMatrix[i].insert(G3.AdjacencyMatrix[i].end(), row.begin(), row.end());
	}
	for (int i = n1; i < n1 + n2; i++)
	{
		vector <int> row(n1, 0);
		G3.AdjacencyMatrix.push_back(row);
		G3.AdjacencyMatrix[i].insert(G3.AdjacencyMatrix[i].end(), G2.AdjacencyMatrix[i - n1].begin(), G2.AdjacencyMatrix[i - n1].end());
	}
	G3.Vertices = n1 + n2;
	G3.Edges = G1.Edges + G2.Edges;
	G3.Directed = G1.Directed;
	return G3;
}

Graph Graph::Addition(Graph G1, Graph G2)
{
	Graph G3;
	int n1 = G1.Vertices, n2 = G2.Vertices;
	for (int i = 0; i < n1; i++)
	{
		G3.AdjacencyMatrix.push_back(G1.AdjacencyMatrix[i]);
		vector <int> row(n2, 1);
		G3.AdjacencyMatrix[i].insert(G3.AdjacencyMatrix[i].end(), row.begin(), row.end());
	}
	for (int i = n1; i < n1 + n2; i++)
	{
		vector <int> row(n1, 1);
		G3.AdjacencyMatrix.push_back(row);
		G3.AdjacencyMatrix[i].insert(G3.AdjacencyMatrix[i].end(), G2.AdjacencyMatrix[i - n1].begin(), G2.AdjacencyMatrix[i - n1].end());
	}
	G3.Vertices = n1 + n2;
	G3.Edges = G1.Edges + G2.Edges + n1 * n2;
	G3.Directed = G1.Directed;
	return G3;
}

Graph Graph::Multiplication(Graph G1, Graph G2)
{
	Graph G3;
	int n1 = G1.Vertices, n2 = G2.Vertices;
	G3.Vertices = n1 * n2;
	G3.Edges = G1.Edges * n2 + G2.Edges * n1;
	G3.Directed = G1.Directed;

	for (int a1 = 0; a1 < n1; a1++)
	{
		for (int b1 = 0; b1 < n2; b1++)
		{
			vector <int> row;
			for (int a2 = 0; a2 < n1; a2++)
			{
				for (int b2 = 0; b2 < n2; b2++)
				{
					if (((a1 == a2) && (G2.AdjacencyMatrix[b1][b2] > 0)) || ((b1 == b2) && (G1.AdjacencyMatrix[a1][a2] > 0)))
						row.push_back(1);
					else
						row.push_back(0);
				}
			}
			G3.AdjacencyMatrix.push_back(row);
		}
	}

	return G3;
}

void Graph::VertexDoubling(int n)
{
	if (n > Vertices)
		return;
	AddVertex();
	for (int j = 0; j < Vertices - 1; j++)
	{
		bool eq = (j == n - 1);
		for (int k = 0; k < AdjacencyMatrix[n - 1][j] >> eq; k++)
			AddEdge(Vertices, (eq ? Vertices : j + 1));
	}
}

void Graph::VertexReproduction(int n)
{
	if (n > Vertices)
		return;
	VertexDoubling(n);
	AddEdge(Vertices, n);
}

void Graph::EdgeSubdivision(int i, int j)
{
	if (i > Vertices || j > Vertices)
		return;
	AddVertex();
	RemoveEdge(i, j);
	AddEdge(i, Vertices);
	AddEdge(Vertices, j);
}

void Graph::VerticesIdentification(int n, int m)
{
	if (n > Vertices || m > Vertices)
		return;
	AddVertex();
	for (int i = 0; i < Vertices - 1; i++)
		AdjacencyMatrix[i][Vertices - 1] = max(AdjacencyMatrix[i][n - 1], AdjacencyMatrix[i][m - 1]);
	for (int j = 0; j < Vertices - 1; j++)
		AdjacencyMatrix[Vertices - 1][j] = max(AdjacencyMatrix[n - 1][j], AdjacencyMatrix[m - 1][j]);
	RemoveVertex(max(n, m));
	RemoveVertex(min(n, m));
}

void Graph::PrintSpanningTreeUsingDFS()
{
	Graph O;
	O.Vertices = Vertices;
	O.Edges = O.Vertices - 1;
	O.Directed = Directed;

	vector <int> row(O.Vertices, 0);
	O.AdjacencyMatrix.assign(O.Vertices, row);

	vector <bool> visited(Vertices, false);

	stack <pair<int, int>> S;
	S.push(make_pair(-1, 0));

	while (!S.empty())
	{
		int u0 = S.top().first;
		int u1 = S.top().second;
		S.pop();
		if (visited[u1])
			continue;
		visited[u1] = true;
		if (u0 != -1)
		{
			O.AdjacencyMatrix[u0][u1]++;
			if (!Directed)
				O.AdjacencyMatrix[u1][u0]++;
		}
		for (int v = Vertices - 1; v >= 0; v--)
		{
			if (AdjacencyMatrix[u1][v] <= 0)
				continue;
			if (!visited[v])
				S.push(make_pair(u1, v));
		}
	}

	O.PrintAdjacencyMatrix();
}

void Graph::PrintSpanningTreeUsingBFS()
{
	Graph O;
	O.Vertices = Vertices;
	O.Edges = O.Vertices - 1;
	O.Directed = Directed;

	vector <int> row(O.Vertices, 0);
	O.AdjacencyMatrix.assign(O.Vertices, row);

	vector <bool> visited(Vertices, false);

	queue <pair<int, int>> Q;
	Q.push(make_pair(-1, 0));

	while (!Q.empty())
	{
		int u0 = Q.front().first;
		int u1 = Q.front().second;
		Q.pop();

		if (visited[u1])
			continue;
		visited[u1] = true;

		if (u0 != -1)
		{
			O.AdjacencyMatrix[u0][u1]++;
			if (!Directed)
				O.AdjacencyMatrix[u1][u0]++;
		}

		for (int v = 0; v < Vertices; v++)
		{
			if (AdjacencyMatrix[u1][v] <= 0)
				continue;
			if (!visited[v])
				Q.push(make_pair(u1, v));
		}
	}

	O.PrintAdjacencyMatrix();
}

void Graph::PrintMinimumWeightedSpanningTree()
{
	Graph MinSpanningTree;
	MinSpanningTree.Vertices = Vertices;
	MinSpanningTree.Edges = MinSpanningTree.Vertices - 1;
	MinSpanningTree.Directed = Directed;

	vector <int> row(MinSpanningTree.Vertices, 0);
	MinSpanningTree.AdjacencyMatrix.assign(MinSpanningTree.Vertices, row);

	vector <pair<int, pair<int, int>>> g;
	for (int i = 0; i < Vertices; i++)
	{
		for (int j = (Directed ? 0 : i + 1); j < Vertices; j++)
		{
			if (AdjacencyMatrix[i][j] != 0)
				g.push_back(make_pair(AdjacencyMatrix[i][j], make_pair(i, j)));
		}
	}
	sort(g.begin(), g.end());

	int cost = 0;

	vector <int> tree_id(Vertices);
	for (int i = 0; i < Vertices; ++i)
		tree_id[i] = i;

	for (int i = 0; i < g.size(); ++i)
	{
		int a = g[i].second.first, b = g[i].second.second, l = g[i].first;
		if (tree_id[a] != tree_id[b])
		{
			cost += l;
			MinSpanningTree.AdjacencyMatrix[a][b]++;
			if (!Directed)
				MinSpanningTree.AdjacencyMatrix[b][a]++;
			int old_id = tree_id[b], new_id = tree_id[a];
			for (int j = 0; j < Vertices; ++j)
				if (tree_id[j] == old_id)
					tree_id[j] = new_id;
		}
	}

	cout << "Weight: " << cost << endl;
	MinSpanningTree.PrintAdjacencyMatrix();
}

void Graph::GetAllDistancesUsingFloyd()
{
	vector <vector <int>> D, P;
	D = AdjacencyMatrix;
	for (int i = 0; i < Vertices; i++)
	{
		vector <int> row;
		for (int j = 0; j < Vertices; j++)
		{
			if (D[i][j] != inf && D[i][j] != 0)
				row.push_back(i + 1);
			else
				row.push_back(0);
		}
		P.push_back(row);
	}

	for (int k = 0; k < Vertices; k++)
	{
		for (int i = 0; i < Vertices; i++)
		{
			for (int j = 0; j < Vertices; j++)
			{
				if (D[i][j] > D[i][k] + D[k][j])
				{
					D[i][j] = D[i][k] + D[k][j];
					P[i][j] = k + 1;
				}
			}
		}
	}
	cout << "Distances:" << endl;
	for (int i = 0; i < Vertices; i++)
	{
		for (int j = 0; j < Vertices; j++)
		{
			cout << D[i][j] << " ";
		}
		cout << endl;
	}
	cout << "Parents:" << endl;
	for (int i = 0; i < Vertices; i++)
	{
		for (int j = 0; j < Vertices; j++)
		{
			cout << P[i][j] << " ";
		}
		cout << endl;
	}
}

void Graph::GetDistancesInContourlessNetwork(int v1)
{
	vector <int> D(Vertices), Previous(Vertices), sorted;
	D[v1] = 0;
	Previous[v1] = 0;

	sorted = top_sort();

	for (int v = 0; v < Vertices; v++)
	{
		if (v == v1)
			continue;
		D[v] = inf;
	}
	for (int v = 0; v < Vertices; v++)
	{
		for (int w = 0; w < Vertices; w++)
		{
			if (AdjacencyMatrix[sorted[w]][sorted[v]] == inf || w == v)
				continue;
			if (D[sorted[w]] + AdjacencyMatrix[sorted[w]][sorted[v]] < D[sorted[v]])
			{
				D[sorted[v]] = D[sorted[w]] + AdjacencyMatrix[sorted[w]][sorted[v]];
				Previous[sorted[v]] = sorted[w] + 1;
			}
		}
	}

	cout << "Distances:" << endl;
	for (int i = 0; i < Vertices; i++)
	{
		if (D[i] == inf || D[i] == inf - 1)
			cout << "inf ";
		else
			cout << D[i] << " ";
	}
	cout << endl << "Parents:" << endl;
	for (int i = 0; i < Vertices; i++)
	{
		cout << Previous[i] << " ";
	}
	cout << endl;
}

void Graph::FindMaxFlow(int s, int t)
{
	vector <int> f0(Vertices, 0);
	vector <vector <int>> f(Vertices, f0);

	vector <int> d;
	int maxFlow = 0;
	vector <int> p;
	while (BFSForMaximalFlow(s, t, f, d))
	{
		p.assign(Vertices, 0);
		int flow = DFSForMaximalFlow(s, inf, s, t, p, f, d);
		while (flow != 0)
		{
			maxFlow += flow;
			flow = DFSForMaximalFlow(s, inf, s, t, p, f, d);
		}
	}
	for (int i = 0; i < Vertices; i++)
	{
		for (int j = 0; j < Vertices; j++)
		{
			cout << f[i][j] << " ";
		}
		cout << endl;
	}
	cout << "Maximal flow: " << maxFlow << endl;
}

int Graph::GetDistance(int n, int m)
{
	vector <int> distances = BFS(n - 1, m - 1, false);
	return distances[m - 1];
}

int Graph::GetEccentricity(int n)
{
	vector <int> distances = BFS(n - 1, -1, true);
	int max = 0;
	for (int v : distances)
		if (v > max)
			max = v;
	return max;
}

int Graph::GetRadius()
{
	int min = int_max;
	vector <int> eccentricities;
	for (int i = 1; i <= Vertices; i++)
		eccentricities.push_back(GetEccentricity(i));
	for (int e : eccentricities)
		if (e < min)
			min = e;
	return min;
}

int Graph::GetDiameter()
{
	int max = 0;
	vector <int> eccentricities;
	for (int i = 1; i <= Vertices; i++)
		eccentricities.push_back(GetEccentricity(i));
	for (int e : eccentricities)
		if (e > max)
			max = e;
	return max;
}

void Graph::PrintPeripheralVertices()
{
	int diameter = GetDiameter();
	cout << "Peripheral:" << endl << "{ ";
	for (int i = 1; i <= Vertices; i++)
		if (GetEccentricity(i) == diameter)
			cout << i << " ";
	cout << "}" << endl;
}

void Graph::PrintCentralVertices()
{
	int radius = GetRadius();
	cout << "Central:" << endl << "{ ";
	for (int i = 1; i <= Vertices; i++)
		if (GetEccentricity(i) == radius)
			cout << i << " ";
	cout << "}" << endl;
}

int Graph::GetDegree(int n)
{
	if (Directed)
		return -1;
	if (n > Vertices)
		return -1;
	int degree = 0;
	for (int j = 0; j < Vertices; j++)
	{
		degree += AdjacencyMatrix[n - 1][j] / (n - 1 == j ? 2 : 1);
	}
	return degree;
}

int Graph::GetInSemidegree(int n)
{
	if (!Directed)
		return -1;
	if (n > Vertices)
		return -1;
	int degree = 0;
	for (int i = 0; i < Vertices; i++)
	{
		degree += AdjacencyMatrix[i][n - 1];
	}
	return degree;
}

int Graph::GetOutSemidegree(int n)
{
	if (!Directed)
		return -1;
	if (n > Vertices)
		return -1;
	int degree = 0;
	for (int j = 0; j < Vertices; j++)
	{
		degree += AdjacencyMatrix[n - 1][j];
	}
	return degree;
}

int Graph::GetNumberOfVertices()
{
	return Vertices;
}

int Graph::GetNumberOfEdges()
{
	return Edges;
}

bool Graph::GetOrientation()
{
	return Directed;
}