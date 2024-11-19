#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <ctime>
#include <random>
#include <queue>
#include <chrono>
#include <bits/stdc++.h>

using namespace std;

// GLOBAL VARIABLES
long long numNodes = 0, numEdges = 0;
vector<vector<long long>> matrix;
vector<vector<pair<long long, long long>>> graph;
vector<vector<long long>> vector_matrix;
vector<long long> d;
vector<bool> visited;

// FUNCTIONS

// LOAD FROM FILE
void LoadFile(string filename)
{
    matrix.clear();
    graph.clear();
    vector_matrix.clear();
    ifstream file(filename);
    if (!file)
    {
        cout << "File not found." << endl;
    }
    file >> numNodes >> numEdges;
    cout << endl
         << "Number of Nodes: " << numNodes << " Number of Edges: " << numEdges << endl;
    // cout << "NoE | Start | End | Weight" << endl;
    vector<long long> row;
    for (long long i = 0; i < numEdges + 1; i++)
    {
        long long node1, node2, weight; // inserting to matrix vector
        file >> node1 >> node2 >> weight;
        row.push_back(i);
        row.push_back(node1);
        row.push_back(node2);
        row.push_back(weight);
        matrix.push_back(row);
        row.clear();
    }

    // EDGES DISPLAY
    // for (long long i = 0; i < numEdges; i++)
    // {
    //     cout.width(3);
    //     cout << matrix[i][0] << " | " << matrix[i][1] << " | " << matrix[i][2] << " | " << matrix[i][3] << endl;
    // }
    // END EDGES DISPLAY

    // GRAPH LOAD
    graph.resize(numNodes);
    for (long long i = 0; i < numEdges + 1; ++i)
    {
        long long start = matrix[i][1];
        long long end = matrix[i][2];
        long long weight = matrix[i][3];
        bool updated = false;
        // Sprawdzanie i aktualizowanie istniejących krawędzi
        for (auto &edge : graph[start])
        {
            if (edge.second == end)
            {
                if (weight < edge.first)
                {
                    edge.first = weight; // Aktualizowanie na mniejszą wagę
                }
                updated = true;
                break;
            }
        }
        if (!updated)
        {
            graph[start].push_back(make_pair(weight, end));
        }
    }

    // END GRAPH LOAD
    // VECTOR MATRIX LOAD
    vector_matrix.resize(numNodes + 1, vector<long long>(numNodes, -1));
    for (long long i = 0; i < numNodes; i++)
    {
        for (long long j = 0; j < graph[i].size(); j++)
        {
            long long start = i;
            long long end = graph[i][j].second;
            long long weight = graph[i][j].first;

            if (weight != -1 && (vector_matrix[start][end] == -1 || weight < vector_matrix[start][end]))
            {
                vector_matrix[start][end] = weight;
            }
        }
    }
    // END VECTOR MATRIX LOAD

    file.close();
}
// END LOAD FROM FILE

unsigned long long GenerateRandom(long long max_val)
{
    unsigned long long random;
    random = rand() % max_val;

    return random;
}

void GenerateGraph(long long Density, long long numNodes, long long numEdges)
{
    string generated = "generated.txt";
    ofstream file(generated);
    if (!file)
    {
        cout << "File not found." << endl;
    }
    file << numNodes << " " << numEdges << endl;
    // Generating Edges
    for (long long i = 0; i < numEdges; i++)
    {
        vector<long long> row;
        long long node1, node2, weight; // inserting to matrix vector

        node1 = GenerateRandom(numNodes);
        node2 = GenerateRandom(numNodes);
        while (node2 == node1)
        {
            node2 = GenerateRandom(numNodes);
        }
        weight = GenerateRandom(1000000);
        file << node1 << " " << node2 << " " << weight << endl;
        row.push_back(i);
        row.push_back(node1);
        row.push_back(node2);
        row.push_back(weight);
        matrix.push_back(row);
        row.clear();
    }
    // LOAD FROM FILE
    LoadFile(generated);
    // END LOAD FROM FILE
}

// DIJKSTRA'S ALGORITHM ON LSIT

void Dijkstra(vector<long long> &d, vector<bool> &visited, long long start, long long end)
{
    d.assign(numNodes + 1, numeric_limits<long long>::max());
    visited.assign(numNodes + 1, false);
    priority_queue<pair<long long, long long>, vector<pair<long long, long long>>, greater<pair<long long, long long>>> train;

    long long dest = 0, here = 0, waga = 0;

    d[start] = 0;
    train.push({0, start});
    while (!train.empty())
    {
        here = train.top().second; // second is destination, first is weight
        train.pop();
        if (!visited[here])
        {
            visited[here] = true;
            for (long long i = 0; i < graph[here].size(); i++)
            {
                waga = graph[here][i].first;
                dest = graph[here][i].second;
                if (d[here] + waga < d[dest])
                {
                    d[dest] = d[here] + waga;
                }
                train.push({d[dest], dest});
            }
        }
    }
}

// END DIJKSTRA'S ALGORITHM ON LIST
// DIJKSTRA'S ALGORITHM ON MATRIX

void DijkstraMatrix(vector<long long> &d, vector<bool> &visited, long long start, long long end)
{
    d.assign(numNodes + 1, numeric_limits<long long>::max());
    visited.assign(numNodes + 1, false);
    priority_queue<pair<long long, long long>, vector<pair<long long, long long>>, greater<pair<long long, long long>>> train;

    long long dest = 0, here = 0, waga = 0;

    d[start] = 0;
    train.push({0, start});
    while (!train.empty())
    {
        here = train.top().second; // second is destination, first is weight
        train.pop();
        if (!visited[here])
        {
            visited[here] = true;
            for (long long i = 0; i < numNodes; i++)
            {
                if (vector_matrix[here][i] != -1)
                {
                    waga = vector_matrix[here][i];
                    dest = i;
                    if (d[here] + waga < d[dest])
                    {
                        d[dest] = d[here] + waga;
                    }
                    train.push({d[dest], dest});
                }
            }
        }
    }
}

// END DIJKSTRA'S ALGORITHM ON MATRIX

double doDijkstra()
{
    d.clear();
    visited.clear();
    long long start = 0, end = numNodes - 1;

    chrono::time_point<std::chrono::high_resolution_clock> start_clock, end_clock; // zmienne mierzonego czasu
    start_clock = chrono::high_resolution_clock::now();
    //Dijkstra(d, visited, start, end);
    DijkstraMatrix(d, visited, start, end);
    end_clock = chrono::high_resolution_clock::now();
    chrono::duration<double> result = end_clock - start_clock;
    return result.count();
}

int main()
{
    ofstream wynik, licznik;
    wynik.open("wynik.txt", ios::app);
    licznik.open("licznik.txt", ios::app);
    srand(time(0));
    long long choice = 1;
    while (choice != 0)
    {

        // MENU
        cout << "Testing Effectiveness of Graph Algorythms" << endl;
        cout << "----------------------------------------" << endl;
        cout << "1. Load from file" << endl;
        cout << "2. Generate (generated.txt)" << endl;
        cout << "3. Display (list/matrix)" << endl;
        cout << "4. Dijkstra's algorithm" << endl;
        cout << "5. Gathering results (long wait)" << endl;
        cout << "6. Quit" << endl;
        cout << "Enter your choice: ";
        cin >> choice;
        if (choice == 6)
        {
            return 0;
        }
        // END MENU

        // LOAD FROM FILE
        if (choice == 1)
        {
            string filename;
            cout << "Enter filename: ";
            cin >> filename;
            LoadFile(filename);
        }
        // GENERATE RANDOM GRAPH
        else if (choice == 2)
        {
            long long Density;
            cout << "Insert number of Nodes: " << endl;
            cin >> numNodes;
            cout << "Insert Density of this graph (in %): " << endl;
            cin >> Density;

            numEdges = ((Density * 0.01) * numNodes * (numNodes - 1));
            GenerateGraph(Density, numNodes, numEdges);
        }
        // END GENERATE RANDOM GRAPH

        // DISPLAY
        else if (choice == 3)
        {
            if (numNodes == 0)
            {
                cout << "No graph loaded." << endl;
                continue;
            }
            cout << "Displaying graph..." << endl;
            cout << "1 - as list" << endl;
            cout << "2 - as matrix" << endl;
            cout << "Enter choice: ";
            cin >> choice;
            if (choice == 1)
            {
                vector<long long> order;
                cout << "Graph as list:" << endl;
                for (long long i = 0; i < numNodes; i++)
                {
                    cout << "Node " << i << ": ";
                    for (long long j = 0; j < graph[i].size(); j++)
                    {
                        order.push_back(graph[i][j].second);
                    }
                    order.erase(unique(order.begin(), order.end()), order.end());
                    sort(order.begin(), order.end());
                    for (long long k = 0; k < order.size(); k++)
                    {
                        cout.width(4);
                        cout << order[k] << " ";
                    }
                    order.clear();
                    cout << endl;
                }
                cout << endl;
            }
            else if (choice == 2)
            {
                cout << "Graph as matrix:" << endl;
                for (long long i = 0; i < numNodes; i++)
                {
                    cout.width(7);
                    cout << "Node " << i << ": ";
                    for (long long j = 0; j < numNodes; j++)
                    {
                        cout.width(10);
                        cout << vector_matrix[i][j] << " | ";
                    }
                    cout << endl;
                }
                cout << endl;
            }
            else
            {
                cout << "Invalid choice" << endl;
            }
        }
        // END DISPLAY
        // DIJKSTRA'S ALGORITHM
        else if (choice == 4)
        {

            unsigned long long start = numNodes + 1, end = numNodes + 1;
            cout << "Dijkstra's algorithm on..." << endl;
            cout << "...list - 1" << endl;
            cout << "...matrix - 2" << endl;
            cin >> choice;
            if (choice == 1)
            {

                while (start > numNodes)
                {
                    cout << "Insert starting node: " << endl;
                    cin >> start;
                }
                cout << "Insert destination node: " << endl;
                cin >> end;

                if (end == start)
                {
                    cout << "The same node was chosen." << endl;
                    cout << "Distance is 0." << endl;
                }
                else
                {
                    d.clear();
                    visited.clear();

                    chrono::time_point<std::chrono::high_resolution_clock> start_clock, end_clock; // zmienne mierzonego czasu
                    start_clock = chrono::high_resolution_clock::now();
                    Dijkstra(d, visited, start, end);
                    end_clock = chrono::high_resolution_clock::now();
                    chrono::duration<double> result = end_clock - start_clock;
                    licznik << result.count() << " " << endl;
                }
            }
            else if (choice == 2)
            {
                while (start > numNodes)
                {
                    cout << "Insert starting node: " << endl;
                    cin >> start;
                }
                cout << "Insert destination node: " << endl;
                cin >> end;

                if (end == start)
                {
                    cout << "The same node was chosen." << endl;
                    cout << "Distance is 0." << endl;
                }
                else
                {
                    d.clear();
                    visited.clear();

                    chrono::time_point<std::chrono::high_resolution_clock> start_clock, end_clock; // zmienne mierzonego czasu
                    start_clock = chrono::high_resolution_clock::now();
                    DijkstraMatrix(d, visited, start, end);
                    end_clock = chrono::high_resolution_clock::now();
                    chrono::duration<double> result = end_clock - start_clock;
                    licznik << result.count() << " " << endl;
                }
            }
            // DISTANCE FROM START TO END
            if (d[end] == numeric_limits<long long>::max())
            {
                cout << "No path from " << start << " to " << end << endl;
            }
            else
            {
                for (int i = 0; i < numNodes; i++)
                {
                    if (i == start)
                    {
                        wynik << "0 ";
                        continue;
                    }
                    else if (d[i] == numeric_limits<long long>::max())
                    {
                        wynik << "-1 ";
                        cout << "No path from " << start << " to " << i << endl;
                    }
                    else
                    {
                        cout << "Distance from " << start << " to " << i << " is: " << d[i] << endl;
                    }

                    if (d[i] != numeric_limits<long long>::max())
                    {
                        wynik << d[i] << " ";
                        // wynik << start<< " do "<<i<<": "<<d[i] << " "<<endl;
                    }
                }
            }
            wynik << endl;
            cout << endl;
            // END DISTANCE FROM START TO END
        }

        else if (choice == 5) // TESTING MODE
        {
            long long Density;
            for (numNodes = 250; numNodes <= 250 * 6; numNodes += 250)
            {
                for (Density = 25; Density <= 100; Density += 25)
                {
                    for (int i = 0; i < 20; i++)
                    {
                        numEdges = ((Density * 0.01) * numNodes * (numNodes - 1));
                        GenerateGraph(Density, numNodes, numEdges);
                        licznik << doDijkstra() << " ";
                    }
                    licznik << endl;
                }
                licznik << endl;
            }
        }
    }
    return 0;
}
