#include <iostream>
#include <vector>
#include <queue>
#include <cassert>
#include <limits>
#include <algorithm>

using namespace std;

/**
 * @brief Computes the shortest distances from a starting node to all other nodes in a graph using Dijkstra's algorithm.
 * @note Time Complexity: O(E log V), where E is the number of edges and V is the number of vertices.
 * @note Space Complexity: O(V + E), where V is the number of vertices and E is the number of edges.
 * 
 * @param start The index of the starting node.
 * @param adj The adjacency list representing the graph, where each element is a vector of pairs.
 *            Each pair (v, w) in adj[u] represents an edge from u to v with weight w.
 * @return A vector of shortest distances from the starting node to each other node.
 */
vector<int> dijkstra(int start, const vector<vector<pair<int, int>>>& adj) {
    int n = adj.size();
    // Initialize distances to infinity
    vector<int> dist(n, numeric_limits<int>::max());
    // Min-heap to store nodes to visit, ordered by distance
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;

    dist[start] = 0;
    pq.push({0, start});

    while (!pq.empty()) {
        int d = pq.top().first;
        int u = pq.top().second;
        pq.pop();

        if (d > dist[u]) continue;

        for (const auto& edge : adj[u]) {
            int v = edge.first;
            int weight = edge.second;
            if (dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                pq.push({dist[v], v});
            }
        }
    }

    return dist;
}

void testDijkstra() {
    // Test case 1: Basic graph
    vector<vector<pair<int, int>>> adj1 = {
        {{1, 4}, {2, 2}},
        {{3, 5}},
        {{1, 1}, {3, 8}},
        {}
    };
    vector<int> expected1 = {0, 3, 2, 8};
    vector<int> result1 = dijkstra(0, adj1);
    assert(result1 == expected1);

    // Test case 2: Disconnected graph
    vector<vector<pair<int, int>>> adj2 = {
        {{1, 1}},
        {},
        {{3, 1}},
        {}
    };
    vector<int> expected2 = {0,1, numeric_limits<int>::max(),numeric_limits<int>::max()};
    vector<int> result2 = dijkstra(0, adj2);
    assert(result2 == expected2);

    // Test case 3: Non-zero start node
    vector<vector<pair<int, int>>> adj3 = {
        {{1, 4}, {2, 2}},
        {{3, 5}},
        {{1, 1}, {3, 1}},
        {{4, 2}},
        {{0, 1}}
    };
    vector<int> expected3 = {4, 1, 0, 1, 3};
    vector<int> result3 = dijkstra(2, adj3);
    assert(result3 == expected3);

}

void runDijkstraSample() {
    // Sample graph
    vector<vector<pair<int, int>>> adj = {
        {{1, 4}, {2, 2}},
        {{3, 5}},
        {{1, 1}, {3, 8}},
        {}
    };
    int startNode = 0;
    
    // Run Dijkstra's algorithm
    vector<int> distances = dijkstra(startNode, adj);
    
    // Output the results
    cout << "Shortest distances from node " << startNode << ":" << endl;
    for (size_t i = 0; i < distances.size(); ++i) {
        cout << "To node " << i << ": ";
        if (distances[i] == numeric_limits<int>::max()) {
            cout << "unreachable";
        } else {
            cout << distances[i];
        }
        cout << endl;
    }
}

int main() {
    testDijkstra();
    runDijkstraSample();
    return 0;
}
