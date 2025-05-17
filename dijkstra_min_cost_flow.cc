#include <cassert>
#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>

using namespace std;

const int INF = numeric_limits<int>::max();

struct Edge {
    /** @brief The destination vertex of the edge. */
    int to;
    /** @brief The maximum capacity of the edge. */
    int capacity;
    /** @brief The cost per unit of flow on this edge. */
    int cost;
    /** @brief The current flow on this edge. */
    int flow;
    /** @brief The index of the reverse edge in the adjacency list of the destination vertex. */
    int rev;

};

class MinCostFlow {

public:
    /** @brief Number of nodes in the network. */
    int n;
    /** @brief Adjacency list to represent the graph. */
    vector<vector<Edge>> adj;
    /** @brief Potential of each node for Dijkstra's algorithm. */
    vector<int> potential;
    /** @brief Distance to each node from the source in Dijkstra's. */
    vector<int> dist;
    /** @brief Previous node in shortest path from Dijkstra's. */
    vector<int> prevv, preve;

    MinCostFlow(int n) : n(n), adj(n), potential(n, 0), dist(n, INF), prevv(n), preve(n) {}

    /**
     * @brief Adds a directed edge to the flow network.
     * 
     * @param u The source vertex of the edge.
     * @param v The destination vertex of the edge.
     * @param capacity The capacity of the edge.
     * @param cost The cost per unit of flow through the edge.
     * 
     * Note that this also adds a reverse edge with a cost of -cost and a capacity of 0.
     */
    void add_edge(int u, int v, int capacity, int cost) {
        adj[u].push_back({v, capacity, cost, 0, (int)adj[v].size()});
        adj[v].push_back({u, 0, -cost, 0, (int)adj[u].size() - 1});
    }

    /**
     * @brief Finds the shortest path from the source to the sink using Dijkstra's algorithm with potential.
     * 
     * This function is a helper to find the shortest path from the source `s` to the sink `t` in the residual graph.
     * It uses Dijkstra's algorithm with a potential function to handle negative edge weights.
     * 
     * @param s The source vertex.
     * @param t The sink vertex.
     * @return true if a path is found, false otherwise.
     */
    bool dijkstra(int s, int t) {
        dist.assign(n, INF);
        prevv.assign(n, -1);
        preve.assign(n, -1);
        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;

        dist[s] = 0;
        pq.push({0, s});

        while (!pq.empty()) {
            pair<int, int> current = pq.top();
            pq.pop();
            int u = current.second;
            if (dist[u] < current.first) continue;

            for (int i = 0; i < adj[u].size(); ++i) {
                Edge& e = adj[u][i];
                if (e.capacity - e.flow > 0 && dist[e.to] > dist[u] + e.cost + potential[u] - potential[e.to]) {
                    dist[e.to] = dist[u] + e.cost + potential[u] - potential[e.to];
                    prevv[e.to] = u;
                    preve[e.to] = i;
                    pq.push({dist[e.to], e.to});
                }
            }
        }

        if (dist[t] == INF) return false;
        for (int v = 0; v < n; ++v) {
            if (dist[v] != INF) {
                potential[v] += dist[v];
            }
        }
        return true;
    }

    /**
     * @brief Computes the minimum cost maximum flow from source 's' to sink 't' with a maximum flow of 'f'.
     * 
     * This function iteratively finds the shortest path from the source 's' to the sink 't' in the residual graph using the dijkstra helper function.
     * It augments the flow along this path until the required flow 'f' is reached or no more path exist.
     * The algorithm uses the concept of potential to handle the case of negative cost edges.
     * 
     * Time Complexity: O(f * E * log(V)), where f is the maximum flow, E is the number of edges, and V is the number of vertices.
     *   - dijkstra is called f times (in the worst case).
     *   - dijkstra complexity is O(E * log(V)) using a priority queue.
     * Space Complexity: O(V + E) to store the graph and the supporting structures.
     * @param s The source vertex.
     * @param t The sink vertex.
     * @param f The maximum flow required.
     * @return A pair containing the achieved flow and the minimum cost.
     */
    pair<int, int> min_cost_flow(int s, int t, int f) {
        int flow = 0;
        int cost = 0;
        potential.assign(n, 0);

        while (f > 0 && dijkstra(s, t)) {
            int d = f;
            for (int v = t; v != s; v = prevv[v]) {
                d = min(d, adj[prevv[v]][preve[v]].capacity - adj[prevv[v]][preve[v]].flow);
            }
            if (d == 0) break;
            f -= d;
            flow += d;
            for (int v = t; v != s; v = prevv[v]) {
                Edge& e = adj[prevv[v]][preve[v]];
                e.flow += d;
                adj[v][e.rev].flow -= d;
                cost += d * e.cost;
            }
        }
        return {flow, cost};
    }
};

void test_min_cost_flow() {
    // Test case 1: Basic test with positive costs
    MinCostFlow mcf(4);
    mcf.add_edge(0, 1, 10, 2);
    mcf.add_edge(0, 2, 2, 4);
    mcf.add_edge(1, 3, 7, 6);
    mcf.add_edge(2, 3, 5, 1);
    pair<int, int> result = mcf.min_cost_flow(0, 3, 5);
    assert(result.first == 5);
    assert(result.second == 34);

    // Test case 2: Another basic test with different network structure
    MinCostFlow mcf2(6);
    mcf2.add_edge(0, 1, 10, 1);
    mcf2.add_edge(0, 2, 2, 2);
    mcf2.add_edge(1, 3, 6, 3);
    mcf2.add_edge(2, 4, 6, 4);
    mcf2.add_edge(3, 5, 10, 1);
    mcf2.add_edge(4, 5, 10, 1);
    pair<int, int> result_test2 = mcf2.min_cost_flow(0, 5, 5);
    assert(result_test2.first == 5);
    assert(result_test2.second == 25);

    // Test case 3: Test with INF flow, to test max flow computation
    MinCostFlow mcf3(6);
    mcf3.add_edge(0, 1, 2, 2);
    mcf3.add_edge(0, 2, 1, 3);
    mcf3.add_edge(1, 3, 2, 1);
    mcf3.add_edge(2, 3, 1, 2);
    mcf3.add_edge(3, 4, 3, 2);
    mcf3.add_edge(4, 5, 3, 1);
    pair<int, int> result_test3 = mcf3.min_cost_flow(0, 5, INF);
    assert(result_test3.first == 3);
    assert(result_test3.second == 20);

    // Test case 4: Test with negative cost edges
    MinCostFlow mcf4(4);
    mcf4.add_edge(0, 1, 3, -2);
    mcf4.add_edge(0, 2, 2, 4);
    mcf4.add_edge(1, 3, 2, 3);
    mcf4.add_edge(2, 3, 2, 1);
    pair<int, int> result_test4 = mcf4.min_cost_flow(0, 3, 3);
    assert(result_test4.first == 3);
    assert(result_test4.second == 7);

}

void run_min_cost_flow_sample(){
    MinCostFlow mcf(4);
    mcf.add_edge(0, 1, 10, 2);
    mcf.add_edge(0, 2, 2, 4);
    mcf.add_edge(1, 3, 7, 6);
    mcf.add_edge(2, 3, 5, 1);
    pair<int, int> result = mcf.min_cost_flow(0, 3, 5);
    cout << "Flow: " << result.first << ", Cost: " << result.second << endl;
}

int main() {
    test_min_cost_flow();
    run_min_cost_flow_sample();
    return 0;
}

