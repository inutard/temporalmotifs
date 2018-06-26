#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <map>
#include <random>
#include <cassert>

using namespace std;

struct full_edge {
	int src, dst, t;
	full_edge(int src, int dst, int t) : src(src), dst(dst), t(t) {};
	bool operator<(const full_edge& o) const {
		if (o.t != t) return t < o.t;
		return make_pair(src, dst) < make_pair(o.src, o.dst);
	}
};

struct half_edge {
	int dst, t;
};

typedef map<int, vector<half_edge>> Graph;

double count_2tmotifs(const vector<pair<int, int>>& subgraph, int delta, int window) {	
	// lets count motifs with 3 fwd edges
	//dp[num fwd][curr edge]
	double res = 0;
	for (int i = 0; i < subgraph.size(); i++) {
		if (subgraph[i].second == 0) continue;

		int nfwd = 1, nbk = 0;
		for (int j = i+1; j < subgraph.size(); j++) {
			if (subgraph[j].first - subgraph[i].first > delta) {
				break;
			}
			nfwd += subgraph[j].second;
			nbk += !subgraph[j].second;

			if (subgraph[j].second == 1) {
				int dt = subgraph[j].first - subgraph[i].first;
				res += window * nbk / double(window - dt);
			}
		}
	}

	for (int i = 0; i < subgraph.size(); i++) {
		if (subgraph[i].second == 1) continue;

		int nfwd = 1, nbk = 0;
		for (int j = i+1; j < subgraph.size(); j++) {
			if (subgraph[j].first - subgraph[i].first > delta) {
				break;
			}
			nfwd += !subgraph[j].second;
			nbk += subgraph[j].second;

			if (subgraph[j].second == 0) {
				int dt = subgraph[j].first - subgraph[i].first;
				res += window * nbk / double(window - dt);
			}
		}
	}
	return res;
}

// counts the weighted tmotifs, weighted by 1/(tf - ti).
double count_tmotifs(Graph& g, int delta, int window) {
	// counts 2 node motifs
	// pool up all the edges as (u, v, t) with u < v
	map<pair<int, int>, vector<pair<int, int>>> subgraphs;
	for (auto& kv : g) {
		int src = kv.first;
		for (auto e : kv.second) {
			int dst = e.dst, t = e.t;
			if (src < dst) {
				subgraphs[make_pair(src, dst)].push_back(make_pair(t, 1));
			} else {
				subgraphs[make_pair(dst, src)].push_back(make_pair(t, 0));
			}
		}
	}

	double res = 0;
	for (auto& subgraphKV : subgraphs) {
		sort(subgraphKV.second.begin(), subgraphKV.second.end());
		res += count_2tmotifs(subgraphKV.second, delta, window);
	}
	return res;
}

int main(int argc, char* argv[]) {
	ios::sync_with_stdio(0);
	cin.tie(0);

	freopen(argv[1], "r", stdin);
	int delta = atoi(argv[2]);

	vector<full_edge> edges;
	int u, v, t;
	while (cin >> u >> v >> t) {
		if (u == v) continue;
		edges.push_back({u, v, t});
	}
	sort(edges.begin(), edges.end());
	cerr << "Number of edges: " << edges.size() << endl;

	default_random_engine generator;
	uniform_real_distribution<double> distribution(0.0, 1.0);

	Graph g;
	srand(time(0));
	double st = clock();

	const int totaltime = edges.back().t - edges[0].t;
	const int window = 2 * delta;
	int n_innertrials = 24000, n_outertrials = 10;
	
	int offset = rand() % window;
	double tot_estimate = 0, tot_var = 0;
	for (int c_outer = 0; c_outer < n_outertrials; c_outer++) {
		double inner_est = 0;
		for (int c_inner = 0; c_inner < n_innertrials; c_inner++) {
			// pick a random delta offset
			int npos = 0;
			// pick a random edge and a random delta window around it
			int ei = rand() % edges.size();

			int ti = edges[ei].t;
			int st = ((ti-offset)/window) * window + offset + window;

			int ne = 0;
			int idx = lower_bound(edges.begin(), edges.end(), full_edge(-1, -1, st - window)) - edges.begin();
			while (idx < edges.size() && edges[idx].t <= st) {
				auto e = edges[idx];
				g[e.src].push_back({e.dst, e.t});
				idx++;
				ne++;
			}

			double pi = ne / double(edges.size());
			double res = count_tmotifs(g, delta, window) / pi;
			if (res > 0) npos++;
			inner_est += res;
			g.clear();
		}
		double est = double(inner_est) / n_innertrials;
		tot_var += est * est;
		tot_estimate += est;
		cerr << "Done trial " << c_outer << endl;
		cerr << "Estimate: " << est << endl;
	}

	double tvar = 0, tans = 0;
	int nedges = 0, nxt = edges[0].t + offset, npos = 0;
	vector<full_edge> reserve;
	for (const auto& e : edges) {
		int src, dst, t;
		tie(src, dst, t) = make_tuple(e.src, e.dst, e.t);
		if (nxt <= t) {
			double pi = double(nedges) / edges.size();
			double p = distribution(generator);

			for (const auto& e : reserve) {
				g[e.src].push_back({e.dst, e.t});
			}

			double tres = count_tmotifs(g, delta, window);
			tvar += tres * tres / pi;
			tans += tres;

			g.clear();
			reserve.clear();
			nedges = 0;
			nxt = ((t-offset)/window) * window + window + offset;
		}
		reserve.push_back({src, dst, t});
		nedges++;
	}

	double pi = double(nedges) / edges.size();
	double p = distribution(generator);

	for (const auto& e : reserve) {
		g[e.src].push_back({e.dst, e.t});
	}

	double tres = count_tmotifs(g, delta, window);
	tvar += tres * tres / pi;
	tans += tres;

	cerr << "Time (s): " << (clock() - st) / CLOCKS_PER_SEC << endl;

	double mean = tot_estimate / n_outertrials;
	double var = (tot_var / n_outertrials - mean * mean);
	cerr << "(Mean, var): " << mean << " " << var << endl;

	tvar -= tans * tans;
	cerr << "Theoretical ans: " << tans << endl;
	cerr << "Theoretical var: " << tvar / n_innertrials << endl;
	return 0;
}