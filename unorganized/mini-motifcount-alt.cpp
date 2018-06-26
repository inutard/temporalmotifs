#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <map>
#include <random>

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

/*
struct MotifCounter {
	unordered_map<int, Graph> gt;
	MotifCounter(Graph g, int delta) {
	}

	int count_tmotifs() {

	}
};
*/

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
	const int window = delta;
	int ntrials = 10000;
	
	int offset = rand() % window;
	int nxt = edges[0].t + offset;
	int cnt = 0;
	int nt = 1;
	double mean = 0, var = 0;
	for (auto e : edges) {
		while (e.t >= nxt) {
			mean += cnt;
			var += cnt * cnt;
			cnt = 0;
			nxt += window;
			nt++;
		}
		cnt++;
	}
	mean += cnt;
	var += cnt * cnt;

	mean = mean / nt;
	var = var / nt - mean * mean;
	//cerr << mean << " " << var << endl;
	//double m_x = 47.504;
	double c_t = -0.005;//-4.2;

	double tot_estimate = 0, tot_var = 0, tot_edges = 0;
	double cov = 0;
	for (int c = 0; c < ntrials; c++) {
		int w_idx = rand() % nt;
		int st = w_idx * window + edges[0].t + offset;

		int idx = lower_bound(edges.begin(), edges.end(), full_edge(-1, -1, st - window)) - edges.begin();
		vector<full_edge> reserve;
		while (idx < edges.size() && edges[idx].t <= st) {
			reserve.push_back(edges[idx]);
			idx++;
		}

		if (reserve.size()) {
			for (const auto& e : reserve) {
				g[e.src].push_back({e.dst, e.t});
			}
			tot_estimate += count_tmotifs(g, delta, window);
		}
		//cov += (count_tmotifs(g, delta, window) - m_x) * (double(reserve.size()) - mean);
		tot_estimate += c_t * (double(reserve.size()) - mean);
		g.clear();
		reserve.clear();
	}
	//cov /= ntrials;
	//cout << cov / var << " " << cov << " " << var << endl;
	cerr << "Time (s): " << (clock() - st) / CLOCKS_PER_SEC << endl;

	//cerr << double(tot_estimate) / ntrials << endl;
	double est = double(tot_estimate) / ntrials * double(totaltime) / window;
	cerr << "Mean: " << est << endl;
	return 0;
}