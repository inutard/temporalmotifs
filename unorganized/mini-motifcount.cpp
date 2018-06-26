#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <map>
#include <random>

using namespace std;

#define STATS_ON 1

struct full_edge {
	int src, dst, t;
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

	const int window = 2 * delta;
	const int samplerate = 2;

	// F[i] = ci / pi with prob pi
	//		   0 otherwise
	// sum pi * F[i] = sum C[i]
	// variance is 
	// sum pi F[i]^2 - (sum C[i])^2
	// sum C[i]^2 / pi - (sum C[i])^2
#ifdef STATS_ON
	double tvar = 0;
	double tsum2 = 0;
	double avgp = 0;
#endif

	int ntrials = 20;
	int offset = rand() % window;
	double tot_estimate = 0, tot_var = 0;
	const double mult = 4e3;
	for (int c = 0; c < ntrials; c++) {
		cerr << "Trial " << c << endl;

#ifdef STATS_ON
		tvar = 0;
		tsum2 = 0;
#endif
		
		int nsamp = 0, ntot = 0, nedges = 0;

		double res = 0, totp = 0;
		int nxt = edges[0].t + offset, npos = 0;
		vector<full_edge> reserve;
		for (const auto& e : edges) {
			int src, dst, t;
			tie(src, dst, t) = make_tuple(e.src, e.dst, e.t);
			if (nxt <= t) {
				ntot++;

				double pi = min(mult * double(nedges) / edges.size(), 1.0);
				double p = distribution(generator);
				totp += pi;
				//cout << pi << " " << p << " " << nedges << endl;
				if (p <= pi) {
					for (const auto& e : reserve) {
						g[e.src].push_back({e.dst, e.t});
					}

					res += count_tmotifs(g, delta, window) / pi;
					if (res > 0) npos++;
					nsamp++;
				}
#ifdef STATS_ON
				if (p > pi) {
					for (const auto& e : reserve) {
						g[e.src].push_back({e.dst, e.t});
					}
				}

				double tres = count_tmotifs(g, delta, window);
				tvar += tres * tres / pi - tres * tres;
				tsum2 += tres * tres;
				/*
				if (tres) {
					cerr << nedges << endl;
					cerr << "Contrib: " << pi << " " << tres * tres / pi - tres * tres << " " << tres << endl;
				}
				*/
#endif

				g.clear();
				reserve.clear();
				nedges = 0;
				nxt = ((t-offset)/window) * window + window + offset;
			}
			reserve.push_back({src, dst, t});
			nedges++;
		}

		ntot++;
		double pi = min(mult * double(nedges) / edges.size(), 1.0);
		double p = distribution(generator);
		totp += pi;
		if (p <= pi) {
			ntot++;
			for (const auto& e : reserve) {
				g[e.src].push_back({e.dst, e.t});
			}

			res += count_tmotifs(g, delta, window) / pi;
			if (res > 0) npos++;
			nsamp++;
		}

#ifdef STATS_ON
		if (p > pi) {
			for (const auto& e : reserve) {
				g[e.src].push_back({e.dst, e.t});
			}
		}

		double tres = count_tmotifs(g, delta, window);
		tvar += tres * tres / pi - tres * tres;
		tsum2 += tres * tres;
		avgp = totp / ntot;
		//cerr << "Contrib: " << pi << " " << tres * tres / pi - tres * tres << endl;
#endif

		g.clear();
		reserve.clear();

		double est = res; //ntot * res / nsamp;
		tot_estimate += est;
		tot_var += est * est;
		cerr << est << endl;
		cerr << "Number sampled: " << nsamp << " " << ntot << endl;
		cerr << "Avg sample prob: " << totp / ntot << endl;
		cerr << "Total non-zero windows: " << ntot << endl;
		cerr << "Estimate: " << est << endl;
	}
	cerr << "Time (s): " << (clock() - st) / CLOCKS_PER_SEC << endl;

	double mean = tot_estimate / ntrials;
	cerr << "(Mean, var): " << mean << " " 
		 << (tot_var / ntrials - mean * mean) << endl;
		 //<< sqrt(tot_var / ntrials - mean * mean) << endl;

#ifdef STATS_ON
	cerr << "Theoretical var: " << tvar << endl;
	//cerr << (1/avgp - 1) << endl;
	cerr << "Sum of squared c_i: " << tsum2 << endl;
#endif
	return 0;
}