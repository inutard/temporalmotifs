// temporalmotifsmain.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include "temporalmotifsfaster.h"
#ifdef USE_OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("Temporalmotifs. build: %s, %s. Time: %s",
			 __TIME__, __DATE__, TExeTm::GetCurTm()));  
  TExeTm ExeTm;  
  Try

  const TStr temporal_graph_filename =
    Env.GetIfArgPrefixStr("-i:", "example-temporal-graph.txt",
			  "Input directed temporal graph file");
  const TStr output = 
    Env.GetIfArgPrefixStr("-o:", "temporal-motif-counts.txt",
			  "Output file in which to write counts");
  const TFlt delta =
    Env.GetIfArgPrefixFlt("-delta:", 4096, "Time window delta");
  const int num_threads =
    Env.GetIfArgPrefixInt("-nt:", 4, "Number of threads for parallelization");

#ifdef USE_OPENMP
  omp_set_num_threads(num_threads);
#endif

  // Count all 2-node and 3-node temporal motifs with 3 temporal edges
  FTempMotifCounter tmc(temporal_graph_filename);
  Counter2D counts;
  tmc.Count3TEdge2NodeFaster(delta, counts);
  int c00 = counts(0,0);
  int c01 = counts(0,1);
  int c10 = counts(1,0);
  int c11 = counts(1,1);
  printf("%d %d %d %d\n", c00, c01, c10, c11);
  
  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(),
	 TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
