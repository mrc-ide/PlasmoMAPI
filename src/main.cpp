
#include "main.h"
#include "utils.h"
#include "probability_v10.h"
#include "sim.Parameters.h"
#include "sim.Dispatcher.h"

#include <chrono>
#include <vector>

using namespace std;


//------------------------------------------------
// assign edges to hexes based on intersection
Rcpp::List assign_map_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
  
  print("Assigning edges to hexes");
  
  // load data and parameters
  vector<double> node_long = rcpp_to_vector_double(args["node_long"]);    //Longitude of data nodes
  vector<double> node_lat = rcpp_to_vector_double(args["node_lat"]);      //Latitude of data nodes
  vector<double> hex_long = rcpp_to_vector_double(args["hex_long"]);      //Longitude of hexes
  vector<double> hex_lat = rcpp_to_vector_double(args["hex_lat"]);        //Latitude of hexes
  double hex_width = rcpp_to_double(args["hex_width"]);                   //Width of hexes
  double eccentricity = rcpp_to_double(args["eccentricity"]);             //Eccentricity of ellipses (see help for details)
  bool report_progress = rcpp_to_bool(args["report_progress"]);           //Whether to update progress bar
  bool pb_markdown = rcpp_to_bool(args["pb_markdown"]);                   //Whether to run in markdown-safe mode
  Rcpp::Function update_progress = args_functions["update_progress"];     //R function for updating progress bar
  
  // get basic properties
  int n_node = node_long.size();
  int n_hex = hex_long.size();
  
  // store list of which edges intersect each hex
  vector<vector<int>> hex_edges(n_hex);
  
  // loop through hexes
  for (int hex = 0; hex < n_hex; ++hex) {
    
    // report progress
    if (report_progress) {
      if ((hex+1) == n_hex) {
        update_progress(args_progress, "pb", hex+1, n_hex);
      } else {
        int remainder = hex % int(ceil(double(n_hex)/100));
        if (remainder == 0 && !pb_markdown) {
          update_progress(args_progress, "pb", hex+1, n_hex);
        }
      }
    }
    
    // loop through pairwise nodes
    int i = 0;
    for (int node1 = 0; node1 < (n_node-1); ++node1) {
      for (int node2 = (node1+1); node2 < n_node; ++node2) {
        i++;
        
        // determine whether ellipse intersects this hex
        bool intersects = collision_test_hex_ellipse(hex_long[hex], hex_lat[hex], hex_width,
                                                     node_long[node1], node_lat[node1],
                                                     node_long[node2], node_lat[node2],
                                                     eccentricity);
        
        // push back edge index if intersects
        if (intersects) {
          hex_edges[hex].push_back(i);
        }
      }
    }
  }
  
  // return list
  return Rcpp::List::create(Rcpp::Named("hex_edges") = hex_edges);
}


//------------------------------------------------
// assign edges to hexes based on intersection - alternate version
Rcpp::List assign_map2_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {

  int i, j, node1, node2, node_loc1, node_loc2, n_dups1, n_dups2;
  double node_long1, node_lat1, hex_long1, hex_lat1, dev;
  
  print("Assigning edges to hexes");
  
  // load data and parameters
  vector<double> node_long = rcpp_to_vector_double(args["node_long"]);    //Longitude of data nodes
  vector<double> node_lat = rcpp_to_vector_double(args["node_lat"]);      //Latitude of data nodes
  vector<double> hex_long = rcpp_to_vector_double(args["hex_long"]);      //Longitude of hexes
  vector<double> hex_lat = rcpp_to_vector_double(args["hex_lat"]);        //Latitude of hexes
  double hex_width = rcpp_to_double(args["hex_width"]);                   //Width of hexes
  double eccentricity = rcpp_to_double(args["eccentricity"]);             //Eccentricity of ellipses (see help for details)
  bool report_progress = rcpp_to_bool(args["report_progress"]);           //Whether to update progress bar
  bool pb_markdown = rcpp_to_bool(args["pb_markdown"]);                   //Whether to run in markdown-safe mode
  Rcpp::Function update_progress = args_functions["update_progress"];     //R function for updating progress bar
  
  // get basic properties
  int n_node = node_long.size();
  int n_hex = hex_long.size();

  // Check through nodes for ones with same locations
  int n_node_locs = 0;
  vector<int> flagDuplicate(n_node);
  vector<int> node_loc_list(n_node);
  vector<vector<int>> ellipse_list(n_node, vector<int>(n_node));
  int ellipse = 0;
  for(node1 = 0; node1 < (n_node-1); node1++)
  { 
    flagDuplicate[node1]=0; 
    for (node2 = (node1+1); node2 < n_node; node2++) 
    {
        ellipse++;
        ellipse_list[node1][node2]=ellipse;
        ellipse_list[node2][node1]=ellipse;
	}
  }
  
  flagDuplicate[n_node-1]=0;

  for (node1 = 0; node1 < (n_node-1); node1++) 
  {
    if(flagDuplicate[node1]==0)
    {
        node_loc_list[n_node_locs]=node1;
        node_long1=node_long[node1];
        node_lat1=node_lat[node1];
        for (node2 = (node1+1); node2 < n_node; node2++) 
        {
            dev=abs(sq(node_long1-node_long[node2])+sq(node_lat1-node_lat[node2]));
            if(dev==0.0)
            { // Duplicate found    
                flagDuplicate[node2]=1;
	    	}
	    }    
        n_node_locs++;
	}
  }

  //Rcpp::Rcout << "\nNumber of nodes = " << n_node << "\tNumber of non-duplicate nodes = " << n_node_locs << "\nNode location list:";
  //for(node1 = 0; node1 < n_node_locs; node1++){ Rcpp::Rcout << " " << node_loc_list[node1];}
  //R_FlushConsole();

  vector<vector<int>> DuplicateLists(n_node_locs);
  vector<int> nDuplicates(n_node_locs);
  vector<double> loc_long(n_node_locs);
  vector<double> loc_lat(n_node_locs);

  for(node_loc1 = 0; node_loc1 < n_node_locs; node_loc1++)
  {
      node1 = node_loc_list[node_loc1];
      DuplicateLists[node_loc1].push_back(node1);
      node_long1 = node_long[node1];
      node_lat1 = node_lat[node1];
      for(node2 = (node1+1); node2 < n_node; node2++)
      {
        dev=abs(sq(node_long1-node_long[node2])+sq(node_lat1-node_lat[node2]));
        if(dev==0.0)
        { // Duplicate found  
            DuplicateLists[node_loc1].push_back(node2);
	    }      
	  }
      nDuplicates[node_loc1]=DuplicateLists[node_loc1].size();
      loc_long[node_loc1]=node_long1;
      loc_lat[node_loc1]=node_lat1;
  }
    
  // store list of which edges intersect each hex
  vector<vector<int>> hex_edges(n_hex);
  
  // loop through hexes
  for (int hex = 0; hex < n_hex; ++hex) {
    
    // report progress
    if (report_progress) {
      if ((hex+1) == n_hex) {
        update_progress(args_progress, "pb", hex+1, n_hex);
      } else {
        int remainder = hex % int(ceil(double(n_hex)/100));
        if (remainder == 0 && !pb_markdown) {
          update_progress(args_progress, "pb", hex+1, n_hex);
        }
      }
    }

    hex_long1=hex_long[hex];
    hex_lat1=hex_lat[hex];    

    // loop through node locations
    for (node_loc1 = 0; node_loc1 < (n_node_locs-1); node_loc1++)
    {
      node_long1=loc_long[node_loc1];
      node_lat1=loc_lat[node_loc1];
      n_dups1=nDuplicates[node_loc1];
      for (node_loc2 = (node_loc1+1); node_loc2 < n_node_locs; node_loc2++)
      {
          bool intersects = collision_test_hex_ellipse(hex_long1, hex_lat1, hex_width,
                                                     node_long1, node_lat1,
                                                     loc_long[node_loc2], loc_lat[node_loc2],
                                                     eccentricity);
          if(intersects)
          {
            n_dups2=nDuplicates[node_loc2];
            for(i = 0; i < n_dups1; i++)
            {
                node1=DuplicateLists[node_loc1][i];
                for(j = 0; j < n_dups2; j++)
                {
                    node2=DuplicateLists[node_loc2][j];
                    ellipse=ellipse_list[node1][node2];
                    hex_edges[hex].push_back(ellipse);
				}
			}
		  }          
	  }
	}
  }
  
  // return list
  return Rcpp::List::create(Rcpp::Named("hex_edges") = hex_edges, Rcpp::Named("loc_long") = loc_long, Rcpp::Named("loc_lat") = loc_lat,
                            Rcpp::Named("nDuplicates") = nDuplicates);
}

//------------------------------------------------
// run main analysis
Rcpp::List pm_analysis_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
	
	// start timer
	//chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();
  
	// ------------------------------------------------------------------------------------------------
	// Convert Rcpp arguments to native c++ arguments
  
  vector<int> perm_group = rcpp_to_vector_int(args["perm_group"]);              // The permutation group of each observed edge
	vector<vector<double>> perm_list = rcpp_to_matrix_double(args["perm_list"]);  // The set of values in each permutation group
	vector<vector<int>> hex_edges = rcpp_to_matrix_int(args["hex_edges"]);        // The edges that intersect each hex
	int n_perms = rcpp_to_int(args["n_perms"]);                                   // Number of permutations to run
	bool report_progress = rcpp_to_bool(args["report_progress"]);                 // Whether to update progress bar
	bool pb_markdown = rcpp_to_bool(args["pb_markdown"]);                         //Whether to run in markdown-safe mode
	Rcpp::Function update_progress = args_functions["update_progress"];           // R function for updating progress bar
  
  
  // ------------------------------------------------------------------------------------------------
  // Derived values
  
  // numbers of edges and hexes
  int n_edge = int(perm_group.size());
  int n_hex = int(hex_edges.size());
  int n_breaks = int(perm_list.size());
  
  // size of list elements
  vector<int> perm_list_size(n_breaks);
  for (/*unsigned */int i = 0; i < n_breaks; ++i) {
    perm_list_size[i] = perm_list[i].size();
  }
  vector<int> hex_edges_size(n_hex);
  for (/*unsigned */int i = 0; i < n_hex; ++i) {
    hex_edges_size[i] = hex_edges[i].size();
  }
  
  // objects for storing results
  vector<double> edge_values(n_edge);
  vector<double> hex_values(n_hex);
  vector<double> ret_sum(n_hex);
  vector<vector<double>> ret_sum_sq(n_hex, vector<double>(n_hex));
  
  
  // ------------------------------------------------------------------------------------------------
  // Carry out permutation test
  
  print("Carrying out permutation test");
  
  // loop through permutations
  for (int perm = 0; perm < n_perms; ++perm) {
    
    // report progress
    if (report_progress) {
      if ((perm+1) == n_perms) {
        update_progress(args_progress, "pb", perm+1, n_perms);
      } else {
        int remainder = perm % int(ceil(double(n_perms)/100));
        if (remainder == 0 && !pb_markdown) {
          update_progress(args_progress, "pb", perm+1, n_perms);
        }
      }
    }
    
    // resample edge values
    for (int i = 0; i < n_edge; ++i) {
      int pg = perm_group[i] - 1;
      int rnd_index = sample2(1, perm_list_size[pg]) - 1;
      edge_values[i] = perm_list[pg][rnd_index];
    }
    
    // recalculate hex values
    fill(hex_values.begin(), hex_values.end(), 0.0);
    for (int h = 0; h < n_hex; ++h) {
      
      // skip hexes with no edges
      if (hex_edges_size[h] == 0) {
        continue;
      }
      
      // recalculate hex value
      for (int i = 0; i < hex_edges_size[h]; ++i) {
        hex_values[h] += edge_values[hex_edges[h][i] - 1];
      }
      hex_values[h] /= double(hex_edges_size[h]);
      
      // update running sums
      ret_sum[h] += hex_values[h];
      for (int j = 0; j < (h+1); ++j) {
        ret_sum_sq[h][j] += hex_values[j]*hex_values[h];
      }
      
    }
    
  }  // end loop over permutations
  
  // return list
  return Rcpp::List::create(Rcpp::Named("ret_sum") = ret_sum,
                            Rcpp::Named("ret_sum_sq") = ret_sum_sq);
}

//------------------------------------------------
// simulate from simple individual-based model
Rcpp::List sim_falciparum_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // extract model parameters into separate class
  Parameters parameters(args);
  //parameters.print_summary();
  
  // R functions
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // create simulation dispatcher object
  Dispatcher dispatcher(parameters, update_progress, args_progress);
  
  // carry out simulation
  dispatcher.simulate();
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
  print("simulation completed in", time_span.count(), "seconds\n");
  
  // return list
  return Rcpp::List::create(Rcpp::Named("daily_values") = dispatcher.daily_values.arr,
                            Rcpp::Named("genotypes") = dispatcher.genotypes.arr,
                            Rcpp::Named("indlevel_data") = dispatcher.indlevel_data.arr);
}
