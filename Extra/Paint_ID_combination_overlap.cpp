//
//  combination_overlap.cpp
//  
//  
//  Created by Nathalie Stroeymeyt on May 19th 2018
//  Copyright 2018 __University of Lausanne__. All rights reserved.
//  
//  DO NOT MODIFY

using namespace std;
#include <random>
#include <chrono>
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
DataFrame define_edge_list (List all_comb_list, int nb_locations, int nb_pairs){
  // declare variables that will be used in the function
  NumericVector colours; int nb_comb; int count; vector <int> comb_i; vector <int> comb_j;
  vector <int> nb_matches (nb_pairs,0);
  
  for (int location(0);location<nb_locations; location++){
    count = 0;
    colours = all_comb_list[location];
    nb_comb = colours.size();
    for (int comb1 (0);comb1<(nb_comb-1);comb1++){
      for (int comb2 (comb1+1);comb2<nb_comb;comb2++){
        if (location==0){
          comb_i.push_back(comb1+1);comb_j.push_back(comb2+1);
        }
        if (colours[comb1]==colours[comb2]){
          nb_matches[count]++;
        }
        count ++;
      }
    }
  }
  return DataFrame::create(_["comb_1"]          = comb_i             ,
                           _["comb_2"]          = comb_j             ,
                           _["nb_matches"]      = nb_matches         
  );
  
};

// [[Rcpp::export]]
NumericVector subsample_combinations (DataFrame edge_list, int max_common_spots, int nb_rand, vector <int> combs_ori){
  // define random engine generator
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();//set a seed linked with the time at which the program is runs; ensures random number sequence will always be differemt
  std::default_random_engine generator (seed); //set the random number generator using the desited seed

  //define variables that will be used in the function
  NumericVector comb_1;  NumericVector comb_2     ;  NumericVector nb_matches ;  int          nb_pairs  ; vector <int> combs;
  vector <int> candidate_comb1;vector <int> candidate_comb2; vector <int> candidate_matches; vector <int> candidate_combs;
  NumericVector best_node_list; vector <int> current_node_list; int best_list_size; int current_list_size;
  vector <int> to_remove; bool remove; bool in_list1; bool in_list2;
  int random;
  int new_node;
  
  // now perform random sampling of the list
  for (int i (0); i < nb_rand; i++){
    if ((i+1)%10==0){cout << i+1 << "th randomisation..." << endl;}
    // reinitialise comb_1, comb_2, nb_matches and nb_pairs
    comb_1     = edge_list["comb_1"];
    comb_2     = edge_list["comb_2"];
    nb_matches = edge_list["nb_matches"];
    nb_pairs   = comb_1.size();
    combs = combs_ori;
    // clear current_node_list
    current_node_list.clear();

    while(nb_pairs>0){
      // define random number generator that will draw a random node fropm vector combs
      std::uniform_int_distribution<int> distribution(0,(combs.size()-1));

      // draw random number
      random = distribution(generator);
      // use random number to define new focal node 
      new_node = combs[random];
      // remove that focal node from combs so it does not get redrawn later
      combs.erase(combs.begin() + random);
      // add focal node to current_node_list
      current_node_list.push_back(new_node);
      
      // Trim edge list: remove lines where both nodes are part of current_node_list, which are not informative
      candidate_comb1.clear();candidate_comb2.clear(); candidate_matches.clear();

      for (int pair (0); pair<nb_pairs; pair ++ ){
        in_list1 =0; in_list2=0;
        for (unsigned int curr_list_idx(0);curr_list_idx<current_node_list.size(); curr_list_idx++){
          if (comb_1[pair]==current_node_list[curr_list_idx]){in_list1=1;}
          if (comb_2[pair]==current_node_list[curr_list_idx]){in_list2=1;}
        }
        if (! ( (in_list1) && (in_list2))){
          candidate_comb1.push_back(comb_1[pair]);candidate_comb2.push_back(comb_2[pair]); candidate_matches.push_back(nb_matches[pair]);
        }
      }
      comb_1     = candidate_comb1;
      comb_2     = candidate_comb2;
      nb_matches = candidate_matches;
      nb_pairs   = comb_1.size();
      
      //Exclude all nodes that are too similar to focal node; i.e., nodes that have more than max_common_spots in common with focal node
      // First identify nodes to remove; 
      to_remove.clear();
      for (int pair (0); pair<nb_pairs; pair ++ ){
        if ((comb_1[pair]==new_node)&&(nb_matches[pair]>max_common_spots)){to_remove.push_back(comb_2[pair]);}
        if ((comb_2[pair]==new_node)&&(nb_matches[pair]>max_common_spots)){to_remove.push_back(comb_1[pair]);}
      }     
      // Second remove the "to_remove" nodes from combs
      candidate_combs.clear();  
      for (unsigned int combs_idx(0);combs_idx<combs.size();combs_idx++){
        remove=0;
        for (unsigned int to_remove_idx(0);to_remove_idx<to_remove.size();to_remove_idx++){
          if (combs[combs_idx]==to_remove[to_remove_idx]){
            remove=1;
          }
        }
        if (!remove){
          candidate_combs.push_back(combs[combs_idx]);
        }
      }
      combs = candidate_combs;
      // Third remove the "to_remove" nodes from edge_list 
      candidate_comb1.clear();candidate_comb2.clear(); candidate_matches.clear();
      for (int pair (0); pair<nb_pairs; pair ++ ){
        remove =0; //initialise: a priori edge does not need to be removed
        for (unsigned int j(0); j<to_remove.size(); j++){// then check if one of the nodes of the edge must be removed
          if ((comb_1[pair]==to_remove[j])||(comb_2[pair]==to_remove[j])){
            remove = 1;
          }
        }
        if (!remove){
          candidate_comb1.push_back(comb_1[pair]);candidate_comb2.push_back(comb_2[pair]); candidate_matches.push_back(nb_matches[pair]);
        }
      }
      comb_1     = candidate_comb1;
      comb_2     = candidate_comb2;
      nb_matches = candidate_matches;
      nb_pairs   = comb_1.size();
    }
    
    // Check if current list is better than previous list
    current_list_size = current_node_list.size(); best_list_size = best_node_list.size();
    if (current_list_size>best_list_size){
      best_node_list = current_node_list;
    }
  }
  return best_node_list;
};

// [[Rcpp::export]]
DataFrame sort_combinations(List selected_combination_lists,DataFrame edge_list){
  // read DataFrame
  NumericVector comb_1     =  edge_list["comb_1"];  
  NumericVector comb_2     =  edge_list["comb_2"];  
  NumericVector nb_matches =  edge_list["nb_matches"]; 
  int           nb_pairs   =  comb_1.size();
  
  // declare variables which will be used in the function
  bool already_sorted; 
  int nb_list = selected_combination_lists.size();
  NumericVector sorted_list; NumericVector current_list; vector <int> temp; vector <double> mean_overlaps; 
  double mean_overlap; double N_for_mean; double previous_overlap; int new_index;
  for (int list(0);list<nb_list;list++){
    // for all lists, make sure each combination is added sequentially to ensure minimum average overlap between all combinations
    // read current list
    current_list = selected_combination_lists[list];
    
    if (list==0){// if this is the first list read, initialise sorted_list
      sorted_list.push_back(current_list[0]); mean_overlaps.push_back(0);
    }
    
      // reduce current list to exclude elements already in sorted list
      temp.clear();
      for (unsigned int current (0); current < current_list.size(); current++){
        already_sorted = 0;
        for (unsigned int sorted (0); sorted< sorted_list.size(); sorted++){
          if (current_list[current]==sorted_list[sorted]){already_sorted = 1;}
        }
        if (!already_sorted){temp.push_back(current_list[current]);}
      }
      current_list = temp;
      // Now sort elements from current_list in order to maximise distance
      while ( current_list.size()>0){
        previous_overlap = -1; mean_overlap = 0;N_for_mean=0;
        for (unsigned int current (0); current < current_list.size(); current++){
          for (int pair (0); pair<nb_pairs;pair++){
            if (comb_1[pair]==current_list[current]){
              for (unsigned int sorted (0); sorted< sorted_list.size(); sorted++){
                if (comb_2[pair]==sorted_list[sorted]){
                  mean_overlap = mean_overlap + nb_matches[pair];N_for_mean++;   
                }
              }
            }else if (comb_2[pair]==current_list[current]){
              for (unsigned int sorted (0); sorted< sorted_list.size(); sorted++){
                if (comb_1[pair]==sorted_list[sorted]){
                  mean_overlap = mean_overlap + nb_matches[pair];N_for_mean++;
                }
              }
            }
          }
          mean_overlap = mean_overlap/N_for_mean;
          if (previous_overlap==-1){previous_overlap=mean_overlap;}
          if (mean_overlap<=previous_overlap){previous_overlap=mean_overlap;new_index=current;}
        }
        sorted_list.push_back(current_list[new_index]);
        mean_overlaps.push_back(previous_overlap);
        current_list.erase(current_list.begin() + new_index);
      }
    
  }
  return DataFrame::create(_["sorted_indices"]                    = sorted_list       ,
                           _["mean_overlap"]                      = mean_overlaps
  );
};
