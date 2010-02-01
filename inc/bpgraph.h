#ifndef BP_GRAPH_H_
#define BP_GRAPH_H_
#include <iostream>
#include <string>
#include <boost/config.hpp>
#include <iostream>                      // for std::cout
#include <utility>                       // for std::pair
#include <algorithm>                     // for std::for_each
#include <vector>                        // 
#include <boost/utility.hpp>             // for boost::tie
#include <boost/graph/graph_traits.hpp>  // for boost::graph_traits
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include "node.h"
using namespace std;
using namespace boost;
 
struct GraphProperty {
    string name;
    int numVal; 
};
typedef pair<LD, LD> EdgeProperty; 
typedef adjacency_list<setS, setS, bidirectionalS, property<vertex_name_t, Node>, EdgeProperty, GraphProperty> Graph; 


class BPGraph {
      
};
#endif /*BP_GRAPH_H_*/

