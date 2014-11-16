#include <iostream>
#include <Eigen/Eigenvalues> 
#include "graph_manip/graph_data_structures.hpp"
void print_help()
{
  std::cout <<
    "calc_sis: Utility for calculating the SIS effective strenth of a virus" <<
    " in a static network" << std::endl <<
    "Usage:" << std::endl << std::endl <<
    "  calc_sis <graph_file> <transmission_prob> <healing_prob>" << std::endl <<
    "   <graph_file> path to file that contains the graph to test" << std::endl <<
    "   <transmission_prob> Probability of transmission (between 0 and 1)" << std::endl <<
    "   <healing_prob> Probability of healing (between 0 and 1)" << std::endl;
}
 

int main(int argc, char** argv)
{
  using namespace csc791;
  using namespace Eigen;
  if (argc >= 2)
  {
    if (std::string(argv[1])=="--help" || 
        std::string(argv[1])=="-h")
    {
      print_help();
      return 0;
    }
  }

  if (argc < 4)
  {
    std::cerr << "Error, you must have 4 arguments. See the following usage" << std::endl;
    print_help();
    return -1;
  }
  std::string infile = argv[1];
  std::stringstream s(argv[2]);
  float xmission_prob,healing_prob;
  s >> xmission_prob;
  s.clear();
  s.str(std::string(argv[3]));
  s >> healing_prob;

  if (xmission_prob > 1 || xmission_prob < 0)
  {
    std::cerr << "Error, xmission prob must be between 0 and 1, entered as "
	      << xmission_prob << std::endl;
    return -2;
  }
  if (healing_prob > 1 || healing_prob < 0)
  {
    std::cerr << "Error, healing prob must be between 0 and 1, entered as "
	      << healing_prob << std::endl;
    return -3;
  }
  EdgeList el = read_edge_list(infile.c_str());
  AdjacencyMatrix am = edge_list2adjacency_matrix(el,false);
  //std::cout << am.cols() << "," << am.rows() << std::endl;

  size_t n = am.cols();

  EigenSolver<MatrixXd> es(MatrixXd(am),false);
  auto lambda = es.eigenvalues();
  
  //find the max
  double max_eigen = -std::numeric_limits<double>::infinity();
  for (int i = 0; i <n ; i++)
  {
    auto eigen_val = abs(lambda[i]);
    
    max_eigen = std::max(max_eigen,eigen_val);
  }

  double effective_strength = max_eigen*xmission_prob/healing_prob;
  std::cout << effective_strength << std::endl;
}
