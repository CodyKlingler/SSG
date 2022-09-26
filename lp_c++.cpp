#include "include/lp_c++.h"
/*
#include "SSG.h"
#include "include/gurobi_c++.h"

using namespace std;


void reconstruct_strategy(SSG& gg, std::vector<bool> &strategy, std::vector<double> p){
  for(int i = 0; i<gg.n; i++){
    auto type = gg.type[i];
    if(type == vertex_type::min){
      int j = gg.outgoing_edge[i][0];
      int k = gg.outgoing_edge[i][1];
      strategy[i] = p[j] > p[k]; // selects smaller probability.
    }
  }
}

vector<double> optimize_min_LP(SSG& gg, std::vector<bool> &strategy)
{
  // returns prob
  vector<double> prob(0);

  try {
    // Create an envert[i]ronment
    GRBEnv env = GRBEnv(true);
    env.set(GRB_IntParam_OutputFlag,0);
    env.start();

    // Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables
    GRBVar vert[gg.n];
    for(int i = 0; i<gg.n; i++){
      vert[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
    }

    // Add constraints
    for(int i = 0; i<gg.n; i++){
      auto type = gg.type[i];
      int j = gg.outgoing_edge[i][0];
      int k = gg.outgoing_edge[i][1];
      int selected = gg.outgoing_edge[i][strategy[i]];

      double prob_move_next = (1-gg.beta);

      if(type == vertex_type::max){
        model.addConstr(vert[i],GRB_EQUAL,vert[selected]*prob_move_next);
      }
      else if(type == vertex_type::min){
        model.addConstr(vert[i],GRB_LESS_EQUAL,vert[j]*prob_move_next );
        model.addConstr(vert[i],GRB_LESS_EQUAL,vert[k]*prob_move_next );
      }
      else if(type == vertex_type::ave){
        model.addConstr(vert[i],GRB_EQUAL,.5 *prob_move_next * vert[k] + .5 * prob_move_next  * vert[j]);
      }
      else if(type == vertex_type::sink_max){
        model.addConstr(vert[i],GRB_EQUAL,1);
      }
      else if(type == vertex_type::sink_min){
        model.addConstr(vert[i],GRB_EQUAL,0);
      }
    }

    // find coefficient for each vertex in the objective function.
    // this may not be necessary.
    vector<double> coeffs(gg.n, 1.0);

    // set objective function
    GRBLinExpr objective_fn;
    objective_fn.addTerms(coeffs.data(),vert,gg.n);
    model.setObjective(objective_fn, GRB_MAXIMIZE);

    // solve model
    model.optimize();

    // build probability vector
    for(int i = 0; i<gg.n; i++){
      prob.push_back(vert[i].get(GRB_DoubleAttr_X));
    }
  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }
  
  return prob;
}
*/