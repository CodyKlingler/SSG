/* Copyright 2022, Gurobi Optimization, LLC */

/* This example formulates and solves the following simple MIP model:

     maximize    x +   y + 2 z
     subject to  x + 2 y + 3 z <= 4
                 x +   y       >= 1
                 x, y, z binary
*/


#include "include/gurobi_c++.h"
#include "include/lp_c++.h"
#include "SSG.h"

using namespace std;

int main2(SSG& gg, std::vector<bool> strategy)
{
  try {

    // Create an environment
    GRBEnv env = GRBEnv(true);
    env.set("LogFile", "mip1.log");
    env.start();

    // Create an empty model
    GRBModel model = GRBModel(env);


    // Set objective: maximize x + y + 2 z

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

      if(type == vertex_type::max){
        model.addConstr(vert[i],GRB_EQUAL,vert[selected]);
      }
      else if(type == vertex_type::min){
        model.addConstr(vert[i],GRB_LESS_EQUAL,vert[j]*(1-gg.beta));
        model.addConstr(vert[i],GRB_LESS_EQUAL,vert[k]*(1-gg.beta));
      }
      else if(type == vertex_type::ave){
        model.addConstr(vert[i],GRB_LESS_EQUAL,.5 *(1-gg.beta)* vert[k] + .5 * (1-gg.beta) * vert[j]);
      }
      else if(type == vertex_type::sink_max){
        model.addConstr(vert[i],GRB_EQUAL,1);
      }
      else if(type == vertex_type::sink_min){
        model.addConstr(vert[i],GRB_EQUAL,0);
      }
    }

    // find coefficient for each vertex in the objective function
    double coeffs[gg.n] = {1.0};

    for(int i = 0; i<gg.n; i++){
      if(gg.type[i] == vertex_type::max)
        coeffs[i] = 0.0;
    }

    GRBLinExpr objective_fn;
    objective_fn.addTerms(coeffs,vert,gg.n);
    
    model.setObjective(objective_fn, GRB_MAXIMIZE);

    // Optimize model
    model.optimize();

    cout << flush;

    for(int i = 0; i<gg.n; i++){
      cout << vert[i].get(GRB_DoubleAttr_X) << " ";
    } cout << endl;

    cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }

  return 0;
}

int gerby()
{
  try {

    // Create an environment
    GRBEnv env = GRBEnv(true);
    env.set("LogFile", "mip1.log");
    env.start();

    // Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables
    GRBVar X = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "X");
    GRBVar M = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "M");
    GRBVar A = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "A");
    //GRBVar XS = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "XS");
    //GRBVar MS = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "XM");

    #define XS 1
    #define MS 0

    model.setObjective(X + M + A + XS + MS, GRB_MINIMIZE);
    //model.addConstr(XS == 1, "max_sink");
    //model.addConstr(MS == 0, "min_sink");

    model.addConstr(X >= A, "c0");
    model.addConstr(X >= MS);

    model.addConstr(A == 0.5*(MS + XS));
    model.addConstr(M == A);



    // Optimize model
    model.optimize();

    cout << "X: " << X.get(GRB_DoubleAttr_X) << endl;
    cout << "M: " << M.get(GRB_DoubleAttr_X) << endl;
    cout << "A: " << A.get(GRB_DoubleAttr_X) << endl;

    cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }

  return 0;
}
