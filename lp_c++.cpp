#pragma once
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

int main2()
{
  try {

    // Create an environment
    GRBEnv env = GRBEnv(true);
    env.set("LogFile", "mip1.log");
    env.start();

    // Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables
    GRBVar x = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x");
    GRBVar y = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y");
    GRBVar z = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z");

    // Set objective: maximize x + y + 2 z
    model.setObjective(x + y + 2 * z, GRB_MAXIMIZE);

    // Add constraint: x + 2 y + 3 z <= 4
    model.addConstr(x + 2 * y + 3 * z <= 4, "c0");

    // Add constraint: x + y >= 1
    model.addConstr(x + y >= 1, "c1");

    // Optimize model
    model.optimize();

    cout << x.get(GRB_StringAttr_VarName) << " "
         << x.get(GRB_DoubleAttr_X) << endl;
    cout << y.get(GRB_StringAttr_VarName) << " "
         << y.get(GRB_DoubleAttr_X) << endl;
    cout << z.get(GRB_StringAttr_VarName) << " "
         << z.get(GRB_DoubleAttr_X) << endl;

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


std::vector<bool> optimize_min_lp(std::vector<bool> s, SSG g){

  #define S_TYPE GRB_CONTINUOUS //GRB_BINARY

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

    std::vector<GRBVar> vars(g.n);

    for(int i = 0; i<g.n; i++){
      GRBVar current = model.addVar(0.0, 1.0, 0.0, GRB_STA)
      vars.push_back
    }


    //GRBVar XS = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "XS");
    //GRBVar MS = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "XM")

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

}


