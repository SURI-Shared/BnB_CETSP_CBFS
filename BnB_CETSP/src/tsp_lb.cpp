#include "tsp_lb.h"

ILOSTLBEGIN

void subgraph_tsp_solver::solve_cplex_miqcp()
{
    int sizeProblem = m_size_subgraph;
    // CPLEX environment. Takes care of everything, including memory management for CPLEX objects.
    IloEnv env;

    // CPLEX model. We put variables and constraints in it!
    IloModel model(env);

    IloNumArray xCoord(env);
    IloNumArray yCoord(env);
    IloNumArray zCoord(env);
    IloNumVarArray x(env, sizeProblem, -IloInfinity, IloInfinity); 
    IloNumVarArray y(env, sizeProblem, -IloInfinity, IloInfinity); 
    IloNumVarArray z(env, sizeProblem, -IloInfinity, IloInfinity);

    //variables
    //add variable f ( head of cones )
    //================================================================
    IloNumVarArray f(env, sizeProblem, 0, IloInfinity);
    char var[100];
    for (int i = 0; i < sizeProblem; i++) {
        sprintf(var, "f(%d)", i );
        f[i].setName(var);
        //model.add(f[ i ]);
    }
    model.add(f);
    //getchar();
    //add variables corresponding to coordinates x, y and z
    //================================================================
    //variable x
    for (int i = 0; i < sizeProblem; i++) {
        sprintf( var, "x(%d)", i );
        x[i].setName( var );
        //model.add( x[ i ] );
    }
    model.add(x);

    //variable y
    for (int i = 0; i < sizeProblem; i++) {
        sprintf(var, "y(%d)", i );
        y[i].setName(var);
        //model.add( y[ i ] );
    }
    model.add(y);

    //variable z
    for (int i = 0; i < sizeProblem; i++) {
        sprintf(var, "z(%d)", i );
        z[i].setName(var);
        //model.add( z[ i ] );
    }
    model.add(z);

    //add auxiliary variables w, u, v, s, t and q
    //================================================================
    IloNumVarArray w(env, sizeProblem, -IloInfinity, IloInfinity);
    for (int i = 0; i < sizeProblem; i++) {
        sprintf(var, "w(%d)", i );
        w[i].setName(var);
        //model.add(w[ i ]);
    }
    model.add(w);

    IloNumVarArray u(env, sizeProblem, -IloInfinity, IloInfinity);
    for (int i = 0; i < sizeProblem; i++) {
        sprintf(var, "u(%d)", i );
        u[i].setName(var);
        //model.add(u[ i ]);
    }
    model.add(u);

    IloNumVarArray v(env, sizeProblem, -IloInfinity, IloInfinity);
    for (int i = 0; i < sizeProblem; i++) {
        sprintf(var, "v(%d)", i );
        v[i].setName(var);
        //model.add(v[ i ]);
    }
    model.add(v);

    IloNumVarArray s(env, sizeProblem, -IloInfinity, IloInfinity);
    for (int i = 0; i < sizeProblem; i++) {
        sprintf(var, "s(%d)", i );
        s[i].setName(var);
        //model.add(s[ i ]);
    }
    model.add(s);

    IloNumVarArray t(env, sizeProblem, -IloInfinity, IloInfinity);
    for (int i = 0; i < sizeProblem; i++) {
        sprintf(var, "t(%d)", i );
        t[i].setName(var);
        //model.add(t[ i ]);
    }
    model.add(t);

    IloNumVarArray q(env, sizeProblem, -IloInfinity, IloInfinity);
    for (int i = 0; i < sizeProblem; i++) {
        sprintf(var, "q(%d)", i );
        q[i].setName(var);
        //model.add(q[ i ]);
    }
    model.add(q);

    // TSP variables
    // Create variable t[0] and fix it to value 1
    // This breaks symmetry, because it fixes node 0 as the starting node of the tour
    //verord[0] = IloNumVar(env, 1, 1, IloNumVar::Int, "verord(0)");

    // Create variables t[1], ..., t[n]
    IloArray<IloNumVarArray> verord(env, sizeProblem);
    for(int i = 0; i < sizeProblem; ++i) 
    {
        verord[i] = IloNumVarArray(env, sizeProblem);
        for (int j = 0; j < sizeProblem; ++j)
        {
            sprintf(var, "verord(%d)(%d)", i, j);
            //if (i == 0 && j == 0)
                //verord[i][j] = IloNumVar(env, 1, 1, IloNumVar::Bool, var);
            verord[i][j] = IloNumVar(env, 0, 1, IloNumVar::Bool, var);
        }
        model.add(verord[i]);
    }

    // Flow type veriables
    // IloInt numArcs  = numNodes * numNodes;
    // IloInt vNumVars = (numNodes-1) * numArcs;
    // IloNumVarArray flowvar(env, vNumVars, 0, IloInfinity);
    // for (k = 1; k < numNodes; ++k) {
    //     for(i = 0; i < numNodes; ++i) {
    //         for(j = 0; j < numNodes; ++j) {
    //             sprintf(var, "flowvar.%d.%d.%d", (int) k, (int) i, (int) j); 
    //             flowvar[(k-1)*numArcs + i*numNodes + j].setName(varName);
    //         }
    //     }
    // }
    // model.add(flowvar);

    // Create variables x
    IloArray<IloNumVarArray> edge(env, sizeProblem);
    for(int i = 0; i < sizeProblem; ++i) 
    {
        edge[i] = IloNumVarArray(env, sizeProblem);
        for(int j = 0; j < sizeProblem; ++j) 
        {
            sprintf(var, "edge(%d)(%d)", i, j);
            edge[i][j] = IloNumVar(env, 0, 1, IloNumVar::Bool, var);
        }
        model.add(edge[i]);
    }

    // Connection variables
    IloNumVarArray pos_coord_x(env, sizeProblem, -IloInfinity, IloInfinity);
    for (auto i = 0; i < sizeProblem; ++i)
    {
        sprintf(var, "pos_coord_x(%d)", i);
        pos_coord_x[i].setName(var);
        model.add(pos_coord_x[i]);
    }

    IloNumVarArray pos_coord_y(env, sizeProblem, -IloInfinity, IloInfinity);
    for (auto i = 0; i < sizeProblem; ++i)
    {
        sprintf(var, "pos_coord_y(%d)", i);
        pos_coord_y[i].setName(var);
        model.add(pos_coord_y[i]);
    }

    IloNumVarArray pos_coord_z(env, sizeProblem, -IloInfinity, IloInfinity);
    for (auto i = 0; i < sizeProblem; ++i)
    {
        sprintf(var, "pos_coord_z(%d)", i);
        pos_coord_z[i].setName(var);
        model.add(pos_coord_z[i]);
    }

    IloNumVarArray pos_radius(env, sizeProblem, 0, IloInfinity);
    for (auto i = 0; i < sizeProblem; ++i)
    {
        sprintf(var, "pos_radius(%d)", i);
        pos_radius[i].setName(var);
        model.add(pos_radius[i]);
    }

    //Set objetive function
    //=============================================
    IloExpr FO(env);
    for (int i = 0; i < sizeProblem; i++) {
        FO += f[ i ];
    }

    model.add(IloMinimize(env, FO));

    //set SOC constraints
    //===============================================
    for (int j = 0; j < sizeProblem; j++) {
        IloRange r = ( - f[ j ]*f[ j ] + w[ j ]*w[ j ] + u[ j ]*u[ j ] + v[ j ]*v[ j ] <= 0);
        char c[100];
        sprintf(c, "cone%d", j );
        r.setName(c);
        model.add(r);
    }

    for (int i = 0; i < sizeProblem; i++) {
        IloRange r = ( s[ i ]*s[ i ] + t[ i ]*t[ i ] + q[ i ]*q[ i ] - pos_radius[i]*pos_radius[i]<= 0 );
        char c[100];
        sprintf(c, "q%d", i );
        r.setName(c);
        model.add(r);
    }

    //radius constraints	
    //=============================================================
    for (int i = 0; i < sizeProblem; i++) {
        IloRange r = ( s[ i ] + x[ i ] - pos_coord_x[i] == 0);
        char c[100];
        sprintf(c, "b%d", i );
        r.setName(c);
        model.add(r);
    }

    for (int i = 0; i < sizeProblem; i++) {
        IloRange r = ( t[ i ] + y[ i ] - pos_coord_y[i] == 0);
        char c[100];
        sprintf(c, "b%d", i + sizeProblem);
        r.setName(c);
        model.add(r);
    }

    for (int i = 0; i < sizeProblem; i++) {
        IloRange r = ( q[ i ] + z[ i ] - pos_coord_z[i] == 0);
        char c[100];
        sprintf(c, "b%d", i + 2*sizeProblem);
        r.setName(c);
        model.add(r);
    }

    //#########################################################################################################################
    //separating first constraint
    IloRange isol_1 = ( w[ 0 ] - x[ sizeProblem - 1 ] + x[ 0 ] == 0 );
    char c[100];
    sprintf(c, "c%d", 0);
    isol_1.setName(c);
    model.add( isol_1 );

    for (int j = 1; j < sizeProblem; j++) {
        IloRange r = ( w[ j ] - x[ j - 1 ] + x[ j ] == 0 );
        char c[100];
        sprintf(c, "c%d", j);
        r.setName(c);
        model.add(r);
    }

    IloRange isol_2 = ( u[ 0 ] - y[ sizeProblem - 1 ] + y[ 0 ] == 0 );
    char c2[100];
    sprintf(c2, "c%d", sizeProblem );
    isol_2.setName(c2);
    model.add( isol_2 );

    for (int j = 1; j < sizeProblem; j++) {
        IloRange r = ( u[ j ] - y[ j - 1 ] + y[ j ] == 0 );
        char c[100];
        sprintf(c, "c%d", j + sizeProblem);
        r.setName(c);
        model.add(r);
    }

    IloRange isol_3 = ( v[ 0 ] - z[ sizeProblem - 1 ] + z[ 0 ] == 0 );
    char c3[100];
    sprintf(c3, "c%d", sizeProblem );
    isol_3.setName(c3);
    model.add( isol_3 );

    for (int j = 1; j < sizeProblem; j++) {
        IloRange r = ( v[ j ] - z[ j - 1 ] + z[ j ] == 0 );
        char c[100];
        sprintf(c, "c%d", j + 2*sizeProblem);
        r.setName(c);
        model.add(r);
    }

    // Create inbound constraints
    //IloRangeArray inbound_arcs(env, sizeProblem)
    // for(int i = 0; i < sizeProblem; ++i) 
    // {
    //     IloExpr expr(env);
    //     for (int j = 0; j < i; ++j) expr += edge[j][i];
    //     for (int j = i + 1; j < sizeProblem; ++j) expr += edge[j][i];
    //     //char c[100];
    //     //sprintf(c, "inbound_%d", i);
    //     //nbound_arcs[i] = IloRange(env, 1, expr, 1, c);
    //     model.add(expr == 1);
    //     expr.clear(); // Clean expr
    // }
    //model.add(inbound_arcs);

    // Create outbound constraints
    //IloRangeArray outbound_arcs(env, sizeProblem);
    // for(int i = 0; i < sizeProblem; ++i) 
    // {
    //     IloExpr expr(env);
    //     for (int j = 0; j < i; ++j) expr += edge[i][j];
    //     for (int j = i + 1; j < sizeProblem; ++j) expr += edge[i][j];
    //     //char c[100];
    //     //sprintf(c, "outbound_%d", i);
    //     //outbound_arcs[i] = IloRange(env, 1, expr, 1, c);
    //     model.add(expr == 1);
    //     expr.clear(); // Clean expr
    // }
    //model.add(outbound_arcs);

    // Create flow constraints
    // int flowpara = 0;
    // for (int k = 1; k < sizeProblem; ++k)
    // {
    //     for (int i = 0; i < sizeProblem; ++i)
    //     {
    //         IloExpr expr(env);
    //         if (i == 0) flowpara = 1;
    //         if (k == i) flowpara = -1;
    //         for (int j = 0; j < sizeProblem; ++j)
    //         {
    //             expr += flowvar[(k-1)*numArcs + i*numNodes + j] - 
    //                 flowvar[(k-1)*numArcs + j*numNodes + i];
    //         }
    //         model.add(expr == flowpara);
    //         expr.clear();
    //     }
    // }
    // for (int k = 1; k < sizeProblem; ++k)
    // {
    //     for (int i = 0; i < sizeProblem; ++i)
    //     {
    //         for (int j = 0; j < sizeProblem; ++j)
    //         {
    //             IloExpr expr(env);
    //             expr += flowvar[(k-1)*numArcs + i*numNodes + j] - edge[i][j];
    //             model.add(expr <= 0);
    //             expr.clear();
    //         }
    //     }
    // }

    //mtz[0] = IloRangeArray(env);
    // We then continue normally for all other i > 0
    // for(int i = 1; i < sizeProblem; ++i) 
    // {
    //     //mtz[i] = IloRangeArray(env, sizeProblem);
    //     for(int j = 1; j < sizeProblem; ++j) 
    //     {
    //         IloExpr expr(env);
    //         for (int k = 1; k < sizeProblem; ++k)
    //         {
    //             expr += (verord[i][k] - verord[j][k])*k;
    //         }
    //         expr += (static_cast<int>(sizeProblem) - 1) * edge[i][j] + (static_cast<int>(sizeProblem) - 3) * edge[j][i];
    //         //name << "mtz_" << i << "_" << j;
    //         //mtz[i][j] = IloRange(env, -IloInfinity, expr, sizeProblem - 1, name.str().c_str());
    //         //name.str(""); // Clean name
    //         model.add(expr <= static_cast<int>(sizeProblem) - 2);
    //         expr.clear(); // Clean expr
    //     }
    // }

    // position constraint
    for (int i = 0; i < 1; ++i)
    {
        IloRange r = (verord[0][0] == 1);
        model.add(r);
        for (int j = 1; j < sizeProblem; ++j)
        {
            IloRange r = (verord[0][j] == 0);
            model.add(r);
        }
    }

    // vertex order variable constraint (1 vertex takes 1 position)
    for (int i = 1; i < sizeProblem; ++i)
    {
        IloExpr expr(env);
        for (int j = 1; j < sizeProblem; ++j) expr += verord[i][j];
        model.add(expr == 1);
        expr.clear();
    }
    for (int i = 1; i < sizeProblem; ++i)
    {
        IloExpr expr(env);
        for (int j = 1; j < sizeProblem; ++j) expr += verord[j][i];
        model.add(expr == 1);
        expr.clear();
    }

    if (m_precedence_sequence.size() > 1)
    {
        int from, to;
        for (int i = 0; i < m_precedence_sequence.size(); i++)
        {
            from = m_precedence_sequence[i];
            int prec_dist = -1;
            for (int j = i+1; j < m_precedence_sequence.size(); j++)
            {
                to = m_precedence_sequence[j];
                IloExpr expr(env);
                for (int k = 1; k < sizeProblem; k++) 
                {
                    expr += (verord[from][k] - verord[to][k]) * k;
                }
                model.add(expr <= prec_dist);
                expr.clear();
                prec_dist--;
            }
        }
    }

    // Connection constraints
    for (int i = 0; i < sizeProblem; i++)
    {
        IloExpr expr(env);
        for (int j = 0; j < sizeProblem; j++)
        {
            expr += m_data->getCoordx(j) * verord[i][j];
        }
        model.add(expr == pos_coord_x[i]);
        expr.clear();
    }
    for (int i = 0; i < sizeProblem; i++)
    {
        IloExpr expr(env);
        for (int j = 0; j < sizeProblem; j++)
        {
            expr += m_data->getCoordy(j) * verord[i][j];
        }
        model.add(expr == pos_coord_y[i]);
        expr.clear();
    }
    for (int i = 0; i < sizeProblem; i++)
    {
        IloExpr expr(env);
        for (int j = 0; j < sizeProblem; j++)
        {
            expr += m_data->getCoordz(j) * verord[i][j];
        }
        model.add(expr == pos_coord_z[i]);
        expr.clear();
    }
    for (int i = 0; i < sizeProblem; i++)
    {
        IloExpr expr(env);
        for (int j = 0; j < sizeProblem; j++)
        {
            expr += m_data->getRadius(j) * verord[i][j];
        }
        model.add(expr == pos_radius[i]);
        expr.clear();
    }

    //expr.end();

    IloCplex cpx(model);
    cpx.setParam( IloCplex::Threads, 8 );

    // cpx.setParam(IloCplex::Param::OptimalityTarget, 1);
    // cpx.setParam(IloCplex::Param::MIP::Cuts::Gomory, 2);
    // cpx.setParam( IloCplex::Param::MIP::Cuts::LiftProj, 3);

    cpx.setParam(IloCplex::Param::MIP::Limits::GomoryCand, 10000);
    cpx.setParam(IloCplex::Param::MIP::Limits::GomoryPass, 10);

    // cpx.setParam(IloCplex::Param::Tune::Measure, CPX_TUNE_AVERAGE);

    //cpx.setParam(IloCplex::Param::Benders::Strategy,
    //                  IloCplex::BendersFull);
    // Write out the auto-generated annotation.
    //cpx.writeBendersAnnotation("benders.ann");

    bool solved = false;
    try 
    {
        // Try to solve with CPLEX (and hope it does not raise an exception!)
        solved = cpx.solve();
        // IloInt tunestat = cpx.tuneParam();
        // if ( tunestat == IloCplex::TuningComplete)
        //     cout << "Tuning complete." << endl;
        // else if ( tunestat == IloCplex::TuningAbort)
        //     cout << "Tuning abort." << endl;
        // else if ( tunestat == IloCplex::TuningTimeLim)
        //     cout << "Tuning time limit." << endl;
        // else
        //     cout << "Tuning status unknown." << endl;

        // cpx.writeParam("param.out");
    } 
    catch(const IloException& e) 
    {
        cerr << "\n\nCPLEX Raised an exception:\n";
        cerr << e << "\n";
        env.end();
        throw;
    }
    m_solution = cpx.getBestObjValue();
}

void subgraph_tsp_solver::solve_cplex_weak()
{
    int n = m_size_subgraph;
    if (n <= 4) return;
    
    // CPLEX environment. Takes care of everything, including memory management for CPLEX objects.
    IloEnv env;

    // CPLEX model. We put variables and constraints in it!
    IloModel model(env);

    // Model:
    //
    // BINARY VARIABLE x[i][j]    For all i,j = 0, ..., n - 1
    //    x[i][j] == 1            If arc (i,j) is selected
    //    x[i][j] == 0            Otherwise
    //
    // INTEGER VARIABLE t[i]      For all i = 0, ..., n - 1
    //    t[i] == k               Iff node i is the k-th node in the tour
    //    t[0] == 1
    //    t[i] in [2, ..., n]     For all i = 1, ... n - 1
    //
    // OBJECTIVE FUNCTION
    //    MIN sum((i,j), c[i][j] * x[i][j])
    //
    // CONSTRAINTS
    //    1) sum(j, x[j][i]) == 1                    For all i
    //    2) sum(j, x[i][j]) == 1                    For all i
    //    3) t[i] - t[j] + 1 <= n * (1 - x[i][j])    For all i,j = 1, ..., n - 1
    //       Can be written as:
    //       t[i] - t[j] + n * x[i][j] <= n - 1

    // Variables
    IloArray<IloNumVarArray> x(env, n);
    IloNumVarArray t(env, n);

    // Constraints
    IloRangeArray inbound_arcs(env, n);  // Constraints 1)
    IloRangeArray outbound_arcs(env, n); // Constraints 2)
    IloArray<IloRangeArray> mtz(env, n); // Constraints 3)

    // We use this stringstream to create variable and constraint names
    stringstream name;

    // Create variable t[0] and fix it to value 1
    // This breaks symmetry, because it fixes node 0 as the starting node of the tour
    t[0] = IloNumVar(env, 1, 1, IloNumVar::Int, "t_0");

    // Create variables t[1], ..., t[n]
    for(auto i = 1u; i < n; ++i) 
    {
        name << "t_" << i;
        t[i] = IloNumVar(env, 2, n, IloNumVar::Int, name.str().c_str());
        name.str(""); // Clean name
    }

    // Create variables x
    for(auto i = 0u; i < n; ++i) 
    {
        x[i] = IloNumVarArray(env, n);
        for(auto j = 0u; j < n; ++j) 
        {
            name << "x_" << i << "_" << j;
            x[i][j] = IloNumVar(env, 0, 1, IloNumVar::Bool, name.str().c_str());
            name.str(""); // Clean name
        }
    }

    IloExpr expr(env);

    // Create constraints 1)
    for(auto i = 0u; i < n; ++i) 
    {
        for(auto j = 0u; j < n; ++j) 
        {
            expr += x[j][i];
        }

        name << "inbound_" << i;
        inbound_arcs[i] = IloRange(env, 1, expr, 1, name.str().c_str());
        name.str(""); // Clean name
        expr.clear(); // Clean expr
    }

    // Add constraints 1) to the model
    model.add(inbound_arcs);

    // Create constraints 2)
    for(auto i = 0u; i < n; ++i) 
    {
        for(auto j = 0u; j < n; ++j) 
        {
            expr += x[i][j];
        }

        name << "outbound_" << i;
        outbound_arcs[i] = IloRange(env, 1, expr, 1, name.str().c_str());
        name.str(""); // Clean name
        expr.clear(); // Clean expr
    }

    // Add constraints 2) to the model
    model.add(outbound_arcs);

    // Create constraints 3)
    // The constraint is for i = 1,...,n and therefore we add empty constraints for i == 0
    mtz[0] = IloRangeArray(env);
    // We then continue normally for all other i > 0
    for(auto i = 1u; i < n; ++i) 
    {
        mtz[i] = IloRangeArray(env, n);
        for(auto j = 1u; j < n; ++j) 
        {
            expr = t[i] - t[j] + static_cast<int>(n) * x[i][j];

            name << "mtz_" << i << "_" << j;
            mtz[i][j] = IloRange(env, -IloInfinity, expr, n - 1, name.str().c_str());
            name.str(""); // Clean name
            expr.clear(); // Clean expr
        }
        // Add constraints 3)[i] to the model
        model.add(mtz[i]);
    }

    // Create objective function
    for(auto i = 0; i < n; ++i) 
    {
        for(auto j = 0; j < n; ++j) 
        {
            expr += m_subgraph->get_relaxed_dist(i, j) * x[i][j];
        }
    }
    IloObjective obj(env, expr, IloObjective::Minimize);

    // Add the objective function to the model
    model.add(obj);

    // Free the memory used by expr
    expr.end();

    // Create the solver object
    IloCplex cplex(model);
    cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::TiLim,m_time_lim);
    cplex.setParam( IloCplex::Threads, 1 );
    //cplex.setParam( IloCplex::ParallelMode, 1 );

    // Export model to file (useful for debugging!)
    //cplex.exportModel("model.lp");

    bool solved = false;

    try 
    {
        // Try to solve with CPLEX (and hope it does not raise an exception!)
        solved = cplex.solve();
    } 
    catch(const IloException& e) 
    {
        cerr << "\n\nCPLEX Raised an exception:\n";
        cerr << e << "\n";
        env.end();
        throw;
    }

    m_solution = cplex.getBestObjValue();

    env.end();
}

void subgraph_tsp_solver::solve_concorde()
{
    if (m_size_subgraph <= 4) return;

    vector<vector<int>> int_dist_matrix;
    int_dist_matrix.resize(m_size_subgraph);
    for (int i = 0; i < m_size_subgraph; i++)
    {
        int_dist_matrix[i].resize(m_size_subgraph, 0);
    }

    for (int i = 0; i < m_size_subgraph; i++)
    {
        for (int j = i+1; j < m_size_subgraph; j++)
        {
            int_dist_matrix[i][j] = m_subgraph->get_relaxed_dist(i,j) * m_int_factor;
            if (int_dist_matrix[i][j] < 0)
                int_dist_matrix[i][j] = 0;
            int_dist_matrix[i][j] += m_dist_addon;
            int_dist_matrix[j][i] = int_dist_matrix[i][j];
        }
    }

    //TSP for more than 4 cities
    int rval = 0; //Concorde functions return 1 if something fails
    double szeit; //Measure cpu time
    double optval; //Value of the optimal tour
    double *in_val = (double *) NULL; //Can be used to specify an initial upperbound (it can be NULL)
    double *timebound = (double *) NULL;; //Run time limit
    double time_lim = m_time_lim;
    int success; //1 if the run finished normally, and set to 0 if the search was terminated early (by hitting some predefined limit) 
    int foundtour; //1 if a tour has been found (if success is 0, then it may not be the optimal tour)   
    int hit_timebound; //1 if timebound was reached
    int *in_tour = (int *) NULL; //Gives a starting tour in node node node format (it can be NULL)
    int *out_tour = (int *) NULL; //Optimal tour (it can be NULL, if it is not NULL then it should point to an array of length at least ncount).  
    char *name = (char *) NULL; //Specifes a char string that will be used to name various files that are written during the branch and bound search
    static int silent = 1; //Suppress most output if set to a nonzero value
    double lbval = -1;

    CCrandstate rstate;
    int seed = rand();
    CCutil_sprand(seed, &rstate); //Initialize the portable random number generator
    int ncount = int_dist_matrix.size(); //Number of nodes (cities)
    int ecount = (ncount * (ncount - 1)) / 2; //Number of edges
    int *elist = new int[ecount * 2]; //Array giving the ends of the edges (in pairs)
    int *elen = new int[ecount]; //Array giving the weights of the edges
    int edge = 0;
    int edgeWeight = 0;
    for (int i = 0; i < ncount; i++) 
    {
        for (int j = i + 1; j < ncount; j++) 
        {
            if (i != j) 
            {
                elist[edge] = i;
                elist[edge + 1] = j;
                elen[edgeWeight] = int_dist_matrix[i][j];
                if (elen[edgeWeight] <= 0)
                    elen[edgeWeight] = 1;
                edgeWeight++;
                edge = edge + 2;
            }
        }
    }

    out_tour = CC_SAFE_MALLOC (ncount, int);
    name = CCtsp_problabel("_");
    CCdatagroup dat;

    // Populate initial solution if available
    if (!m_init_sequence.empty())
    {
        int temp_seq_array[m_init_sequence.size()];
        for (int i = 0; i < m_init_sequence.size(); i++)
        {
            temp_seq_array[i] = m_init_sequence[i];
        }
        in_tour = temp_seq_array;
    }

    //Initialize a CCdatagroup
    CCutil_init_datagroup (&dat);

    //Convert a matrix of edge lengths to a CCdatagroup
    rval = CCutil_graph2dat_matrix (ncount, ecount, elist, elen, 0, &dat);

    //Solves the TSP over the graph specified in the datagroup
    try
    {
        rval = CCtsp_solve_dat (ncount, &dat, in_tour, out_tour, in_val, &optval, &success, &foundtour, name, timebound, &hit_timebound, silent, &rstate);
        //rval = CCtsp_solve_sparse (ncount, ecount, elist, elen, in_tour, out_tour, in_val, &optval, &success, &foundtour, name, timebound, &hit_timebound, silent, &rstate);
        //rval = CCtsp_solve_dat_wlb (ncount, &dat, in_tour, out_tour, in_val, &optval, &success, &foundtour, name, &time_lim, &hit_timebound, silent, &rstate, &lbval);
    }
    catch (int n)
    {
        cout << "Contour exception." << endl;
    }

    if (success == 1)
    {
       m_solution = (optval - ncount * m_dist_addon) / m_int_factor;
    }
    //cout << "TEST LB VALUE: " << lbval << endl;
    else if (lbval > 0)
    {
        m_solution = (lbval - ncount * m_dist_addon) / m_int_factor;
    }

    for (int i = 0; i < m_size_subgraph; i++)
    {
        m_sequence[i] = out_tour[i];
    }
    //szeit = CCutil_zeit();
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CC_IFFREE (out_tour, int);
    CC_IFFREE (name, char);
}

void subgraph_tsp_solver::solve_concorde(vector<vector<double>>& dist_mat)
{
    m_size_subgraph = dist_mat.size();
    if (m_size_subgraph <= 4) return;
    vector<vector<int>> int_dist_matrix;
    int_dist_matrix.resize(m_size_subgraph);
    for (int i = 0; i < m_size_subgraph; i++)
    {
        int_dist_matrix[i].resize(m_size_subgraph, 0);
    }

    for (int i = 0; i < m_size_subgraph; i++)
    {
        for (int j = i+1; j < m_size_subgraph; j++)
        {
            int_dist_matrix[i][j] = static_cast<int>(dist_mat[i][j] * m_int_factor);
            if (int_dist_matrix[i][j] < 0)
                int_dist_matrix[i][j] = 0;
            int_dist_matrix[i][j] += m_dist_addon;
            int_dist_matrix[j][i] = int_dist_matrix[i][j];
        }
    }

    //TSP for more than 4 cities
    int rval = 0; //Concorde functions return 1 if something fails
    double szeit; //Measure cpu time
    double optval; //Value of the optimal tour
    double *in_val = (double *) NULL; //Can be used to specify an initial upperbound (it can be NULL)
    double *timebound = (double *) NULL;; //Run time limit
    double time_lim = m_time_lim;
    int success; //1 if the run finished normally, and set to 0 if the search was terminated early (by hitting some predefined limit) 
    int foundtour; //1 if a tour has been found (if success is 0, then it may not be the optimal tour)   
    int hit_timebound; //1 if timebound was reached
    int *in_tour = (int *) NULL; //Gives a starting tour in node node node format (it can be NULL)
    int *out_tour = (int *) NULL; //Optimal tour (it can be NULL, if it is not NULL then it should point to an array of length at least ncount).  
    char *name = (char *) NULL; //Specifes a char string that will be used to name various files that are written during the branch and bound search
    static int silent = 1; //Suppress most output if set to a nonzero value
    double lbval = -1;

    CCrandstate rstate;
    int seed = rand();
    CCutil_sprand(seed, &rstate); //Initialize the portable random number generator
    int ncount = m_size_subgraph; //Number of nodes (cities)
    int ecount = (ncount * (ncount - 1)) / 2; //Number of edges
    int *elist = new int[ecount * 2]; //Array giving the ends of the edges (in pairs)
    int *elen = new int[ecount]; //Array giving the weights of the edges
    int edge = 0;
    int edgeWeight = 0;
    for (int i = 0; i < ncount; i++) 
    {
        for (int j = i + 1; j < ncount; j++) 
        {
            if (i != j) 
            {
                elist[edge] = i;
                elist[edge + 1] = j;
                elen[edgeWeight] = int_dist_matrix[i][j];
                if (elen[edgeWeight] <= 0)
                    elen[edgeWeight] = 1;
                edgeWeight++;
                edge = edge + 2;
            }
        }
    }
    out_tour = CC_SAFE_MALLOC (ncount, int);
    name = CCtsp_problabel("_");
    CCdatagroup dat;
    // Populate initial solution if available
    // if (!m_init_sequence.empty())
    // {
    //     int temp_seq_array[m_init_sequence.size()];
    //     for (int i = 0; i < m_init_sequence.size(); i++)
    //     {
    //         temp_seq_array[i] = m_init_sequence[i];
    //     }
    //     in_tour = temp_seq_array;
    // }
    //Initialize a CCdatagroup
    CCutil_init_datagroup (&dat);
    //Convert a matrix of edge lengths to a CCdatagroup
    rval = CCutil_graph2dat_matrix (ncount, ecount, elist, elen, 0, &dat);
    //Solves the TSP over the graph specified in the datagroup
    try
    {
        rval = CCtsp_solve_dat (ncount, &dat, in_tour, out_tour, in_val, &optval, &success, &foundtour, name, timebound, &hit_timebound, silent, &rstate);
        // rval = CCtsp_solve_sparse (ncount, ecount, elist, elen, in_tour, out_tour, in_val, &optval, &success, &foundtour, name, timebound, &hit_timebound, silent, &rstate);
    }
    catch (int n)
    {
        cout << "Contour exception." << endl;
    }
    m_sequence.clear();
    if (success == 1)
    {
        m_solution = static_cast<double>(optval - ncount * m_dist_addon) / m_int_factor;
        for (int i = 0; i < m_size_subgraph; i++)
        {
            m_sequence.push_back(out_tour[i]);
        }
    }
    
    //szeit = CCutil_zeit();
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CC_IFFREE (out_tour, int);
    CC_IFFREE (name, char);
}

// void subgraph_tsp_solver::solve_concorde_lk()
// {
//     if (m_size_subgraph <= 4) return;

//     vector<vector<int>> int_dist_matrix;
//     int_dist_matrix.resize(m_size_subgraph);
//     for (int i = 0; i < m_size_subgraph; i++)
//     {
//         int_dist_matrix[i].resize(m_size_subgraph, 0);
//     }

//     for (int i = 0; i < m_size_subgraph; i++)
//     {
//         for (int j = i+1; j < m_size_subgraph; j++)
//         {
//             int_dist_matrix[i][j] = m_subgraph->get_relaxed_dist(i,j) * m_int_factor;
//             if (int_dist_matrix[i][j] < 0)
//                 int_dist_matrix[i][j] = 0;
//             int_dist_matrix[i][j] += m_dist_addon;
//             int_dist_matrix[j][i] = int_dist_matrix[i][j];
//         }
//     }

//     //TSP for more than 4 cities
//     int rval = 0; //Concorde functions return 1 if something fails
//     double tour_val = 0;

//     CCrandstate rstate;
//     int seed = rand();
//     CCutil_sprand(seed, &rstate); //Initialize the portable random number generator
//     int ncount = int_dist_matrix.size(); //Number of nodes (cities)
//     int ecount = (ncount * (ncount - 1)) / 2; //Number of edges
//     int *elist = new int[ecount * 2]; //Array giving the ends of the edges (in pairs)
//     int *elen = new int[ecount]; //Array giving the weights of the edges
//     int edge = 0;
//     int edgeWeight = 0;
//     for (int i = 0; i < ncount; i++) 
//     {
//         for (int j = i + 1; j < ncount; j++) 
//         {
//             if (i != j) 
//             {
//                 elist[edge] = i;
//                 elist[edge + 1] = j;
//                 elen[edgeWeight] = int_dist_matrix[i][j];
//                 if (elen[edgeWeight] <= 0)
//                     elen[edgeWeight] = 1;
//                 edgeWeight++;
//                 edge = edge + 2;
//             }
//         }
//     }
//     CCdatagroup dat;

//     //Initialize a CCdatagroup
//     CCutil_init_datagroup (&dat);

//     //Convert a matrix of edge lengths to a CCdatagroup
//     rval = CCutil_graph2dat_matrix (ncount, ecount, elist, elen, 0, &dat);

//     //Solves the TSP over the graph specified in the datagroup
//     try
//     {
//         rval = CCtsp_solve_dat_lktour(ncount, &dat, &tour_val, &rstate);
//     }
//     catch (int n)
//     {
//         cout << "Contour exception." << endl;
//     }

//     m_solution = (tour_val - ncount * m_dist_addon) / m_int_factor;

//     //szeit = CCutil_zeit();
//     CC_IFFREE (elist, int);
//     CC_IFFREE (elen, int);
// }

void subgraph_tsp_solver::get_tour_sequence(vector<int>& sequence)
{
    sequence = m_sequence;
}

void subgraph_tsp_solver::set_init_solution(vector<int>& init_seq)
{
    m_init_sequence = init_seq;
}