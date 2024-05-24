#include "SolveSocpCplex.h"

ILOSTLBEGIN // import namespace std

SolveSocpCplex::SolveSocpCplex( Data *dataObject, vector< int >& sequence )
: objectData( dataObject ), model( env ), SOCP( model ), xCoord( env ), yCoord( env ), zCoord( env ), x(env, sequence.size(), -IloInfinity, IloInfinity), y(env, sequence.size(), -IloInfinity, IloInfinity), z(env, sequence.size(), -IloInfinity, IloInfinity)
{			
   sizeProblem = sequence.size();
   f_value = 0;		
}

SolveSocpCplex::SolveSocpCplex(Data* data, int size_seq) : objectData(data), sizeProblem(size_seq), model(env), SOCP(model), xCoord(env), yCoord(env), zCoord(env)
{
   x = IloNumVarArray(env, sizeProblem, -IloInfinity, IloInfinity);
   y = IloNumVarArray(env, sizeProblem, -IloInfinity, IloInfinity);
   z = IloNumVarArray(env, sizeProblem, -IloInfinity, IloInfinity);
   m_allf = IloNumVarArray(env, sizeProblem, 0, IloInfinity);
   m_allw = IloNumVarArray(env, sizeProblem, -IloInfinity, IloInfinity);
   m_allu = IloNumVarArray(env, sizeProblem, -IloInfinity, IloInfinity);
   m_allv = IloNumVarArray(env, sizeProblem, -IloInfinity, IloInfinity);
   m_alls = IloNumVarArray(env, sizeProblem, -IloInfinity, IloInfinity);
   m_allt = IloNumVarArray(env, sizeProblem, -IloInfinity, IloInfinity);
   m_allq = IloNumVarArray(env, sizeProblem, -IloInfinity, IloInfinity);

   m_dists_to_c = IloRangeArray(env);
   m_coord_x = IloRangeArray(env);
   m_coord_y = IloRangeArray(env);
   m_coord_z = IloRangeArray(env);

   // initialize_init_point();

   f_value = 0;
}

// destrutor
SolveSocpCplex::~SolveSocpCplex()
{
   env.end();
}

int SolveSocpCplex::setSizeProblem( vector< int >& sequence )
{
   sizeProblem = sequence.size();
   return sizeProblem;
}

void SolveSocpCplex::setF_value()
{
   f_value = SOCP.getObjValue();
}

double SolveSocpCplex::getF_value()
{
   return f_value;
}
void SolveSocpCplex::printF_value()
{
   cout << "Value of objective function: "<< f_value << endl;
}
void SolveSocpCplex::printSolution( vector<int>& sequence )
{
   vector< IloNumArray > solution;

   solution.push_back( xCoord );
   solution.push_back( yCoord );
   solution.push_back( zCoord );

   cout << "Solution: " << endl;

   for ( int j = 0; j < sizeProblem; j++ ){
      cout << "( "<< solution[ 0 ][ j ] << ", " << solution[ 1 ][ j ] << ", " << solution[ 2 ][ j ] << " )" << endl;
   }

   cout<< endl;
}

double SolveSocpCplex::getSolutionX( int i )
{
   return xCoord[ i ];
}

double SolveSocpCplex::getSolutionY( int i )
{
   return yCoord[ i ];
}

double SolveSocpCplex::getSolutionZ( int i )
{
   return zCoord[ i ];
}

void SolveSocpCplex::solveSOCP( vector< int >& sequence )
{
   cout<<"Solve w/ CPLEX"<<endl;
   static pthread_mutex_t cs_mutex = PTHREAD_MUTEX_INITIALIZER;

   createModel( sequence );
   SOCP.setOut( env.getNullStream() );
   SOCP.setParam( IloCplex::BarQCPEpComp,  1e-12 );
   SOCP.setParam( IloCplex::Threads, 1 );
   SOCP.setParam( IloCplex::ParallelMode, 1 );
   SOCP.solve();
   SOCP.getStatus();
   violation = SOCP.getQuality(IloCplex::SumScaledPrimalInfeas);
   setF_value();
   SOCP.getValues( xCoord, x );		
   SOCP.getValues( yCoord, y );
   SOCP.getValues( zCoord, z );

   pthread_mutex_unlock(&cs_mutex);
}

void SolveSocpCplex::finishSOCP()
{
   SOCP.end();
}

void SolveSocpCplex::createModel( vector< int >& sequence )
{
   //variables
   //add variable f ( head of cones )
   //================================================================
   IloNumVarArray f(env, sizeProblem, 0, IloInfinity);
   char var[100];
   for (int i = 0; i < sizeProblem; i++) {
      sprintf(var, "f(%d)", i );
      f[i].setName(var);
      model.add(f[ i ]);
   }
   //getchar();
   //add variables corresponding to coordinates x, y and z
   //================================================================
   //variable x
   for (int i = 0; i < sizeProblem; i++) {
      sprintf( var, "x(%d)", i );
      x[i].setName( var );
      model.add( x[ i ] );
   }

   //variable y
   for (int i = 0; i < sizeProblem; i++) {
      sprintf(var, "y(%d)", i );
      y[i].setName(var);
      model.add( y[ i ] );
   }

   //variable z
   for (int i = 0; i < sizeProblem; i++) {
      sprintf(var, "z(%d)", i );
      z[i].setName(var);
      model.add( z[ i ] );
   }

   //add auxiliary variables w, u, v, s, t and q
   //================================================================
   IloNumVarArray w(env, sizeProblem, -IloInfinity, IloInfinity);
   for (int i = 0; i < sizeProblem; i++) {
      sprintf(var, "w(%d)", i );
      w[i].setName(var);
      model.add(w[ i ]);
   }

   IloNumVarArray u(env, sizeProblem, -IloInfinity, IloInfinity);
   for (int i = 0; i < sizeProblem; i++) {
      sprintf(var, "u(%d)", i );
      u[i].setName(var);
      model.add(u[ i ]);
   }

   IloNumVarArray v(env, sizeProblem, -IloInfinity, IloInfinity);
   for (int i = 0; i < sizeProblem; i++) {
      sprintf(var, "v(%d)", i );
      v[i].setName(var);
      model.add(v[ i ]);
   }

   IloNumVarArray s(env, sizeProblem, -IloInfinity, IloInfinity);
   for (int i = 0; i < sizeProblem; i++) {
      sprintf(var, "s(%d)", i );
      s[i].setName(var);
      model.add(s[ i ]);
   }

   IloNumVarArray t(env, sizeProblem, -IloInfinity, IloInfinity);
   for (int i = 0; i < sizeProblem; i++) {
      sprintf(var, "t(%d)", i );
      t[i].setName(var);
      model.add(t[ i ]);
   }

   IloNumVarArray q(env, sizeProblem, -IloInfinity, IloInfinity);
   for (int i = 0; i < sizeProblem; i++) {
      sprintf(var, "q(%d)", i );
      q[i].setName(var);
      model.add(q[ i ]);
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
      IloRange r = ( s[ i ]*s[ i ] + t[ i ]*t[ i ] + q[ i ]*q[ i ] <= objectData->getRadius( sequence[ i ] )*objectData->getRadius( sequence[ i ] ) );
      char c[100];
      sprintf(c, "q%d", i );
      r.setName(c);
      model.add(r);
   }

   //radius constraints	
   //=============================================================
   for (int i = 0; i < sizeProblem; i++) {
      IloRange r = ( s[ i ] + x[ i ] == objectData->getCoordx( sequence[ i ] ) );
      char c[100];
      sprintf(c, "b%d", i );
      r.setName(c);
      model.add(r);
   }

   for (int i = 0; i < sizeProblem; i++) {
      IloRange r = ( t[ i ] + y[ i ] == objectData->getCoordy( sequence[ i ] ) );
      char c[100];
      sprintf(c, "b%d", i + sizeProblem);
      r.setName(c);
      model.add(r);
   }

   for (int i = 0; i < sizeProblem; i++) {
      IloRange r = ( q[ i ] + z[ i ] == objectData->getCoordz( sequence[ i ] ) );
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
}

void SolveSocpCplex::initialize_model()
{
   model.add(m_allf);
   model.add(x); model.add(y); model.add(z);
   model.add(m_allw); model.add(m_allu); model.add(m_allv);
   model.add(m_allq); model.add(m_alls); model.add(m_allt);

   IloExpr FO(env);
   for (int i = 0; i < sizeProblem; i++)
   {
      FO += m_allf[i];
   }
   model.add(IloMinimize(env, FO));
   //set SOC constraints
   //===============================================
   for (int j = 0; j < sizeProblem; j++) 
   {
      IloRange r = ( - m_allf[ j ]*m_allf[ j ] + m_allw[ j ]*m_allw[ j ] + m_allu[ j ]*m_allu[ j ] + m_allv[ j ]*m_allv[ j ] <= 0);
      // char c[100];
      // sprintf(c, "cone%d", j );
      // r.setName(c);
      model.add(r);
   }
   //separating first constraint
   IloRange isol_1 = ( m_allw[ 0 ] - x[ sizeProblem - 1 ] + x[ 0 ] == 0 );
   // char c[100];
   // sprintf(c, "c%d", 0);
   // isol_1.setName(c);
   model.add( isol_1 );

   for (int j = 1; j < sizeProblem; j++) 
   {
      IloRange r = ( m_allw[ j ] - x[ j - 1 ] + x[ j ] == 0 );
      // char c[100];
      // sprintf(c, "c%d", j);
      // r.setName(c);
      model.add(r);
   }

   IloRange isol_2 = ( m_allu[ 0 ] - y[ sizeProblem - 1 ] + y[ 0 ] == 0 );
   // char c2[100];
   // sprintf(c2, "c%d", sizeProblem );
   // isol_2.setName(c2);
   model.add( isol_2 );

   for (int j = 1; j < sizeProblem; j++) 
   {
      IloRange r = ( m_allu[ j ] - y[ j - 1 ] + y[ j ] == 0 );
      // char c[100];
      // sprintf(c, "c%d", j + sizeProblem);
      // r.setName(c);
      model.add(r);
   }

   IloRange isol_3 = ( m_allv[ 0 ] - z[ sizeProblem - 1 ] + z[ 0 ] == 0 );
   // char c3[100];
   // sprintf(c3, "c%d", sizeProblem );
   // isol_3.setName(c3);
   model.add( isol_3 );

   for (int j = 1; j < sizeProblem; j++) 
   {
      IloRange r = ( m_allv[ j ] - z[ j - 1 ] + z[ j ] == 0 );
      // char c[100];
      // sprintf(c, "c%d", j + 2*sizeProblem);
      // r.setName(c);
      model.add(r);
   }
 
   SOCP.setOut( env.getNullStream() );
   SOCP.setParam( IloCplex::BarQCPEpComp,  1e-12 );
   SOCP.setParam( IloCplex::Threads, 1 );
   SOCP.setParam( IloCplex::ParallelMode, 1 );
}

void SolveSocpCplex::populate_removable_constraints(vector<int>& sequence)
{
   m_dists_to_c.clear();
   m_coord_x.clear();
   m_coord_y.clear();
   m_coord_z.clear();
   for (int i = 0; i < sizeProblem; i++) {
      IloRange r = (m_alls[i]*m_alls[i] + m_allt[i]*m_allt[i] + m_allq[i]*m_allq[i] <= objectData->getRadius( sequence[ i ] )*objectData->getRadius( sequence[ i ] ) );
      // char c[100];
      // sprintf(c, "q%d", i );
      // r.setName(c);
      m_dists_to_c.add(r);
      // model.add(r);
   }

   //radius constraints	
   //=============================================================
   for (int i = 0; i < sizeProblem; i++) {
      IloRange r = (m_alls[i] + x[i] == objectData->getCoordx(sequence[i]));
      // char c[100];
      // sprintf(c, "b%d", i );
      // r.setName(c);
      m_coord_x.add(r);
      // model.add(r);
   }

   for (int i = 0; i < sizeProblem; i++) {
      IloRange r = (m_allt[i] + y[i] == objectData->getCoordy(sequence[i]));
      // char c[100];
      // sprintf(c, "b%d", i + sizeProblem);
      // r.setName(c);
      m_coord_y.add(r);
      // model.add(r);
   }

   for (int i = 0; i < sizeProblem; i++) {
      IloRange r = (m_allq[i] + z[i] == objectData->getCoordz(sequence[i]));
      // char c[100];
      // sprintf(c, "b%d", i + 2*sizeProblem);
      // r.setName(c);
      m_coord_z.add(r);
      // model.add(r);
   }

   model.add(m_dists_to_c);
   model.add(m_coord_x);
   model.add(m_coord_y);
   model.add(m_coord_z);
}

void SolveSocpCplex::populate_removable_constraints(vector<int>& sequence, int prev_pos, int curr_pos)
{
   IloRange temp_const;
   for (int i = prev_pos; i <= curr_pos; i++)
   {
      temp_const = (m_alls[i]*m_alls[i] + m_allt[i]*m_allt[i] + m_allq[i]*m_allq[i] <= objectData->getRadius(sequence[i])*objectData->getRadius(sequence[i]));
      m_dists_to_c[i] = temp_const;
      temp_const = (m_alls[i] + x[i] == objectData->getCoordx(sequence[i]));
      m_coord_x[i] = temp_const;
      temp_const = (m_allt[i] + y[i] == objectData->getCoordy(sequence[i]));
      m_coord_y[i] = temp_const;
      temp_const = (m_allq[i] + z[i] == objectData->getCoordz(sequence[i]));
      m_coord_z[i] = temp_const;

      model.add(m_dists_to_c[i]);
      model.add(m_coord_x[i]);
      model.add(m_coord_y[i]);
      model.add(m_coord_z[i]);
   }
}

void SolveSocpCplex::clear_removable_constraints()
{
   model.remove(m_dists_to_c);
   model.remove(m_coord_x);
   model.remove(m_coord_y);
   model.remove(m_coord_z);
}

void SolveSocpCplex::clear_removable_constraints(int prev_pos, int curr_pos)
{
   for (int i = prev_pos; i <= curr_pos; i++)
   {
      model.remove(m_dists_to_c[i]);
      model.remove(m_coord_x[i]);
      model.remove(m_coord_y[i]);
      model.remove(m_coord_z[i]);
   }
}

void SolveSocpCplex::solveSOCP()
{
   cout<<"Solve w/ CPLEX no argument"<<endl;
   static pthread_mutex_t cs_mutex = PTHREAD_MUTEX_INITIALIZER;
   // if (m_num_solves != 0)
   // {
   //    SOCP.setStart(m_allvarvals, m_allvar_redcosts, m_allvars, 0, 0, 0);
   // }
   SOCP.solve();
   SOCP.getStatus();
   violation = SOCP.getQuality(IloCplex::SumScaledPrimalInfeas);
   f_value = SOCP.getObjValue();
   SOCP.getValues( xCoord, x );		
   SOCP.getValues( yCoord, y );
   SOCP.getValues( zCoord, z );
   // set_init_point_values();
   m_num_solves++;

   pthread_mutex_unlock(&cs_mutex);
}

void SolveSocpCplex::initialize_init_point()
{
   m_allvars = IloNumVarArray(env);
   m_allvarvals = IloNumArray(env);
   m_allvar_redcosts = IloNumArray(env);
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvars.add(x[i]);
   }
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvars.add(y[i]);
   }
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvars.add(z[i]);
   }
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvars.add(m_allf[i]);
   }
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvars.add(m_allw[i]);
   }
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvars.add(m_allu[i]);
   }
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvars.add(m_allv[i]);
   }
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvars.add(m_alls[i]);
   }
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvars.add(m_allt[i]);
   }
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvars.add(m_allq[i]);
   }
}

void SolveSocpCplex::set_init_point_values()
{
   m_allvarvals.clear();
   m_allvar_redcosts.clear();

   IloNumArray temp_values(env);
   IloNumArray temp_redcosts(env);

   SOCP.getReducedCosts(temp_redcosts, x);
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvarvals.add(xCoord[i]);
      m_allvar_redcosts.add(temp_redcosts[i]);
   }

   temp_redcosts.clear();
   SOCP.getReducedCosts(temp_redcosts, y);
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvarvals.add(yCoord[i]);
      m_allvar_redcosts.add(temp_redcosts[i]);
   }

   temp_redcosts.clear();
   SOCP.getReducedCosts(temp_redcosts, z);
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvarvals.add(zCoord[i]);
      m_allvar_redcosts.add(temp_redcosts[i]);
   }

   SOCP.getValues(temp_values, m_allf);
   temp_redcosts.clear();
   SOCP.getReducedCosts(temp_redcosts, m_allf);
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvarvals.add(temp_values[i]);
      m_allvar_redcosts.add(temp_redcosts[i]);
   }

   temp_values.clear();
   SOCP.getValues(temp_values, m_allw);
   temp_redcosts.clear();
   SOCP.getReducedCosts(temp_redcosts, m_allw);
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvarvals.add(temp_values[i]);
      m_allvar_redcosts.add(temp_redcosts[i]);
   }

   temp_values.clear();
   SOCP.getValues(temp_values, m_allu);
   temp_redcosts.clear();
   SOCP.getReducedCosts(temp_redcosts, m_allu);
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvarvals.add(temp_values[i]);
      m_allvar_redcosts.add(temp_redcosts[i]);
   }

   temp_values.clear();
   SOCP.getValues(temp_values, m_allv);
   temp_redcosts.clear();
   SOCP.getReducedCosts(temp_redcosts, m_allv);
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvarvals.add(temp_values[i]);
      m_allvar_redcosts.add(temp_redcosts[i]);
   }

   temp_values.clear();
   SOCP.getValues(temp_values, m_alls);
   temp_redcosts.clear();
   SOCP.getReducedCosts(temp_redcosts, m_alls);
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvarvals.add(temp_values[i]);
      m_allvar_redcosts.add(temp_redcosts[i]);
   }

   temp_values.clear();
   SOCP.getValues(temp_values, m_allt);
   temp_redcosts.clear();
   SOCP.getReducedCosts(temp_redcosts, m_allt);
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvarvals.add(temp_values[i]);
      m_allvar_redcosts.add(temp_redcosts[i]);
   }

   temp_values.clear();
   SOCP.getValues(temp_values, m_allq);
   temp_redcosts.clear();
   SOCP.getReducedCosts(temp_redcosts, m_allq);
   for (int i = 0; i < sizeProblem; i++)
   {
      m_allvarvals.add(temp_values[i]);
      m_allvar_redcosts.add(temp_redcosts[i]);
   }
}