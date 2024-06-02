#include "Data.h"

using namespace std;

//constructor initializes sizeGraphMatrix, numberEdges, pNumber, pathMatrix and graphMatrix
Data::Data ( char * instanceName, char * options, double overlapFactor, int argc, char** argv ):ub(INFINITY)
{	

   if ( strcmp( argv[ 2 ], "2D" ) && strcmp( argv[ 2 ], "3D" ) ){
      cout << "There's no option named: " << argv[ 2 ] << endl;
      cout << "CHOOSE 2D OR 3D" << endl;
      exit( 1 );
   }

   //======= WENDA CHANGE =======
   //if ( strcmp( argv[ 3 ], "0.1" ) && strcmp( argv[ 3 ], "0.25" ) && strcmp( argv[ 3 ], "0.5" ) && strcmp( argv[ 3 ], "1.0" ) && strcmp( argv[ 3 ], "1.5" ) ){
   //   cout << "Wrong overlap factor: " << argv[ 3 ] << endl;
   //   cout << "CHOOSE BETWEEN {0.1, 0.25, 0.5, 1.0, 1.5}" << endl;
   //   exit( 1 );
   //}
   //======= END WENDA CHANGE =======

   if ( strcmp( argv[ 5 ], "V1" ) && strcmp( argv[ 5 ], "SB" ) ){
      cout << "There's no option named: " << argv[ 5 ] << endl;
      cout << "CHOOSE BETWEEN:" << endl;
      cout << "Branching Rule 1.0 = [V1]" << endl;
      cout << "Strong Branching = [SB]" << endl;
      exit( 1 );
   }

   if ( strcmp( argv[ 6 ], "DFS" ) && strcmp( argv[ 6 ], "BFS" ) && strcmp( argv[ 6 ], "BeFS" ) && strcmp( argv[ 6 ], "CBFS") ){
      cout << "There's no option named: " << argv[ 6 ] << endl;
      cout << "CHOOSE BETWEEN:" << endl;
      cout << "Depth First Search = [DFS]" << endl;
      cout << "Breadth First Search = [BFS]" << endl;
      cout << "Best First Search = [BeFS]" << endl;
      cout << "Cyclic Best First Search = [CBFS]" << endl;
      exit( 1 );
   }

   if ( strcmp( argv[ 7 ], "1" ) && strcmp( argv[ 7 ], "2" ) && strcmp( argv[ 7 ], "3" ) ){
      cout << "There's no " << argv[ 7 ] << " option" << endl;
      cout << "CHOOSE BETWEEN:" << endl;
      cout << "First Root Selection = [1]" << endl;
      cout << "Second Root Selection = [2]" << endl;
      cout << "Third Root Selection = [3]" << endl;
      exit( 1 );
   }

   overlapRatio = overlapFactor;
   ifstream readFile( instanceName, ios::in );

   // Load data from file
   if ( !readFile ){
      cerr << "Problem for openning instance file!" << endl;
      exit( 1 );
   }

   // identify reader by directory name
   string str = instanceName;
   string DirectoryName;

   string::size_type found = str.find_first_of("/");

   string::size_type loc = str.find_last_of(".", str.size() );
   string::size_type loc2 = str.find_last_of("/", str.size() );

   string folder_path; folder_path.append(instanceName,loc2);
   string::size_type last_slash = folder_path.find_last_of("/");
   DirectoryName.append( instanceName,last_slash+1,folder_path.size()-last_slash-1 );

   if (loc != string::npos) {
      fileName.append(instanceName, loc2+1, loc-loc2-1 );
   }
   else {
      fileName.append(instanceName, loc2+1, str.size() );
   }

   //		######### Mennel Instances ##########
   if ( DirectoryName.compare( "2D" ) == 0 || DirectoryName.compare( "3D" ) == 0 || DirectoryName.compare( "RND" ) == 0 ){
      setAllData( readFile );
      if( DirectoryName.compare( "2D" ) == 0 ){
         setUb( "2D", fileName, options, overlapFactor );
      }
      if( DirectoryName.compare( "3D" ) == 0 ){
         setUb( "3D", fileName, options, overlapFactor );
      }
      if( DirectoryName.compare( "RND" ) == 0 ){
         setUb( "RND", fileName, options, overlapFactor );
      }
   }
   else{

      //			######### Behdani Instances ##########
      if ( DirectoryName.compare( "Behdani" ) == 0 ){
         setSizeInst( instanceName );
         setAllData2( readFile );
         setUb( "Behdani", fileName, options, overlapFactor );
      }
      //			######### Behdani Instances ##########
      else if (DirectoryName.compare( "medium_2D_Behdani_CETSPs" )==0 ){
         setSize_GeordanBehdanis( fileName );
         setAllData2( readFile );
         setUb( "medium_2D_Behdani_CETSPs", fileName, options, overlapFactor);
      }
      else{
         cout << "Something wrong with the directory names!" << endl;
         exit( 1 );
      }


   }

   if ( strcmp( options, "2D" ) == 0 ){
      deleteZCoordinates();
   }

   readFile.close();

   setDistanceMatrix();
   setDepotFarthest();

   //  compute overlap ratio for this instance
   overlap_for_this_inst = 0;
   double sum_radii = 0;
   for ( int i = 1; i < sizeInst; i++ ){
      sum_radii += radius[ i ];
   }

   m_avg_radius = sum_radii/( sizeInst - 1 );

   vector< double > x_aux = coordx;
   vector< double > y_aux = coordy;
   vector< double > z_aux = coordz;
   //		
   sort( x_aux.begin(), x_aux.end() );
   sort( y_aux.begin(), y_aux.end() );
   sort( z_aux.begin(), z_aux.end() );

   vector< double > largest_coord = { x_aux[ sizeInst - 1 ], y_aux[ sizeInst - 1 ], z_aux[ sizeInst - 1 ] };

   sort( largest_coord.begin(), largest_coord.end() );

   overlap_for_this_inst = m_avg_radius/largest_coord[ 2 ];
   //finish set overlap for this instance		
}

// destrutor
Data::~Data()
{

}

void Data::setSizeInst( char * instancePath )
{
   string path=instancePath;
   string::size_type last_slash = path.find_last_of("/");
   string::size_type first_dash=path.find_first_of("-",last_slash);
   string::size_type last_dash =path.find_last_of("-",last_slash);
   
   sizeInst=stoi(path.substr(first_dash+1,last_dash-first_dash));
   sizeInst++;
}

void Data::setSize_GeordanBehdanis( const string& filename){
   string size_string;size_string.append(filename,1,filename.find_first_of("_")-1);
   sizeInst=stoi(size_string);
}

void Data::setUb( string drctName, string instName, char * option, double overlap )
{	
   string tempString;
   cout << "Name: " << instName << endl;
   cout << "directory: " << drctName << endl;
   cout << "option: " << option << endl;

   if( drctName.compare( "2D" ) == 0 ){
      if( strcmp( option, "3D" ) == 0 ){
         cout << "Inconsistent Option For This Instances!" << endl;
         exit( 1 );
      }
      if( strcmp( option, "2D" ) == 0 ){
         if( overlap == 1.0 ){
            ifstream readFile( "Upper_Bounds/2D10.txt", ios::in );
            while ( tempString.compare( instName ) != 0 ){
               readFile >> tempString;
            }
            readFile >> tempString;
            ub = atof( tempString.c_str() );
         }
         else{
            cout << "Wrong Overlap Factor!" << endl;
            exit( 1 );
         }
      }
   }
   if( drctName.compare( "3D" ) == 0 ){
      if( strcmp( option, "2D" ) == 0 ){
         if( overlap == 0.1 ){
            ifstream readFile( "Upper_Bounds/2D01.txt", ios::in );
            while ( tempString.compare( instName ) != 0 ){
               readFile >> tempString;
            }
            readFile >> tempString;
            ub = atof( tempString.c_str() );
         }
         if( overlap == 0.5 ){
            ifstream readFile( "Upper_Bounds/2D05.txt", ios::in );
            while ( tempString.compare( instName ) != 0 ){
               readFile >> tempString;
            }
            readFile >> tempString;
            ub = atof( tempString.c_str() );
         }
         if( overlap == 1.5 ){
            ifstream readFile( "Upper_Bounds/2D15.txt", ios::in );
            while ( tempString.compare( instName ) != 0 ){
               readFile >> tempString;
            }
            readFile >> tempString;
            ub = atof( tempString.c_str() );
         }
         if( overlap != 0.1 && overlap != 0.5 && overlap != 1.5 ){
            ub = 10000001.0;
            //======= WENDA CHANGE =======
            // remove upper bound check
            //cout << "Wrong Overlap Factor!" << endl;
            //exit( 1 );
            //======= END WENDA CHANGE =======
         }
      }
      if( strcmp( option, "3D" ) == 0 ){
         if( overlap == 0.1 ){
            ifstream readFile( "Upper_Bounds/3D01.txt", ios::in );
            while ( tempString.compare( instName ) != 0 ){
               readFile >> tempString;
            }
            readFile >> tempString;
            ub = atof( tempString.c_str() );
         }
         if( overlap == 0.5 ){
            ifstream readFile( "Upper_Bounds/3D05.txt", ios::in );
            while ( tempString.compare( instName ) != 0 ){
               readFile >> tempString;
            }
            readFile >> tempString;
            ub = atof( tempString.c_str() );
         }
         if( overlap == 1.5 ){
            ifstream readFile( "Upper_Bounds/3D15.txt", ios::in );
            while ( tempString.compare( instName ) != 0 ){
               readFile >> tempString;
            }
            readFile >> tempString;
            ub = atof( tempString.c_str() );
         }
         if( overlap != 0.1 && overlap != 0.5 && overlap != 1.5 ){
            ub = 10000001.0;
            //======= WENDA CHANGE =======
            // remove upper bound check
            //cout << "Wrong Overlap Factor!" << endl;
            //exit( 1 );
            //======= END WENDA CHANGE =======
         }
      }
   }

   if( drctName.compare( "RND" ) == 0 ){
      if( strcmp( option, "2D" ) == 0 ){
         if( overlap == 1.0 ){
            ifstream readFile( "Upper_Bounds/2D10.txt", ios::in );
            while ( tempString.compare( instName ) != 0 ){
               readFile >> tempString;
            }
            readFile >> tempString;
            ub = atof( tempString.c_str() );
         }
         else{
            cout << "Wrong Overlap Factor!" << endl;
            exit( 1 );
         }
      }
      if( strcmp( option, "3D" ) == 0 ){
         if( overlap == 1.0 ){
            ifstream readFile( "Upper_Bounds/3D10.txt", ios::in );
            while ( tempString.compare( instName ) != 0 ){
               readFile >> tempString;
            }
            readFile >> tempString;
            ub = atof( tempString.c_str() );
         }
         else{
            cout << "Wrong Overlap Factor!" << endl;
            exit( 1 );
         }
      }
   }

   if( drctName.compare( "Behdani" ) == 0 ){
      if( overlap != 1.0 && overlap != 0.5 && overlap != 0.25 ){
         cout << "Wrong Overlap Factor!" << endl;
         exit( 1 );
      }
      if( strcmp( option, "2D" ) == 0 ){
         ub = 999;
      }
      if( strcmp( option, "3D" ) == 0 ){
         cout << "Inconsistent Option For This Instances!" << endl;
         exit( 1 );
      }			
   }
}

double Data::getUbFromFile()
{
   return ub;
}

double Data::getOverlapRatio()
{
   return overlap_for_this_inst;
}

//	====== SET Functions ======
void Data::setAllData( ifstream &readFile )
{
   sizeInst = 0;	
   string trashString;		

   while ( trashString.compare("//Depot") != 0 && trashString.compare("//Depot:") != 0 ) {
      readFile >> trashString;
   }

   //########### read depot coordinates ############

   //get the "is" out
   if( trashString.compare("//Depot:") != 0 ){
      readFile >> trashString;
   }

   double temp[ 3 ];	
   readFile >> temp[ 0 ] >> trashString >> temp[ 1 ] >> trashString >> temp[ 2 ];
   coordx.push_back(temp[ 0 ]);
   coordy.push_back(temp[ 1 ]);
   coordz.push_back(temp[ 2 ]);
   radius.push_back( 0 );
   demand.push_back( 0 );	
   sizeInst++;

   readFile.seekg( 0 );

   //##########  read customers coordinates	  #################
   while ( trashString.compare("description") != 0 ) {
      readFile >> trashString;
   }

   readFile >> trashString;

   while ( trashString.compare("//Depot") != 0 && trashString.compare("//Depot:") != 0 ) {

      double temp = atof(trashString.c_str());
      coordx.push_back(temp);
      readFile >> temp;
      coordy.push_back(temp);
      readFile >> temp;
      coordz.push_back(temp);
      readFile >> temp; 
      radius.push_back(temp);
      readFile >> temp;
      demand.push_back(temp);
      sizeInst++;
      readFile >> trashString;
   }

   //multiply ratio by the overlapFactor
   for ( int i = 0; i < sizeInst; i++ ){
      radius[ i ] = radius[ i ]*overlapRatio;
   }	
}

void Data::setAllData2( ifstream &readFile )
{
   double x = 0;
   double y = 0;
   double z = 0;

   for ( int i = 0; i < sizeInst; i++ ){
      readFile >> x >> y;
      coordx.push_back( x );
      coordy.push_back( y );
      coordz.push_back( 0 );
      radius.push_back( 1 );
      demand.push_back( 0 );
   }

   radius[ 0 ] = 0;

   for ( int i = 0; i < sizeInst; i++ ){
      radius[ i ] = radius[ i ]*overlapRatio;
   }
}

void Data::deleteZCoordinates()
{
   for ( int i = 0; i < sizeInst; i++ ){
      coordz[ i ] = 0;
   }
   cout << "OPTION 2D ACTIVATED, Z COORDINATES BECAME 0" << endl;
}

void Data::setIntersecMatrix()
{
   long double test = 0;

   intersecMatrix.resize(sizeInst);
   for ( int i = 0; i < sizeInst; i++ ){
      intersecMatrix[i].resize(sizeInst, 0);
   }

   for ( int i = 0; i < sizeInst; i++ ){
      for ( int j = i + 1; j < sizeInst; j++ ){
         test = sqrt( pow( coordx[ i ] - coordx[ j ], 2 ) + pow( coordy[ i ] - coordy[ j ], 2 ) + pow( coordz[ i ] - coordz[ j ], 2 ) ) - ( radius[ i ] + radius[ j ] );
         if ( test <= 0 ){
            intersecMatrix[ i ][ j ] = 1;
            intersecMatrix[ j ][ i ] = 1;
         }
         else{
            intersecMatrix[ i ][ j ] = 0;
            intersecMatrix[ j ][ i ] = 0;
         } 
      }
   }

   for ( int i = 0; i < sizeInst; i++ ){
      for ( int j = 0; j < sizeInst; j++ ){
         countIntersecs.push_back( 0 );
         countIntersecs[ i ] += intersecMatrix[ i ][ j ];
      }
   }
}

bool compare ( pStruct p1, pStruct p2 )
{
   for ( int i = 0; i < p2.size && i != p2.index && p1.size != 0 && p2.size != 0; i++ ){
      if ( p2.intersections[ i ] == p1.index ){
         return true;
      }
   }

   if ( p1.nOfIntersec > p2.nOfIntersec ){
      return true;
   }
   if ( p1.nOfIntersec < p2.nOfIntersec ){
      return false;
   }

   if ( p1.nOfIntersec == p2.nOfIntersec ){
      if ( p1.radius > p1.radius ){
         return true;
      }
      if ( p1.radius <= p1.radius ){
         return false;
      }
   }				
}

void Data::setBranchingRuleList()
{
   pStruct p;
   list < pStruct > priorityList;
   list< pStruct >::iterator itPList;
   list< int >::iterator itBList;

   for ( int i = 0; i < sizeInst; i++ ){
      p.index = i;
      p.nOfIntersec = countIntersecs[ i ];
      p.radius = radius[ i ];
      p.intersections.clear();
      for ( int j = 0; j < sizeInst; j++ ){
         if ( intersecMatrix[ i ][ j ] == 1 ){
            p.intersections.push_back( j );
         }
      }
      p.size = p.intersections.size();

      priorityList.push_back( p );

   }

   itPList = priorityList.begin();
   int k = 8;
   priorityList.sort( compare );
   priorityList.splice( priorityList.begin(), priorityList, itPList);

   for ( itPList = priorityList.begin(); itPList != priorityList.end(); itPList++ ){
      branchingRuleList.push_back( (*itPList).index );
   }
}

void Data::printInstance()
{
   for ( int i = 0; i < sizeInst; i++ ){
      cout << coordx[ i ] << " " << coordy[ i ] << " " << coordz[ i ] << " " << radius[ i ] << " " << demand[ i ]<< endl;
   }
}

void Data::printSizes()
{
   cout<<"N Radii: "<<radius.size()<<" Nx: "<<coordx.size()<<" Ny: "<<coordy.size()<<" Nz: "<<coordz.size()<<" Nd: "<<demand.size()<<endl;
}

int Data::getBranchingRuleList( int i )
{
   iterBranchingRuleList = branchingRuleList.begin();

   for ( int j = 0; j < i; j++ ){
      iterBranchingRuleList++;
   }

   return ( *( iterBranchingRuleList ) );
}

double euclidianDistance( coordinates * p1, coordinates * p2 )
{
   return sqrt( pow( p2->x - p1->x, 2) + pow( p2->y - p1->y, 2) + pow( p2->z - p1->z, 2) );
}

void Data::setDistanceMatrix()
{
   distanceMatrix.resize( sizeInst );
   relaxed_distance_matrix.resize(sizeInst);

   for ( int i = 0; i < sizeInst; i++ ){
      distanceMatrix[ i ].resize( sizeInst );
      relaxed_distance_matrix[i].resize(sizeInst, 0);
   }

   coordinates * p1 = new coordinates;
   coordinates * p2 = new coordinates;

   for ( int i = 0; i < sizeInst; i++ ){
      for ( int j = 0; j < sizeInst; j++ ){
         p1->x = coordx[ i ];
         p1->y = coordy[ i ];
         p1->z = coordz[ i ];

         p2->x = coordx[ j ];
         p2->y = coordy[ j ];
         p2->z = coordz[ j ];				

         distanceMatrix[ i ][ j ] = euclidianDistance( p1, p2 );
         relaxed_distance_matrix[i][j] = distanceMatrix[i][j];
         relaxed_distance_matrix[i][j] = max(relaxed_distance_matrix[i][j] - radius[i] - radius[j], 0.0);
      }
   }

   delete p1;
   delete p2;

}

bool simpleSortVector ( pair< int, double > i, pair< int, double > j) { return ( i.second > j.second ); }

void Data::setDepotFarthest()
{
   pair< int, double > temp;

   for ( int i = 0; i < sizeInst; i++ ){
      temp = make_pair( i, distanceMatrix[ 0 ][ i ] );
      farthest.push_back( temp );
   }

   sort ( farthest.begin(), farthest.end(), simpleSortVector );
}

int Data::getDepotFarthest( int i )
{
   return farthest[ i ].first;
}

int Data::getSizeInst()
{
   return sizeInst;	
}

double Data::getCoordx( int i )
{
   return coordx[ i ];	
}

double Data::getCoordy( int i )
{
   return coordy[ i ];	
}

double Data::getCoordz( int i )
{
   return coordz[ i ];	
}

double Data::getRadius( int i )
{
   return radius[ i ];	
}

double Data::getDemand( int i )
{
   return demand[ i ];	
}

double Data::get_relaxed_dist(int i, int j)
{
   return relaxed_distance_matrix[i][j];
}

void Data::print_relaxed_dists()
{
   ofstream relaxed_dist_out;
   relaxed_dist_out.open("rexdist.dat");
   relaxed_dist_out << "[";
   for (int i = 0; i < sizeInst; i++)
   {
      if (i != 0)
      {
         relaxed_dist_out << "," << endl;
      }
      relaxed_dist_out << "[";
      for (int j = 0; j < sizeInst; j++)
      {
         if (j != 0)
         {
            relaxed_dist_out << ", ";
         }
         if (i == j)
         {
            relaxed_dist_out << "9999";
         }
         else
         {
            relaxed_dist_out << relaxed_distance_matrix[i][j];
         }
      }
      relaxed_dist_out << "]";
   }
   relaxed_dist_out << "]" << endl;
}

bool Data::is_intersect(int i, int j)
{
   if (distanceMatrix[i][j] - radius[i] - radius[j] <= 0)
      return true;
   return false;
}