#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <cstdlib>
#include <stdio.h>
#include <cfloat>
#include <gmp.h>
#include <gmpxx.h>
#include <math.h>

#include"PrintFunctions.h"
#include"util.h"
#include "Data.h"

void printDataToFile ( Data * dataptr, char* option, double overlap )
{
   ofstream outFile ( "Results/results_BnB_CETSP_OPTF.txt", ios::app );

   if ( !outFile ) {
      cout << "file cannot be created\n";
      exit (1);
   }

   if ( !outFile ) {
      cout << "file cannot be created\n";
      exit (1);
   }

   outFile << fixed << setiosflags ( ios::showpoint ) << setprecision( 2 );

   outFile  << dataptr->fileName << " "
      << option << " "
      << overlap << " "
      << dataptr->getOverlapRatio() << " ";
   outFile << endl;
   outFile.close();		
}

void printDataToFile ( Data * dataptr, char* option, double overlap, int sizeInst, double bestKnown, 
                       double best, double best_lb, double gap_root, int count_SOCP_solved, int itCount, 
                       double computationTime, double sbComputationTime, double totalSocpCompTime, 
                       int itToIncum, int numLNodes, int numNodes, int branchingStrategy )
{
   ofstream outFile ( "Results/results_BnB_CETSP_OPTF.txt", ios::app );

   if ( !outFile ) {
      cout << "file cannot be created\n";
      exit (1);
   }

   outFile << fixed << setiosflags ( ios::showpoint ) << setprecision( 5 );

   outFile  << "OPTF" << " "
      << dataptr->fileName << " "
      << option << " "
      << overlap << " "
      << dataptr->getOverlapRatio() << " "
      << sizeInst << " ";

   if (bestKnown >= 10000000.0) outFile << " - ";
   else outFile << bestKnown << " ";

   outFile << best << " ";

   if(best_lb >= 999999999999999) outFile << " - ";
   else outFile << best_lb << " ";

   outFile
      << gap_root << " "
      << count_SOCP_solved << " "
      << itCount << " "
      << computationTime << " "
      << sbComputationTime << " "
      << totalSocpCompTime << " "
      << itToIncum << " "
      << numLNodes << " "
      << numNodes << " "
      << branchingStrategy << " ";
   outFile << endl;
   outFile.close();		
}

void printDataToFile( Data * dataptr, char* option, double overlap, int sizeInst, double bestKnown, 
                      double ub, double best_lb, double gap_real, double gap_lb_bnb, double gap_root, 
                      int count_SOCP_solved, int itCount, double computationTime, double sbComputationTime, 
                      double totalSocpCompTime, int itToIncum, int numLNodes, int numNodes,
                      int branchingStrategy )
{
   ofstream outFile ( "Results/results_BnB_CETSP_noOPTF.txt", ios::app );

   if ( !outFile ) {
      cout << "file cannot be created\n";
      exit (1);
   }

   outFile << fixed << setiosflags ( ios::showpoint ) << setprecision( 5 );

   outFile  << "noOPTF" << " "
      << dataptr->fileName << " "
      << option << " "
      << overlap << " "
      << dataptr->getOverlapRatio() << " "
      << sizeInst << " ";

   if (bestKnown >= 10000000.0) outFile << " - ";
   else outFile << bestKnown << " ";

   outFile << ub << " "
      << best_lb << " "
      << gap_real << " "
      << gap_lb_bnb << " "
      << gap_root << " "
      << count_SOCP_solved << " "
      << itCount << " "
      << computationTime << " "
      << sbComputationTime << " "
      << totalSocpCompTime << " "
      << itToIncum << " "
      << numLNodes << " "
      << numNodes << " "
      << branchingStrategy << " ";
   outFile << endl;
   outFile.close();		
}

void printDataToMatlab( Data * dataptr, int sizeInst, double overlap, double best, vector< int > solucao, vector< vector< double > > solucaoXYZ )
{
   char nomeArquivo[ 100 ] = {"Results/MATLAB/"};
   strcat (nomeArquivo, dataptr->fileName.c_str());
   strcat (nomeArquivo, ".txt");	
   ofstream outFile2 ( nomeArquivo );

   if ( !outFile2 ) {
      cout << "file cannot be created\n";
      exit (1);
   }

   outFile2 << fixed << setiosflags ( ios::showpoint ) << setprecision( 15 );

   outFile2	<< sizeInst << endl
      << overlap << endl
      << best << endl
      << solucao.size() << endl;
   for ( int i = 0; i < solucao.size(); i++ ){
      outFile2 << solucao[ i ] << endl;
   }
   for ( int j = 0; j < solucaoXYZ[ 0 ].size(); j++ ){
      outFile2 << solucaoXYZ[ 0 ][ j ] << " " << solucaoXYZ[ 1 ][ j ] << " " << solucaoXYZ[ 2 ][ j ] << endl;
   }
   outFile2.close();
}

void print_node_info_to_file(node* node, ofstream& out_file)
{
   out_file << node->depth << " " << node->lb;
   for (int item : node->pts)
   {
      out_file << " " << item;
   }
   out_file << " " << -1;
   for (int item : node->notCovered)
   {
      out_file << " " << item;
   }
   out_file << endl;
}

void print_turn_points_to_file(node* node, ofstream& out_file)
{
   out_file << node->depth << " " << node->lb;
   out_file << endl;
   // out_file << " " << -1;
   for (int item : node->notCovered)
   {
      out_file << item << " ";
   }
   out_file << endl;
   for (int item : node->pts)
   {
      out_file << item << " ";
   }
   out_file << endl;
   // for (int i = 0; i < node->solXYZ.size(); i++)
   // {
   //    for (auto item : node->solXYZ[i])
   //    {
   //       out_file << item << " ";
   //    }
   //    out_file << endl;
   // }
}