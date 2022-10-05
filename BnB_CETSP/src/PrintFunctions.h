#ifndef PRINTFUNCTIONS_H
#define PRINTFUNCTIONS_H

#include "structs.h"
#include "Data.h"

void printDataToFile( Data *, char*, double );
void printDataToFile ( Data *, char*, double, int, double, double, double, double, int, int, double, double, double, int, int, int, int );
void printDataToFile( Data *, char*, double, int, double, double, double, double, double, double, int, int, double, double, double, int, int, int, int );
void printDataToMatlab( Data *, int, double, double, vector< int >, vector< vector< double > > );
void print_node_info_to_file(node* node, ofstream& out_file);
void print_turn_points_to_file(node* node, ofstream& out_file);


#endif /* ifndef PRINTFUNCTIONS_H */
