//---------------------------------------------------------------------------
/*
  MeshSpace - executable to access TEdgeSpace, a class to calculate the variables for linear meshspacing

  -------------------
  created : 23-Jul-2004 Markus Hartinger    -  markus.hartinger@imperial.ac.uk
  ------------------
*/

//#include <clx.h>
#pragma hdrstop

//---------------------------------------------------------------------------

#include <iostream>
#include <cstring>
#include <cstdio>
#include <cmath>
#include "TEdgeSpace.h"
#include <cstdlib>


using namespace std;


#pragma argsused
int main(int argc, char *argv[])
{
    double L, arg1, arg2;
    TEdgeSpace ES;

	if (argc != 6) {
		cout << "Wrong syntax!\n";
		cout << argv[0] << " [argname1] [argname2] [Length] [arg1] [arg2]\n";
        cout << "Possible argument names are: 'ds', 'de', 'r', 'k', 'n'\n";
		exit(1);
	}

	L = atof(argv[3]);
	arg1 = atof(argv[4]);
	arg2 = atof(argv[5]);

    ES.calc(argv[1], argv[2], L, arg1, arg2);
//    ES.test();

	cout << "\nMeshspace results by Markus.Hartinger@imperial.ac.uk :\n";
	cout << "Length                  L : " << ES.getL() << "\n";
	cout << "Number of cells         n : " << ES.getn() << "\n";
	cout << "Total expansion ratio   r : " << ES.getr() << "\n";
	cout << "Cell-to-cell expansion  k : " << ES.getk() << "\n";
	cout << "Start cell size         ds: " << ES.getds() << "\n";
	cout << "End cell size           de: " << ES.getde() << "\n\n";

	return 0;
}
//---------------------------------------------------------------------------



