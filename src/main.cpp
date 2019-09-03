///////////////////////////////////////////////////////
// Written by Boguslaw Cyganek, 2019
///////////////////////////////////////////////////////
//
// When using this code please cite the following paper:
//
// "How orthogonal are we? A note on fast and accurate 
// inner product computation in the floating-point arithmetic"
// by Boguslaw Cyganek and Kazimierz Wiatr 
// First International Conference on SOCIETAL AUTOMATION
// September 4-6, 2019, Krakow, Poland
//
///////////////////////////////////////////////////////


#include <iostream>
#include <string>

namespace InnerProducts
{
	void InnerProduct_Test( double );
	void InnerProduct_Test_GeneralExperiment( void );
}


using std::cout, std::endl, std::cin;
using std::string;

string	GetCurrentTime( void );


//////////////

int main( void )
{

	cout << "===========================" << endl;
	cout << "Inner Product Test - let's begin!" << endl;
	cout << GetCurrentTime();
	cout << "===========================" << endl << endl;


	InnerProducts::InnerProduct_Test_GeneralExperiment();

}









