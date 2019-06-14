///////////////////////////////////////////////////////
// BC++ book
///////////////////////////////////////////////////////
// by Boguslaw Cyganek, Wiley, 2019
///////////////////////////////////////////////////////


#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <filesystem>


#include "..\..\ttmath\ttmath.h"




using namespace std;

namespace fs = std::experimental::filesystem;

string	GetCurrentTime( void );

namespace InnerProducts
{
	void InnerProduct_Test( double );
	void InnerProduct_Test_GeneralExperiment( void );
}



//////////////
// ttmath test
// Big<exponent, mantissa>
typedef ttmath::Big<TTMATH_BITS(64), TTMATH_BITS(128)> MyBig;

// See
// D:\Lunert\Precision Arithemtic Libs\TTMATH\ttmath-0.9.4.prerelease-src-2017.03.12\samples

void TT_Test( void )
{
	MyBig	mb { 0 };

	double qpi { 3.14151 };

	mb.FromDouble( qpi );
	//mb = qpi;

	MyBig	sum { - 3.14151 };

	sum += mb;


	std::cout << "mb = " << mb.ToString() << ", qpi = " << qpi << endl;
	std::cout << "sum = " << sum << endl;


}

//////////////

void ES_Test( void )
{
	int ES_Test_main(int argn, char* argc[]);

	char * arr[] = {
						"",				// empty to simulate the environment
						"1000",			// total elems to sum
						"10",			// exponent difference
						"3"				// what-to-do code
	};


	ES_Test_main( sizeof(arr)/sizeof(arr[0]), arr );

}

namespace InnerProducts
{
	namespace ES
	{

		auto DummyTest_908( void ) -> double;
	}
}

//////////////

int paranoia_main( void );

//////////////

int main( void )
{

	cout << "===========================" << endl;
	cout << "Inner Product Test - let's begin!" << endl;
	cout << GetCurrentTime();
	cout << "===========================" << endl << endl;


	//paranoia_main();


	//TT_Test();
	//ES_Test();
	//cout << "ES::DummyTest_908() = " << InnerProducts::ES::DummyTest_908() << endl;

	//InnerProducts::InnerProduct_Test();
	InnerProducts::InnerProduct_Test_GeneralExperiment();

	char c {};
	cin >> c;

	return c;
}









