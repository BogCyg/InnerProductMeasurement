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


#include <string>
#include <iostream>

#include <chrono>
#include <ctime>
#include <cassert>
#include <functional>
#include <random>

#include "range.h"



using namespace std;




// https://en.cppreference.com/w/cpp/chrono


// Returns a string with current date & time
string	GetCurrentTime( void )
{
	using timer = std::chrono::system_clock;
    std::time_t time_point = timer::to_time_t( timer::now() );
	return std::ctime( & time_point );
}


using fun_void_int = std::function< void( int ) >;


// Launches fun( fun_param ) and measures its time
void fun_perform_timer( fun_void_int fun, int fun_param )
{
	using timer = std::chrono::high_resolution_clock;

	// ----------------------------
    auto timer_start = timer::now();
    
	// run function here ...
	fun( fun_param );

    auto timer_end = timer::now();
 	// ----------------------------

    std::chrono::duration< double > comp_period = timer_end - timer_start;

	cout << "Computation time = " << comp_period.count() << " [s]" << endl;

	cout << "Computation time = " <<  std::chrono::duration_cast< std::chrono::seconds >( comp_period ).count() << " [is]" << endl;

	cout << "Computation time = " <<  std::chrono::duration_cast< std::chrono::microseconds >( comp_period ).count() << " [us]" << endl;
}


void arithm_fun( int p )
{
	for( auto a : range( p ) )
		cout << a * a << " ";
}

void fun_perform_timer_test( void )
{
	fun_perform_timer( arithm_fun, 23 );
}


// -------------------------------------------



namespace CppBook
{
	// Nested namespace
	namespace LTimer
	{
		const int kIterations = 10;

		///////////////////////////////////////////////////////////
		// Lambda to measure execution time of a function taking 
		// any set of parameters.
		///////////////////////////////////////////////////////////
		//
		// INPUT:
		//			func - a function to be called func( func_params )
		// OUTPUT:
		//			time in milliseconds per single execution of func
		//
		// REMARKS:
		//			func will be called kIterations times and the 
		//			execution time will be averaged 
		//
		auto measureFuncAvTiming = [ iter = kIterations ]( auto && func, auto && ... func_params )               
		{
			static_assert( kIterations > 0 );	// check at compile time

			using timer = typename std::chrono::high_resolution_clock;
			
			auto time_start = timer::now();
				  
				// Run an algorithm a number of iterations 
				for( auto i : range( iter ) )
					std::forward< decltype( func ) >( func )
					( 
						std::forward< decltype( func_params ) >( func_params ) ... 
					);

			// Compute an average execution time of func in milliseconds
			return 	std::chrono::duration_cast< std::chrono::milliseconds >
						( timer::now() - time_start ).count() 
						/ kIterations;
		};

	}

}



void MathFun_RandTest( int iters, double eps, double & result )
{
	// init Mersenne twister 
	std::mt19937 mtRandomEngine( (random_device())() );	
	std::bernoulli_distribution bern( 0.7 ); // 0 or 1 with p=0.7

	result = 0.0;
	for( int i = 0; i < iters; ++ i )
		result += bern( mtRandomEngine );

	result /= static_cast< double >( iters ) + eps;
}

void TestFunctionTimeMeasurement( void )
{
	auto funTimer { CppBook::LTimer::measureFuncAvTiming };

	double result {};

	auto run_time = funTimer( MathFun_RandTest, 1234567, 1e-5, result );

	cout << "MathFun last result: " << result << endl;
	cout << "MathFun run-time: " << run_time << " [ms]" << endl;
}




// -------------------------------------------



