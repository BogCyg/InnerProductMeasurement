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


#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <random>
#include <limits>
#include <numeric>		// for inner product
#include <valarray>

#include <cassert>

#include <execution>

#include <thread>
#include <future>

#include <chrono>
#include <ctime>

#include <fstream>
#include <iterator>

#include "range.h"

#include "..\..\ttmath\ttmath.h"

#include ".\908\CPP\Src\ExactSum.h"



using namespace std;



namespace InnerProducts
{

	using DVec = vector< double >;
	using DT = DVec::value_type;
	using ST = DVec::size_type;

	using std::inner_product;
	using std::transform;
	using std::accumulate;
	using std::sort;




	/////////////////////////////////////////////////////////////////////
	// Value generators

	class FP_Test_DataSet_Generator
	{

		private:

			// xxxxxxxxxxxxxxxxxxx
			struct _908_sandbox
			{
	
				//int pflag;  //program flag (Data type 1 to 4)
				int pflag { 1 };  //program flag (Data type 1 to 4)
				// random flop generator with exponent differece expo

				int rflag=0;

				double ret_saved;

				double Rand(int expo)
				{
					double ret;
					unsigned int low=0,high=0;
					int i;
					for(i=0;i<32;i++)
					{
						if(rand() & 0x1)
							low=(low>>1) | 0x80000000;
						else
							low=low>>1;
					}
					for(i=0;i<32;i++)
					{
						if(rand() & 0x1)
							high=(high>>1) | 0x80000000;
						else
							high=high>>1;
					}
					*(int*)(&ret)=low;
					*(((int*)&ret)+1)=high;
					if(expo!=0)
						((str_double*)&ret)->exponent=((rand()%(expo))-expo/2)+0x3ff;
					else 
						((str_double*)&ret)->exponent=0x3ff;
					if(pflag==2 || pflag==3)
						return ret;  //pflag=2, random data
					if(pflag==1)
						rflag=0;     //pflag=1, well-conditioned data
					if(rflag==0)  // exact zero
					{
						ret_saved=ret;
						((str_double*)(&ret))->sign=1;
						rflag=1;
					}
					else
					{
						ret=ret_saved;
						((str_double*)(&ret))->sign=0;
						rflag=0;
					}
					return ret;
				}


				void Generate( DVec & inVec, int deltaExp, int _pFlag )
				{
					assert( _pFlag >= 1 && _pFlag <= 4 );		// inherited from 908
					pflag = _pFlag;

					int MAXNUM = inVec.size();


					srand((unsigned)time(NULL));
					int i;

					double *original_list = & inVec[ 0 ];

					double st=0;
	


					for(i=0;i<MAXNUM;i++)
					{
						original_list[i]=(Rand(deltaExp));
						st+=original_list[i];
					}


					// Anderson's ill-conditioned data
					if(pflag==3) 
						for(i=0;i<MAXNUM;i++)
							original_list[i]-=st/MAXNUM;

					//randomly change the order
					for(i=0;i<MAXNUM*2;i++)
					{
						int x=rand()%MAXNUM;
						int y=rand()%MAXNUM;
						while(x==y)
							y=rand()%MAXNUM;
						double temp=original_list[x];
						original_list[x]=original_list[y];
						original_list[y]=temp;
					}




				}
			};

			// xxxxxxxxxxxxxxxxxxx

		public:

			// (1, well-conditioned; 2, random; 3, Anderson's; and
			// 4, exact sum equals zero)
			void Fill_Numerical_Data_No( int flag, DVec & inVec, ST num_of_data, int deltaExp = 10 )
			{
				_908_sandbox	_908_obj;
				inVec.resize( num_of_data );
				_908_obj.Generate( inVec, deltaExp, flag );
			}


			void Fill_Numerical_Data_MersenneUniform( DVec & inVec, ST num_of_data, DT kDataMag )
			{
				inVec.resize( num_of_data );
				assert( inVec.size() == num_of_data );
				mt19937		rand_gen{ random_device{}() };	// Random Mersenne twister
				uniform_real_distribution< double > dist( - kDataMag, + kDataMag );
				std::generate( inVec.begin(), inVec.end(), [&](){ return dist( rand_gen ); } );
			}

			// 
			void Duplicate( DVec & inVec, DT multFactor = 1.0 )
			{
				ST kElems { inVec.size() };
				inVec.resize( 2 * kElems );
				std::generate( inVec.begin() + kElems, inVec.end(), [ n = 0, & inVec, multFactor ] () mutable { return multFactor * inVec[ n ++ ]; } );
			}

			void DuplicateWithNegated( DVec & inVec )
			{
				Duplicate( inVec, -1.0 );
			}

	};


}


/////////////////////////////////////////////////////////////////////
// Inner product




namespace InnerProducts
{

	auto InnerProduct_StdAlg( const DVec & v, const DVec & w )
	{
		// The last argument is an initial value
		return std::inner_product( v.begin(), v.end(), w.begin(), DT() );
	}


	// The transform-reduce parallel version
	auto InnerProduct_TR_Alg( const DVec & v, const DVec & w )
	{
		return std::transform_reduce(	std::execution::par,
										v.begin(), v.end(), w.begin(), DT(),
										[] ( const auto a, const auto b ) { return a + b; },
										[] ( const auto a, const auto b ) { return a * b; }
			);
	}


	auto InnerProduct_SortAlg( const DVec & v, const DVec & w )
	{
		DVec z( std::min( v.size(), w.size() ) );		// Stores element-wise products

		// Elementwise multiplication: c = a .* b

		transform(	std::execution::par, v.begin(), v.end(), w.begin(), 
					z.begin(), 
					[] ( const auto & v_el, const auto & w_el) { return v_el * w_el; } );

		// Sort in descending order.
		std::sort( z.begin(), z.end(), [] ( const DT & p, const DT & q ) { return fabs( p ) < fabs( q ); } );		// Is it magic?

		// The last argument is an initial value
		return accumulate( z.begin(), z.end(), DT() );
	}



	auto Sort_And_Accumulate( DVec & v )
	{
		// Having sort we need to apply a SERIAL accumulate to have GOOD results.
		// This happens because std::reduce will brake the order. Simple.
		std::sort( std::execution::par, v.begin(), v.end(), [] ( const DT & p, const DT & q ) { return fabs( p ) < fabs( q ); } );		// Is it magic?																																		
		return accumulate( v.begin(), v.end(), DT() );		// This is IMPORTANT - we can sort in PARALLEL, but then we must accumulate in SERIAL (not to spoil the order)
	}


	// ACTUALLY IF WE NEED TO SORT, THEN WE DO NOT NEED THE KAHAN ALGORITHM
	// since the sort-and-accumulate is the best
	// In the Kahan algorithm each addition is corrected by a correction
	// factor. In this algorithm the non associativity of FP is used, i.e.:
	// ( a + b ) + c != a + ( b + c )
	// v will be changed
	auto Kahan_Sum( DVec & v )
	{
		DT theSum {};

		// volatile prevents a compiler from applying any optimization
		// on the object since it can be changed by someone else, etc.,
		// in a way that cannot be foreseen by the compiler.

		volatile DT c {};		// a "correction" coefficient

		for( ST i = 0; i < v.size(); ++ i )
		{
			DT y = v[ i ] - c;			// From the summand y subtract the correction factor

			DT t = theSum + y;			// Add corrected summand to the running sum, i.e. theSum
										// But theSum is big, y is small, so its lower bits will be lost

			c = ( t - theSum ) - y;		// Low order bits of y are lost in the summation. High order
										// bits of y are computed in ( t - theSum ). Then, when y
										// is subtracted from this, the low order bits of y are recovered (negative).
										// Algebraically, c should always be 0 (beware of compiler optimization).
			theSum = t;
		}

		return theSum;
	}


	
	// ACTUALLY IF WE NEED TO SORT, THEN WE DO NOT NEED THE KAHAN ALGORITHM
	// since the sort-and-accumulate is the best
	// In the Kahan algorithm each addition is corrected by a correction
	// factor. In this algorithm the non associativity of FP is used, i.e.:
	// ( a + b ) + c != a + ( b + c )
	// v will be changed
	auto Kahan_Sort_And_Sum( DVec & v )
	{
		std::sort( std::execution::par, v.begin(), v.end(), [] ( const DT & p, const DT & q ) { return fabs( p ) < fabs( q ); } );		// Is it magic?
		return Kahan_Sum( v );
	}




	// In the Kahan algorithm each addition is corrected by a correction
	// factor. In this algorithm the non associativity of FP is used, i.e.:
	// ( a + b ) + c != a + ( b + c )
	auto InnerProduct_KahanAlg( const DVec & v, const DVec & w )
	{
		DT theSum {};

		// volatile prevents a compiler from applying any optimization
		// on the object since it can be changed by someone else, etc.,
		// in a way that cannot be foreseen by the compiler.

		volatile DT c {};		// a "correction" coefficient

		const ST kElems = std::min( v.size(), w.size() );

		for( ST i = 0; i < kElems; ++ i )
		{
			DT y = v[ i ] * w[ i ] - c;	// From the summand y subtract the correction factor

			DT t = theSum + y;			// Add corrected summand to the running sum, i.e. theSum
										// But theSum is bit, y is small, so its lower bits will be lost

			c = ( t - theSum ) - y;		// Low order bits of y are lost in the summation. High order
										// bits of y are computed in ( t - theSum ). Then, when y
										// is subtracted from this, the low order bits of y are recovered (negative).
										// Algebraically, c should always be 0 (beware of compiler optimization).
			theSum = t;
		}

		return theSum;
	}


	// Other version of the Kahan algorithms
	auto InnerProduct_KahanAlg( const double * v, const double * w, const size_t kElems )
	{
		DT theSum {};

		// volatile prevents a compiler from applying any optimization
		// on the object since it can be changed by someone else, etc.,
		// in a way that cannot be foreseen by the compiler.

		volatile DT c {};		// a "correction" coefficient


		for( ST i = 0; i < kElems; ++ i )
		{
			DT y = v[ i ] * w[ i ] - c;	// From the summand y subtract the correction factor

			DT t = theSum + y;			// Add corrected summand to the running sum, i.e. theSum
										// But theSum is bit, y is small, so its lower bits will be lost

			c = ( t - theSum ) - y;		// Low order bits of y are lost in the summation. High order
										// bits of y are computed in ( t - theSum ). Then, when y
										// is subtracted from this, the low order bits of y are recovered (negative).
										// Algebraically, c should always be 0 (beware of compiler optimization).
			theSum = t;
		}

		return theSum;
	}





	// Test the two
	auto InnerProduct_Sort_KahanAlg( const DVec & v, const DVec & w )
	{
		DVec z( std::min( v.size(), w.size() ) );		// Stores element-wise products

		// Elementwise multiplication: c = a .* b
		transform(	std::execution::par, v.begin(), v.end(), w.begin(), 
					z.begin(), 
					[] ( const auto & v_el, const auto & w_el) { return v_el * w_el; } );


		std::sort( std::execution::par, z.begin(), z.end(), [] ( const DT & p, const DT & q ) { return fabs( p ) < fabs( q ); } );		// Is it magic?

		// ------------------------



		DT theSum {};

		// volatile prevents a compiler from applying any optimization
		// on the object since it can be changed by someone else, etc.,
		// in a way that cannot be foreseen by the compiler.

		volatile DT c {};		// a "correction" factor

		const ST kElems = z.size();

		for( ST i = 0; i < kElems; ++ i )
		{
			DT y = z[ i ] - c;	// From the summand y subtract the correction factor

			DT t = theSum + y;			// Add corrected summand to the running sum, i.e. theSum
										// But theSum is bit, y is small, so its lower bits will be lost

			c = ( t - theSum ) - y;		// Low order bits of y are lost in the summation. High order
										// bits of y are computed in ( t - theSum ). Then, when y
										// is subtracted from this, the low order bits of y are recovered (negative).
										// Algebraically, c should always be 0 (beware of compiler optimization).
			theSum = t;
		}

		return theSum;
	}
	


	
	// 2nd version
	auto InnerProduct_Sort_KahanAlg(  const double * v, const double * w, const size_t kElems  )
	{
		DVec z;		// Stores element-wise products

		// Elementwise multiplication: c = a .* b
		transform(	v, v + kElems, w, 
					back_inserter( z ), 
					[] ( const auto & v_el, const auto & w_el) { return v_el * w_el; } );

		std::sort( z.begin(), z.end(), [] ( const DT & p, const DT & q ) { return fabs( p ) < fabs( q ); } );		// Is it magic?

		// ------------------------



		DT theSum {};

		// volatile prevents a compiler from applying any optimization
		// on the object since it can be changed by someone else, etc.,
		// in a way that cannot be foreseen by the compiler.

		volatile DT c {};		// a "correction" factor


		for( ST i = 0; i < z.size(); ++ i )
		{
			DT y = z[ i ] - c;	// From the summand y subtract the correction factor

			DT t = theSum + y;			// Add corrected summand to the running sum, i.e. theSum
										// But theSum is bit, y is small, so its lower bits will be lost

			c = ( t - theSum ) - y;		// Low order bits of y are lost in the summation. High order
										// bits of y are computed in ( t - theSum ). Then, when y
										// is subtracted from this, the low order bits of y are recovered (negative).
										// Algebraically, c should always be 0 (beware of compiler optimization).
			theSum = t;
		}

		return theSum;
	}




	//////////////




	namespace PrecLongComp
	{
		//                              exp               mantissa
		using BNum = ttmath::Big< TTMATH_BITS( /*12*//*24*//*128*/384 ), TTMATH_BITS( /*54*//*128*//*256*/512 ) >;
		using ExtraBNum = ttmath::Big< TTMATH_BITS( /*242*//*128*//*256*/512 ), TTMATH_BITS( /*128*//*512*/1024 ) >;

		using InnerProducts::DVec;

		auto InnerProduct_BNum( const DVec & v, const DVec & w )
		{
			// The last argument is an initial value
			//return std::inner_product( v.begin(), v.end(), w.begin(), DT() );

			const DVec::size_type	kSize( std::min( v.size(), w.size() ) );

			ExtraBNum	theSum  = 0;	// The initialization is important here - do not use the initializer list

			BNum	tmp = 0;

			double d {};

			for( auto r : range( kSize ) )
			{
				tmp = v[ r ];

				tmp *= w[ r ];

				theSum = theSum + tmp;
			}

			return theSum;
		}



	}

	//////////////

	// EX == ExactSum
	namespace ES
	{

		auto DummyTest_908( void ) -> double
		{
			DVec z;		// Stores element-wise products

			z.push_back( 0.0 );		// we need that since the guy is starting from index 1, well...

			z.push_back( 1 );
			z.push_back( 2 );
			z.push_back( 3 );
			z.push_back( 4 );
			z.push_back( 5 );
			z.push_back( 6 );
			z.push_back( 7 );
			z.push_back( 8 );
	
			ExactSum mysum;

			mysum.Reset();

			return mysum.OnlineExactSum( & z[ 1 ], z.size() - 1 );
		}

	
		auto InnerProduct_908( const DVec & v, const DVec & w )
		{
			assert( false );			// something is wrong here


			DVec z;		// Stores element-wise products

			z.push_back( 0.0 );		// we need that since the guy is starting from index 1, well...

			// Elementwise multiplication: c = a .* b
			// Add parallelization to transform
			transform(	v.begin(), v.end(), w.begin(), 
						back_inserter( z ), 
						[] ( const auto & v_el, const auto & w_el) { return v_el * w_el; } );
	
			ExactSum mysum;

			mysum.Reset();

			return mysum.OnlineExactSum( & z[ 1 ], z.size() - 1 );
		}


		auto Sum_908( const DVec & v )
		{
			ExactSum mysum;

			mysum.Reset();

			for( auto x : v )
				mysum.AddNumber( x );
			
			return mysum.GetSum();
		}
		auto InnerProduct_908_b( const DVec & v, const DVec & w )
		{
			ExactSum mysum;

			mysum.Reset();

			const auto kSize { std::min( v.size(), w.size() ) };
			for( auto i : range( kSize ) )
				mysum.AddNumber( v[ i ] * w[ i ] );
			
			return mysum.GetSum();
		}


		auto InnerProduct_908_c( const DVec & v, const DVec & w )
		{
			DVec z( v.size() );		// Stores element-wise products

			z.push_back( 0.0 );		// we need that since the guy is starting from index 1, well...

			// Elementwise multiplication: c = a .* b
			// Add parallelization to transform
			transform(	v.begin(), v.end(), w.begin(), 
						z.begin() + 1, 
						[] ( const auto & v_el, const auto & w_el) { return v_el * w_el; } );
	
			ExactSum mysum;

			mysum.Reset();

			mysum.AddArray( & z[ 0 ], v.size() );

			return mysum.GetSum();
		}

		auto InnerProduct_908_par( const double * v, const double * w, const size_t kElems )
		{
			ExactSum mysum;
			mysum.Reset();

			for( auto i : range( kElems ) )
				mysum.AddNumber( v[ i ] * w[ i ] );
			
			return mysum.GetSum();
		}


	}

	//////////////


	// THE BEST PERFORMANCE
	// This is a simple data paralellization of the Kahan algorithm.
	// The input vectors are divided into the chunks which are 
	// then processed in parallel but by the serial Kahan algorithm.
	// The partial sums are then summed up with yet run of the
	// Kahan algorithm.
	auto InnerProduct_KahanAlg_Par( const DVec & v, const DVec & w, const ST kChunkSize = 10000 )
	{
		const auto kMinSize { std::min( v.size(), w.size() ) };

		const auto k_num_of_chunks { ( kMinSize / kChunkSize ) };
		const auto k_remainder { kMinSize % kChunkSize };

		const double * v_data_begin = & v[ 0 ];
		const double * w_data_begin = & w[ 0 ];

		vector< double >	par_sum( k_num_of_chunks + ( k_remainder > 0 ? 1 : 0 ), 0.0 );	

		// The thing is that we wish Kahan because it is much faster than the sort-accum
		auto fun_inter = [] ( const double * a, const double * b, int s ) { return InnerProduct_KahanAlg( a, b, s ); };

		vector< future< double > >		my_thread_poool;

		// Process all equal size chunks of data
		std::decay< decltype( k_num_of_chunks ) >::type i {}; 
		for( i = 0; i < k_num_of_chunks; ++ i )
			my_thread_poool.push_back( async( std::launch::async, fun_inter, v_data_begin + i * kChunkSize, w_data_begin + i * kChunkSize, kChunkSize ) );

		// Process the ramainder, if present
		if( k_remainder > 0 )
			my_thread_poool.push_back( async( std::launch::async, fun_inter, v_data_begin + i * kChunkSize, w_data_begin + i * kChunkSize, k_remainder ) );

		assert( par_sum.size() == my_thread_poool.size() );
		for( i = 0; i < my_thread_poool.size(); ++ i )
			par_sum[ i ] = my_thread_poool[ i ].get();			// get() bocks until the async is done

		return Kahan_Sort_And_Sum( par_sum );			
	}


	auto InnerProduct_SortKahanAlg_Par( const DVec & v, const DVec & w, const size_t kChunkSize = 10000 )
	{
		const auto kMinSize { std::min( v.size(), w.size() ) };

		const auto k_num_of_chunks { ( kMinSize / kChunkSize ) };
		const auto k_remainder { kMinSize % kChunkSize };


		const double * v_data_begin = & v[ 0 ];
		const double * w_data_begin = & w[ 0 ];

		vector< double >	par_sum( k_num_of_chunks + ( k_remainder > 0 ? 1 : 0 ), 0.0 );	

		// The thing is that we wish Kahan because it is much faster than the sort-accum
		auto fun_inter = [] ( const double * a, const double * b, int s ) { return InnerProduct_Sort_KahanAlg( a, b, s ); };

		vector< future< double > >		my_thread_poool;

		// Process all equal size chunks of data
		std::decay< decltype( k_num_of_chunks ) >::type i {};
		for( i = 0; i < k_num_of_chunks; ++ i )
			my_thread_poool.push_back( async( std::launch::async, fun_inter, v_data_begin + i * kChunkSize, w_data_begin + i * kChunkSize, kChunkSize ) );

		// Process the ramainder, if present
		if( k_remainder > 0 )
			my_thread_poool.push_back( async( std::launch::async, fun_inter, v_data_begin + i * kChunkSize, w_data_begin + i * kChunkSize, k_remainder ) );

		assert( par_sum.size() == my_thread_poool.size() );
		for( size_t i = 0; i < my_thread_poool.size(); ++ i )
			par_sum[ i ] = my_thread_poool[ i ].get();


		return Kahan_Sort_And_Sum( par_sum );		
	}




	auto InnerProduct_908_Par( const DVec & v, const DVec & w, const size_t kChunkSize = 10000 )
	{
		const auto kMinSize { std::min( v.size(), w.size() ) };

		const size_t k_num_of_chunks { ( kMinSize / kChunkSize ) };
		const size_t k_remainder { kMinSize % kChunkSize };


		const double * v_data_begin = & v[ 0 ];
		const double * w_data_begin = & w[ 0 ];

		vector< double >	par_sum( k_num_of_chunks + ( k_remainder > 0 ? 1 : 0 ), 0.0 );

		// The lambda for serial summation
		auto fun_inter = [] ( const double * a, const double * b, int s ) { return ES::InnerProduct_908_par( a, b, s ); };

		vector< future< double > >		my_thread_poool;

		// Process all equal size chunks of data
		std::decay< decltype( k_num_of_chunks ) >::type i {};
		for( i = 0; i < k_num_of_chunks; ++ i )
			my_thread_poool.push_back( async( std::launch::async, fun_inter, v_data_begin + i * kChunkSize, w_data_begin + i * kChunkSize, kChunkSize ) );

		// Process the ramainder, if present
		if( k_remainder > 0 )
			my_thread_poool.push_back( async( std::launch::async, fun_inter, v_data_begin + k_num_of_chunks * kChunkSize, w_data_begin + k_num_of_chunks * kChunkSize, k_remainder ) );

		assert( par_sum.size() == my_thread_poool.size() );
		for( auto i : range(  my_thread_poool.size() ) )
			par_sum[ i ] = my_thread_poool[ i ].get();

		return ES::Sum_908( par_sum );		
	}




	void InnerProduct_Test_3( DVec & v, DVec & w )
	{
		assert( v.size() == w.size() );

		// Both dimensions of v and w must be the same
		// and must be an integer multiplication of the kChunkSize
		const int kChunkSize { /*10000*/25000/*24000*/ };


		// The inner product should be close to 0.0, 
		// so let us check the two algorithms.

		using timer = typename std::chrono::high_resolution_clock;
		auto get_duration = [] ( auto prev_time_stamp ) 
							{ return std::chrono::duration_cast< std::chrono::milliseconds >( 
													timer::now() - prev_time_stamp ).count(); };


		DVec	result_errors;
		DVec	result_timing;


		auto comp_error { 0.0 };

		auto ts = timer::now();
		comp_error = fabs( InnerProduct_StdAlg( v, w ) );
		auto tdur = get_duration( ts );
		cout << "Stand alg error = \t"	<< std::setprecision( 8 ) << comp_error << "\t\tT [ms] = " << tdur << endl;
		result_errors.push_back( comp_error );
		result_timing.push_back( tdur );


		ts = timer::now();
		comp_error = fabs( InnerProduct_TR_Alg( v, w ) );
		tdur = get_duration( ts );
		cout << "Parallel Transform-Reduce alg error = \t"	<< std::setprecision( 8 ) << comp_error << "\t\tT [ms] = " << tdur << endl;
		result_errors.push_back( comp_error );
		result_timing.push_back( tdur );
		
		ts = timer::now();
		comp_error = fabs( InnerProduct_SortAlg( v, w ) );
		tdur = get_duration( ts );
		cout << "Sort alg error = \t"	<< std::setprecision( 8 ) << comp_error << "\t\tT [ms] = " << tdur << endl;
		result_errors.push_back( comp_error );
		result_timing.push_back( tdur );
		
		ts = timer::now();
		comp_error = fabs( InnerProduct_KahanAlg( v, w ) );
		tdur = get_duration( ts );
		cout << "Kahan alg error = \t"	<< std::setprecision( 8 ) << comp_error << "\t\tT [ms] = " << tdur << endl;
		result_errors.push_back( comp_error );
		result_timing.push_back( tdur );
		
		ts = timer::now();
		comp_error = fabs( InnerProduct_Sort_KahanAlg( v, w ) );
		tdur = get_duration( ts );
		cout << "Serial Sort-Kahan alg error = \t"	<< std::setprecision( 8 ) << comp_error << "\t\tT [ms] = " << tdur << endl;
		result_errors.push_back( comp_error );
		result_timing.push_back( tdur );
		
	
		ts = timer::now();
		comp_error = fabs( InnerProduct_KahanAlg_Par( v, w, kChunkSize ) );
		tdur = get_duration( ts );
		cout << "Parallel Kahan alg error = \t"	<< std::setprecision( 8 ) << comp_error << "\t\tT [ms] = " << tdur << endl;
		result_errors.push_back( comp_error );
		result_timing.push_back( tdur );


		ts = timer::now();
		comp_error = fabs( InnerProduct_SortKahanAlg_Par( v, w, kChunkSize ) );
		tdur = get_duration( ts );
		cout << "Parallel Sort-Kahan alg error = \t"	<< std::setprecision( 8 ) << comp_error << "\t\tT [ms] = " << tdur << endl;
		result_errors.push_back( comp_error );
		result_timing.push_back( tdur );

		ts = timer::now();
		comp_error = fabs( InnerProduct_908_Par( v, w, kChunkSize ) );
		tdur = get_duration( ts );
		cout << "Parallel 908 alg error = \t"	<< std::setprecision( 8 ) << comp_error << "\t\tT [ms] = " << tdur << endl;
		result_errors.push_back( comp_error );
		result_timing.push_back( tdur );


		// ---------
		// 908
		ts = timer::now();
		comp_error = fabs( ES::InnerProduct_908_b( v, w ) );
		tdur = get_duration( ts );
		cout << "Serial 908 alg error = \t"	<< std::setprecision( 8 ) << comp_error << "\t\tT [ms] = " << tdur << endl;
		result_errors.push_back( comp_error );
		result_timing.push_back( tdur );

		// ---------
		// long precision lib
		ts = timer::now();
		comp_error = fabs( PrecLongComp::InnerProduct_BNum( v, w ).ToDouble() );
		tdur = get_duration( ts );
		cout << "ttmath alg error = \t"	<< std::setprecision( 8 ) << comp_error << "\t\tT [ms] = " << tdur << endl << "- - -" << endl << endl;
		result_errors.push_back( comp_error );
		result_timing.push_back( tdur );





		// -----------------------------------------------
		// Save results
		ofstream res_file( "inner_results.txt", ios::app );
		copy( result_errors.begin(), result_errors.end(), ostream_iterator< double >( res_file, "\t" ) );	res_file << endl;
		copy( result_timing.begin(), result_timing.end(), ostream_iterator< double >( res_file, "\t" ) );	res_file << endl << endl;
		// -----------------------------------------------
	}



	// Run the InnerProduct_Test a number of times
	// with different range of the generated data samples.
	void InnerProduct_Test_GeneralExperiment( void )
	{

		const int kElems = /*100*//*40000000*//*1000000*/20000000/*10000*/;

		vector< int >	deltaExpVec { 10, 30, 50, 100, 300, 500/*, 1000, 2000*/ };

		DVec	v, w;
	
		FP_Test_DataSet_Generator		data_generator;


		enum class FP_TestData_Type { kWellConditioned, kRandom, kAnderson, kExactSumIsZero, kMersenneRand_InnerZero };


		for( auto dtype : range( (int) FP_TestData_Type::kMersenneRand_InnerZero + 1 ) )
		{
			FP_TestData_Type	data_type { dtype };

			for( auto dExp : deltaExpVec )
			{
				switch( data_type )
				{
					// -------------------------------------
					case FP_TestData_Type::kWellConditioned:
						cout << "kWellConditioned" << endl;
						data_generator.Fill_Numerical_Data_No( 1, v, kElems, dExp );
						w.resize( kElems, 1.0 );	
						break;

					// ----------------------------
					case FP_TestData_Type::kRandom:
						cout << "kRandom" << endl;
						data_generator.Fill_Numerical_Data_No( 2, v, kElems, dExp );
						w.resize( kElems, 1.0 );	
						break;
					
					// -----------------------------
					case FP_TestData_Type::kAnderson:
						cout << "kAnderson" << endl;
						data_generator.Fill_Numerical_Data_No( 3, v, kElems, dExp );
						w.resize( kElems, 1.0 );
						break;
						
					// ------------------------------------
					case FP_TestData_Type::kExactSumIsZero:
						cout << "kExactSumIsZero" << endl;
						data_generator.Fill_Numerical_Data_No( 4, v, kElems, dExp );
						w.resize( kElems, 1.0 );
						break;
							
					// --------------------------------------------
					case FP_TestData_Type::kMersenneRand_InnerZero:
						cout << "kMersenneRand_InnerZero" << endl;
						data_generator.Fill_Numerical_Data_MersenneUniform( v, kElems / 2, pow( 2.0, dExp ) );
						data_generator.Duplicate( v, + 1.0 );
						data_generator.Fill_Numerical_Data_MersenneUniform( w, kElems / 2, pow( 2.0, dExp ) );
						data_generator.Duplicate( w, - 1.0 );
						break;

					default: assert( false );
						break;						
			
				}

				cout << "ExpDelta = " << dExp << "\tVecElems = " << v.size() << endl;
				InnerProduct_Test_3( v, w );
			}

		}


	}




}	// end of namespace


