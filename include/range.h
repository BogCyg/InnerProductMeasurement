

#pragma once




#include <vector>
#include <cassert>



// In this version we use the template type deduction 

template < typename T = size_t >
class /*_range*/range 
{
		const T kFrom, kEnd, kStep;

    public:

		///////////////////////////////////////////////////////////
		// Constructor 
		///////////////////////////////////////////////////////////
		//
		// INPUT:
		//		from - Starting number of the sequence.
		//		end - Generate numbers up to, but not including this number.
		//		step -  Difference between each number in the sequence.		
		//
		// REMARKS:
		//		Parameters must be all positive or all negative
		//
		/*_range*/range( const T from, const T end, const T step = 1 ) 
			: kFrom( from ), kEnd( end ), kStep( step ) 
		{
			assert( kStep != 0 );
			//assert( ( kEnd > kFrom && kStep > 0 ) || ( kEnd < kFrom && kStep < 0 ) );
			//assert( ( kFrom >= 0 && kEnd > 0 && kStep > 0 ) || ( kFrom < 0 && kEnd < 0 && kStep < 0 ) );
		}

		// Default from==0, step==1
		/*_range*/range( const T end ) 
			: kFrom( 0 ), kEnd( end ), kStep( 1 ) 
		{
			assert( kEnd > 0 );
		}

    public:

        class /*_range*/range_iter 
		{
            T fVal;
			const T kStep;
		public:
            /*_range*/range_iter( const T v, const T step ) : fVal( v ), kStep( step ) {}
            operator T  () const			{ return fVal; }
            operator const T & ()			{ return fVal; }
            const T operator * () const		{ return fVal; }
            const /*_range*/range_iter & operator ++ ()	{ fVal += kStep; return * this; }
 

			bool operator == ( const /*_range*/range_iter & ri ) const
			{
				return ! operator != ( ri );
			}

			bool operator != ( const /*_range*/range_iter & ri ) const
			{	
				// This is a tricky part - when working with iterators
				// it checks only once for != which must be a hit to stop;
				// However, this does not work if increasing kStart by N times kSteps skips over kEnd
				//return fVal < 0 ? fVal > ri.fVal : fVal < ri.fVal;	
				assert( kStep == ri.kStep );		// This is changed to allow for mixed signs of kFrom and kEnd, such as -3 ... +3
				return kStep > 0 ? fVal < ri.fVal : fVal > ri.fVal ;	
			}      											
		};													

        const /*_range*/range_iter begin()	{ return /*_range*/range_iter( kFrom, kStep ); }
        const /*_range*/range_iter end()		{ return /*_range*/range_iter( kEnd, kStep ); }

    public:

		// Conversion to any vector< T >
		operator std::vector< T > ( void ) 
		{
			std::vector< T > retRange;
			for( T i = kFrom; i < kEnd; i += kStep )
				retRange.push_back( i );
			return retRange;	// use move semantics here
		}
};


#if 0

template < typename T >
class range : public /*_range*/range< T >
{
    public:

		range( const T from, const T end, const T step = 1 ) 
			: /*_range*/range( from, end, step ) {}

		// Default from==0, step==1
		range( const T end ) 
			:/*_range*/range( end ) {}

};


template <>
class range< int > : public /*_range*/range< int >
{
    public:

		range( const int from, const int end, const int step = 1 ) 
			: /*_range*/range( from, end, step ) {}

		// Default from==0, step==1
		range( const int end ) 
			:/*_range*/range( end ) {}

};

template <>
class range< double > : public /*_range*/range< double >
{
    public:

		range( const double from, const double end, const double step = 1 ) 
			: /*_range*/range( from, end, step ) {}

		// Default from==0, step==1
		range( const double end ) 
			:/*_range*/range( end ) {}

};

#endif


// A helper to use pure range meaning /*_range*/range< size_t >
// We need this unitil a proper type deduction is implemented in Microsoft Visual 2017
//typedef /*_range*/range<>	range;


// Let's introduce some template type deduction guides


//range( const size_t ) -> range< const size_t >;
//
//range( const unsigned int ) -> range< const unsigned int >;
//
//range( const int ) -> range< const int >;



