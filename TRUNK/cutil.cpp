#ifndef IONCH_CUTIL_CPP
#define IONCH_CUTIL_CPP 1
/*----------------------------------------------------------------------
  This source file is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
----------------------------------------------------------------------*/


/* --------------------------------------------------------------
 Routines written in C/C++ that can be called from Fortran.

 These methods are not specific to any project.
   - methods to access machine information via the C library
   - bit manipulation methods for testing fuzzy equality
   - bridge to random number generator in external libraries
   - access to C preprocessor definitions
     -- version information for GSL and ATLAS maths libraries
 */

#ifdef _OPENMP
#include "omp.h"
#endif
#ifndef HAVE_CSTDINT
extern "C"
{
#include <stdint.h>
}
#else
#include <cstdint>
#endif
#ifndef HAVE_LLRINT
extern "C"
{
    long long int llrint ();
    long long int llrintl ();
    long long int llrintf ();
    long long int llround ();
    long long int llroundf ();
    long long int llroundl ();
}
# endif
extern "C"
{
#include <sys/stat.h>
}
#include <cerrno>
#include <cstdio>
#include <cstring>

#include <limits>
#include <algorithm>


// Signal handling value and method

#include <csignal>
namespace
{
static int *onterm;
void on_term_signal(int)
{
  *onterm = 1;
}
}

/* --------------------------------------------------------------
  Determine if two floating point numbers are very close to being equal.
*/
int fuzzy_equal(double a, double b);

namespace
{
/* --------------------------------------------------------------
   Generate a random number between 0 and 1
*/
inline void seedrand (int seed);
inline double genrand ();

  void random_shuffle_int(int* array, int sz);

/* --------------------------------------------------------------
  Method that writes details about the computer that
  are difficult (or non-standard) to obtain in fortran
  95:

  hostname and arch type
  os name, type and version
  cache type
*/
  void write_host_info(char *, int*, int*);

/* --------------------------------------------------------------
  Method that generates a UUID string for use in FORTRAN

  This is based on the OSF DCE compatible Universally Unique 
  Identifier library
*/
void ionch_uuid_gen(char *const the_uuid);

/* ------------------------------------------------------------
  Method to convert string of floats into an array of floats

  \param strval : input string containing one or more floats
  \param a_array: output array of doubles
  \param a_size : the size of a_array
  \param a_cnt  : the output number of doubles read
  \param a_istt : error indicator, 0 if no errors
*/
void read_floats(char const* strval, double *a_arry, int a_size, int * a_cnt, int * a_istt);

/* ------------------------------------------------------------
  Method to pretty print a float
  \param a_dbl  : FP number to write
  \param a_result: Result to write.
*/
void pretty_print_number(double *a_dbl, char *a_result, int *sz);

#ifndef USE_MKL
/* --------------------------------------------------------------
  Get version information for the maths library
  for GSL or ATLAS libraries
*/
void mathlib_version(char *const buf);
#endif
}

#ifdef HAVE_NO_IMPLICIT_ISNAN
#include <cmath>
#endif

extern "C"
{
/* -------------------------------------------------------------- 
  C-TO-FORTRAN definitions

  The 'cfortran.h' library is used to allow Fortran code to
use C routines.  The library handles matching argument
and return types between the languages.

  NOTE: These definitions are in an 'extern C' scope to
    avoid C++ name mangling of the Fortran visible names for the
    methods.  The declaration and definitions of the methods are
    contained in an unnamed namespace to make them inaccessible
    except via the Fortran visible names.

  NOTE: 'cfortran.h' requires specification of the compiler to adjust
    for compiler differences.  The Intel and GNU gfortran compilers
    are not directly supported but appear to be compatible with 
    sunFortran in 'cfortran.h' file
 */
#define sunFortran
#include "cfortran.h"
FCALLSCFUN2(LOGICAL,fuzzy_equal,DFEQ,dfeq,DOUBLE,DOUBLE)

FCALLSCSUB1(seedrand,SRANFF,sranff,INT)

FCALLSCSUB1(ionch_uuid_gen,UUIDGN,uuidgn,PPSTRING)

FCALLSCSUB3(write_host_info,MACHINE_INFO,machine_info,PPSTRING,PINT,PINT)

FCALLSCSUB2(random_shuffle_int,SHUFFLE,shuffle,PINT,INT)

FCALLSCSUB5(read_floats,REDFLT,redflt,STRING,PDOUBLE,INT,PINT,PINT)

FCALLSCFUN0(DOUBLE,genrand,RANFF,ranff)

FCALLSCSUB3(pretty_print_number,FMTFL2,fmtfl2,PDOUBLE,PSTRING,PINT)

#ifndef USE_MKL
FCALLSCSUB1(mathlib_version,MATH_LIB_VERSION,math_lib_version,PPSTRING)
#endif

#ifdef HAVE_NO_IMPLICIT_ISNAN
FCALLSCFUN1(LOGICAL,std::isnan,ISNAN,isnan,DOUBLE)
#endif

}

namespace
{

void pretty_print_number(double *a_dbl, char *a_result, int *sz)
{
  char fmt[16];
  fmt[0]='%';
  if (*sz < 9)
    snprintf(&fmt[1],15,"-%d.%dG",*sz,2);
  else
    snprintf(&fmt[1],15,"-%d.%dG",*sz,*sz-7);
  snprintf(a_result,*sz,fmt,*a_dbl);
}
}
/* --------------------------------------------------------------
  Get includes for random numbers depending on whether you
  are using SFMT or C++ TR1
*/
#ifndef USE_SFMT
# include <tr1/random>
# include <memory>
#else
# ifndef HAVE_SSE2
#  ifdef __SSE2__
#   define HAVE_SSE2 1
#  endif
# endif
extern "C"
{
#if ( USE_SFMT == 1 )
# define MEXP 19937
# include "SFMT.c"
#else
# define DSFMT_MEXP 19937
# include "dSFMT.c"
#endif
}
#endif

namespace
{
// Traits class to get integer of same size as floating point type
template < class Float_type >
struct ftoi
{
    typedef int (*int_type)();
};
template<>
struct ftoi<float>
{
    typedef int32_t int_type;
};
template<>
struct ftoi<double>
{
    typedef int64_t int_type;
};

/* --------------------------------------------------------------
 * FUZZY EQUALS.

 As floating point numbers are not exact, two values may not
 compare equal when they should.  This function compares two
 number using '=='.  When this is not true it then tests whether
 the last bits of the floating point numbers are +/- 128
 */
template < class Float_type, int MaxUlps >
inline bool feq (Float_type lhs, Float_type rhs)
{
    // REMOVED if test as it (may) offer performance advantage [icpc: 15/17][g++: 22/20]
    // if (lhs == rhs) { return true; }

    // Make sure maxUlps is non-negative and small enough that the
    // default NAN won't compare as equal to anything.
    typename ftoi<Float_type>::int_type zr (1);
    zr <<= std::numeric_limits<Float_type>::digits;
    const union f_ {
        Float_type f; 
        typename ftoi<Float_type>::int_type i; 
        f_ (Float_type a) : f (a) {}
        f_ (typename ftoi<Float_type>::int_type a) : i(a) {}
    } lhsu (lhs), rhsu (rhs), zero (zr);
    // Make lhs_int lexicographically ordered as a twos-complement int
    const typename ftoi<Float_type>::int_type lhs_int (lhsu.f < 0 ? zero.i - lhsu.i : lhsu.i);

    // Make rhs_uint lexicographically ordered as a twos-complement int
    const typename ftoi<Float_type>::int_type rhs_int (rhsu.f < 0 ? zero.i - rhsu.i : rhsu.i);
    // std::cout << "fULPs = " << std::abs (lhs_int - rhs_int) << " for " << lhs_int << " and " << rhs_int << "\n";
    return MaxUlps >= (lhs_int > rhs_int ? lhs_int - rhs_int : rhs_int - lhs_int);
}

#ifndef USE_SFMT
/* --------------------------------------------------------------
  Use C++ TR1 random numbers
*/
typedef std::tr1::mersenne_twister< uint32_t, 32, 624, 397, 31, 0x9908b0df, 11, 7, 0x9d2c5680, 15, 0xefc60000, 18 > mt19937_t;

static mt19937_t &gen ()
{
#ifndef HAVE_UNIQPTR
    static std::auto_ptr< mt19937_t > gen_ = std::auto_ptr< mt19937_t >();
#else
    static std::unique_ptr< mt19937_t > gen_ = std::unique_ptr< mt19937_t >();
#endif

    if (NULL == gen_.get () )
    {
        gen_.reset (new mt19937_t);
    }
    return *gen_;
}

inline void seedrand (int seed)
{
    mt19937_t &gen_ = gen();
    gen_.seed(static_cast<unsigned long int>(seed));
    // Discard first 4000 numbers
    for (int i = 0; i < 4096; ++i) gen()();
}
inline double genrand ()
{
    mt19937_t &gen_ = gen();
    const double numerator = gen_() - gen_.min();
    const double denominator = gen_.max() - gen_.min();
    return (numerator/denominator);
}
#else
/* --------------------------------------------------------------
  Use SFMT mersenne twister
*/
#if ( USE_SFMT == 2 )
static dsfmt_t data_;
#endif

inline void seedrand (int seed)
{
#if ( USE_SFMT == 1 )
  init_gen_rand(static_cast<uint32_t>(seed));
#else
  dsfmt_init_gen_rand(&data_,static_cast<uint32_t>(seed));
#endif
}

inline double genrand ()
{

#if ( USE_SFMT == 1 )
  return genrand_res53();
#else
  return dsfmt_genrand_close1_open2(&data_) - 1.0;
#endif
}
#endif

/*  --------------------------------------------------------------
  Random shuffle of array
*/
void random_shuffle_int(int * array, int sz)
{
  std::random_shuffle(array, array + sz);
}


} // end namespace


int fuzzy_equal(double a, double b)
{
  return feq< double, 32 >(a,b);
}



// UUID GENERATOR
#include <uuid/uuid.h>

namespace 
{
void ionch_uuid_gen(char *const target)
{
  /* Assume target can hold 32 characters and does not need
     a nul terminator */
  uuid_t val;
  uuid_generate(val);
  char p[3];
  
  for (int i = 0; i != 16; ++i)
  {
    sprintf(&p[0], "%02X", val[i]);
    target[2*i]=p[0];
    target[2*i+1]=p[1];
  }
}
}



#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
extern "C"
{
#include <sys/utsname.h>
#include <sys/types.h>
}


namespace
{
/* --------------------------------------------------------------
  EXECUTION ENVIRONMENT QUERY FUNCTION

  This method gets information about the program environment.

  The key information obtained is:
  - random number library and version
  - cpu type (and cache information if available)
  - hostname
  - operating system

*/
  void write_host_info(char *buf, int* bufsize, int* ontermsig)
{
  class helper_
  {
  public:
    static std::string file_to_string (std::string dirname, const char* fname)
    {
      std::string result;
      dirname.append ("/");
      dirname.append (fname);
      std::ifstream ifs_(dirname.c_str ());
      if (ifs_) std::getline (ifs_, result);
      return result;
    }
    static int file_to_number (std::string dirname, const char* fname)
    {
      int result = -1;
      std::string s_;
      dirname.append ("/");
      dirname.append (fname);
      std::ifstream ifs_(dirname.c_str ());
      if (ifs_)
      {
        ifs_ >> result;
        if (ifs_)
        {
          ifs_ >> s_;
          if (not s_.empty ())
          {
            if ('K' == s_[0] or 'k' == s_[0]) result *= 1024;
            else if ('M' == s_[0] or 'm' == s_[0]) result *= 1024*1024;
          }
        }
      }
      return result;
    }
  };
  onterm = ontermsig;
  utsname uname_;
  uname (&uname_);
  int cache_max = -1;
  const char num_to_char[] = { "0123456789" };
  int cache_sz[3] = { -1, -1, -1 };
  int cache_lnsz[3] = { -1, -1, -1 };
  // RANDOM NUMBER GENERATOR USED:
  std::ostringstream os_;
#ifndef USE_SFMT
  os_ << " RANDOM NUMBER GENERATOR: Standard C++ implmentation\n";
#elif ( USE_SFMT == 1 )
  os_ << " RANDOM NUMBER GENERATOR: " << get_idstring() << "\n";
#elif ( USE_SFMT == 2 )
  os_ << " RANDOM NUMBER GENERATOR: " << dsfmt_get_idstring() << "\n";
#endif
  os_ << "------------------------------------------------------------------------\n";
  if (signal(SIGTERM, &on_term_signal) == SIG_ERR)
    {
      os_ << "Cannot handle SIGTERM!\n"; 
      os_ << "------------------------------------------------------------------------\n";
    }
  if (signal(SIGUSR2, &on_term_signal) == SIG_ERR)
    {
      os_ << "Cannot handle SIGUSR2!\n"; 
      os_ << "------------------------------------------------------------------------\n";
    }
  if (signal(SIGUSR1, &on_term_signal) == SIG_ERR)
    {
      os_ << "Cannot handle SIGUSR1!\n"; 
      os_ << "------------------------------------------------------------------------\n";
    }
// linux specific information
  for (int index = 0; index < 10; ++index)
  {
    std::string dirname ("/sys/devices/system/cpu/cpu0/cache/index_");
    dirname[dirname.size () - 1] = num_to_char[index];
    // std::cerr << "Checking directory \"" << dirname << "\"\n";
    struct stat st;
    if(stat (dirname.c_str (),&st) != 0)
    {
      // std::cerr << "Stat failed\n";
      break;
    }
    const std::string tmp_ (helper_::file_to_string (dirname, "type"));
    // std::cerr << "Got type \"" << tmp_ << "\"\n";
    if (std::string::npos != tmp_.find ("ata")
        or std::string::npos != tmp_.find ("nified"))
    {
      // Have data cache
      int level = helper_::file_to_number (dirname, "level") - 1;
      // std::cerr << "Found level \"" << level << "\"\n";
      if (0 <= level and level > 2) break; // Can`t handle more than three cache levels
      cache_max = std::max (cache_max, level);
      cache_lnsz[level] = helper_::file_to_number (dirname, "coherency_line_size");
      cache_sz[level] = helper_::file_to_number (dirname, "size");
      // std::cerr << "Found cache line \"" << cache_lnsz[level] << "\"\n";
      // std::cerr << "Found cache size \"" << cache_sz[level] << "\"\n";
    }
  }
  /* Print some machine information */
  os_ << " HOST MACHINE\n  NAME: " << uname_.nodename << "\n  TYPE: " << uname_.machine << "\n";
  if (cache_sz[0] != -1)
  {
    os_ << "  L1 CACHE (line,sz): " << cache_sz[0] << ", " << cache_lnsz[0] << "\n";
  }
  else
  {
    os_ << "  L1 CACHE (line,sz): unknown or not present\n";
  }
  if (cache_sz[1] != -1)
  {
    os_ << "  L2 CACHE (line,sz): " << cache_sz[1] << ", " << cache_lnsz[1] << "\n";
  }
  else
  {
    os_ << "  L2 CACHE (line,sz): unknown or not present\n";
  }
  if (cache_sz[2] != -1)
  {
    os_ << "  L3 CACHE (line,sz): " << cache_sz[2] << ", " << cache_lnsz[2] << "\n";
  }
  else
  {
    os_ << "  L3 CACHE (line,sz): unknown or not present\n";
  }
  os_ << " OPERATING SYSTEM\n  NAME: " << uname_.sysname << "\n  VERSION: " << uname_.version << "\n  RELEASE: " << uname_.release << "\n";
#ifdef _OPENMP
  os_ << " OpenMP cores: " << omp_get_max_threads() << "\n";
#else
  os_ << " Code not compiled with OpenMP support.\n";
#endif
  {
    std::string const result_ (os_.str());
    if (*bufsize > static_cast<int>(result_.size()))
      {
	*bufsize = result_.size();
	result_.copy(buf, *bufsize);
      }
    else
      {
	*bufsize = result_.size();
      }
  }

}
} // end namespace

/* --------------------------------------------------------------

*/
#ifdef USE_ATLAS
extern "C"
{
#include <atlas_buildinfo.h>
}

namespace
{
  void mathlib_version(char *const buf)
  {
    size_t len = snprintf(buf, 255, "ATLAS math library version %s: arch %s", ATL_VERS, ATL_ARCH);
    for (;len<256;++len) buf[len]=' ';
  }
} // end namespace
#endif

#ifdef USE_GSL
extern "C"
{
#include <gsl_version.h>
}
namespace 
{
  void mathlib_version(char *const buf)
  {
    size_t len = snprintf(buf, 255, "Gnu Scientific library %s", gsl_version);
    for (;len<256;++len) buf[len]=' ';
  }
} // end namespace
#endif


namespace
{

  // ------------------------------------------------------------
  // Method to convert string of floats into an array of floats
  // 
  //  \param strval : input string containing one or more floats
  //  \param a_array: output array of doubles
  //  \param a_size : the size of a_array
  //  \param a_cnt  : the output number of doubles read
  //  \param a_istt : error indicator, 0 if no errors
  //
  void read_floats(char const* strval, double *a_arry, int a_size, int * a_cnt, int * a_istt)
  {
    std::string input_buffer(strval);
    // Convert any 'd's to 'e's
    std::replace(input_buffer.begin(), input_buffer.end(), 'D', 'E');
    std::replace(input_buffer.begin(), input_buffer.end(), 'd', 'E');
    *a_istt = 0;
    {
      std::stringstream is(input_buffer);
      for (*a_cnt = 1; *a_cnt != a_size + 1; ++(*a_cnt))
	{
	  is >> *(a_arry + *a_cnt - 1);
	  if (not is)
	    {
	      // No more numbers for conversion
	      --(*a_cnt);
	      break;
	    }
	  if (is.eof()) break;
	}
    }
  }
  
}

#endif

