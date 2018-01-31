#include <exception>
#include <stdexcept>
using namespace std;

class myexception: public exception
{
  virtual const char* what() const throw()
    {
      return "Failure in ConfigParser code";
    }
} the_ex;
