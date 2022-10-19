// If parameter is not true, test fails
#define BEGIN_TEST { std::cout << "Starting test " <<  __FUNCTION__ << std::endl; }
#define END_TEST { std::cout << "Ending test " <<  __FUNCTION__ << std::endl; }
#define IS_TRUE(x) { if (!(x)) { std::cout << __FUNCTION__ << " failed on line " << __LINE__ << std::endl; throw "TEST FAILED";} }
#define IS_FALSE(x) { if ((x)) { std::cout << __FUNCTION__ << " failed on line " << __LINE__ << std::endl; throw "TEST FAILED";} }
#define IS_EQUAL(x, expected) { if (x != expected) { std::cout << __FUNCTION__ << " failed equals on line " << __LINE__ <<  std::endl; throw "TEST FAILED"; } }