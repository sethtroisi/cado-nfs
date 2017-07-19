include(CheckCXXSourceCompiles)
CHECK_CXX_SOURCE_COMPILES("
    #include <memory>
    struct foo {
        std::shared_ptr<int> a;
    };
    int main()
    {
        foo f;
        f.a = std::make_shared<int>();
        return 0;
    }
    " HAVE_STD_SHARED_PTR)

if (NOT HAVE_STD_SHARED_PTR)
CHECK_CXX_SOURCE_COMPILES("
    #include <boost/shared_ptr.hpp>
    #include <boost/make_shared.hpp>
    namespace std {
        using boost::shared_ptr;
        using boost::make_shared;
    }
    struct foo {
        std::shared_ptr<int> a;
    };
    int main()
    {
        foo f;
        f.a = std::make_shared<int>();
        return 0;
    }
    " HAVE_BOOST_SHARED_PTR)
endif()

if (NOT HAVE_STD_SHARED_PTR AND NOT HAVE_BOOST_SHARED_PTR)
    message(FATAL_ERROR "Need either a C++11 compiler, or boost libraries installed. You may simply unpack the boost headers at the root of the cado tree, that will probably do. Also, most OS distributions include packages for boost, which can be easily installed (names may vary, e.g.  libboost-dev or boost-devel")
endif()

# See bug #21397
CHECK_CXX_SOURCE_COMPILES("
#include <memory>
class nocopy {
        protected:
                nocopy() {}
                ~nocopy() {}
        private:
                nocopy(const nocopy&);
                nocopy& operator=(const nocopy&);
};

class foo: private nocopy {
        int val;
        public:
        foo(int val) : val(val) {}
        int value() const { return val; }
};

std::shared_ptr<foo> blah()
{
        const int val = 65536;
        return std::make_shared<foo>(val);
}
" HAS_NOT_BUG_21397)
if(NOT HAS_NOT_BUG_21397)
    message(FATAL_ERROR "Error, see bug #21397")
endif()
