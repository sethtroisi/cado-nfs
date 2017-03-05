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
    message(FATAL_ERROR "Need either a C++11 compiler, or boost libraries installed. You may simply unpack the boost headers at the root of the cado tree, that will probably do")
endif()
