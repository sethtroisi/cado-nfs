include(CheckCXXSourceCompiles)

set(tinker_with_glibc_internals "
#include <vector>

struct fred : public std::vector<int> {
    typedef std::vector<int> Base;
    fred(size_type n, int a) : Base(n, a) {}
    bool foo() {
        void * x = Base::_M_impl._M_start;
        void * y = Base::_M_impl._M_finish;
        void * z = Base::_M_impl._M_end_of_storage;
        return x == y || x == z;
    }
    Base::allocator_type * bar() { return &Base::_M_get_Tp_allocator(); };
};

int main() {
    fred F(12, 1);
}
")

CHECK_CXX_SOURCE_COMPILES("${tinker_with_glibc_internals}" HAVE_GLIBC_VECTOR_INTERNALS)
