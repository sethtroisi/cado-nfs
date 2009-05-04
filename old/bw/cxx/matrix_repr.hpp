#ifndef MATRIX_REPR_HPP_
#define MATRIX_REPR_HPP_

template<typename T> inline std::istream_iterator<T> beginof(std::istream& i)
{ return std::istream_iterator<T>(i); }
template<typename T> inline std::istream_iterator<T> endof(std::istream&)
{ return std::istream_iterator<T>(); }


#endif	/* MATRIX_REPR_HPP_ */
