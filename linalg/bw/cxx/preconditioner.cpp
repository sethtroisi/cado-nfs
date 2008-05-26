#include "preconditioner.hpp"

#include <sstream>
#include "parsing_tools.hpp"
#include <iostream>


bool preconditioner_base::init(unsigned int i0, unsigned int i1, std::string const& fn)
{
	// using namespace boost;
	using namespace std;
	ifstream is(fn.c_str());
	if (!is.is_open()) {
		/* Not an error */
		cout << "// No preconditioner file " << fn << endl;
		return true;
	}
	comment_strip cs(is, "//");
	string s;
	combinations.clear();
	while (is) {
		cs.getline(s);
		if (!is)
			break;
		srclist_t src;
		dstlist_t dst;
		istringstream st(s);
		for( ; st ; st >> ws) {
			uint32_t idx;
			int32_t val;
			if (try_match(st, '>'))
				break;
			if (!(st >> idx >> wanted(':') >> noskipws >> val)) {
				is.setstate(ios::failbit);
				break;
			}
			src.push_back(make_pair(idx,val));
		}
		if (src.empty()) {
			is.setstate(ios::failbit);
			break;
		}
		unsigned int some_dst = 0;
		for(st >> skipws ; st ; st >> ws) {
			uint32_t idx;
			if (!(st >> idx)) {
				if (!st.eof()) {
					is.setstate(ios::failbit);
				}
				break;
			}
			some_dst++;
			if (idx >= i0 && idx < i1)
				dst.push_back(idx);
		}
		if (!some_dst) {
			is.setstate(ios::failbit);
			break;
		}
		if (!dst.empty()) {
			combinations.push_back(make_pair(src, dst));
		}
	}
	cout << "// Selected " << combinations.size()
		<< " preconditioning rules from file " << fn << endl;
	return is;
}
