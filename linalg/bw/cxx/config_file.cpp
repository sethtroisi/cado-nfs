#include <string>
#include <fstream>
#include <sstream>
#include "config_file.hpp"
#include "manu.h"

using namespace std;

void read_config_file(config_file_t& cf, std::string const& f)
{
    ifstream is;
    is.open(f.c_str());
    if (!is.is_open()) {
        BUG();
    }
    for(;!is.eof();) {
        string s;
        getline(is, s);
        if (is.eof())
            break;
        string key;
        string value;
        unsigned int c;
        c = s.find_first_of(" :=");
        {
            istringstream st(s.substr(0,c));
            st >> key >> ws;
            BUG_ON(!st.eof());
        }
        c = s.find_first_not_of(" :=", c);
        {
            istringstream st(s.substr(c));
            st >> value >> ws;
            BUG_ON(!st.eof());
        }

        cf.insert(make_pair(key, value));
    }
}
