#ifndef CONFIG_FILE_HPP_
#define CONFIG_FILE_HPP_

#include <map>
#include <string>
#include <sstream>

typedef std::map<std::string, std::string> config_file_t;

void read_config_file(config_file_t&, std::string const& f, bool mandatory);
void write_config_file(config_file_t const & cf, std::string const& f);

template <class T>
bool get_config_value(config_file_t const & cf, std::string& k, T & x)
{
    config_file_t::const_iterator ptr = cf.find(k);
    if (ptr == cf.end()) {
        return false;
    }

    std::istringstream  st(ptr->second);
    return (st >> x);
}

template <class T>
bool get_config_value(config_file_t const & cf, const char * k, T & x)
{
    std::string ks(k);
    return get_config_value(cf, ks, x);
}

template <class T>
void set_config_value(config_file_t & cf, std::string& k, T const & x)
{
    std::ostringstream  st;
    st << x;
    config_file_t::iterator ptr = cf.find(k);
    if (ptr == cf.end()) {
        cf.insert(make_pair(k, st.str()));
    } else {
        ptr->second = st.str();
    }
}
template <class T>
void set_config_value(config_file_t & cf, const char * k, T const & x)
{
    std::string ks(k);
    return set_config_value(cf, ks, x);
}

#endif	/* CONFIG_FILE_HPP_ */
