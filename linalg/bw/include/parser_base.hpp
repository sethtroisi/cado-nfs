#ifndef PARSER_BASE_HPP_
#define PARSER_BASE_HPP_

#ifdef	WITH_GMP
#include <gmp.h>
#include <gmpxx.h>
#endif

#include <string>
#include <stack>
#include <sstream>

struct parse_error {
	std::string msg;
	/*
	parse_error(const std::string& omsg, const char * buf, const char * ptr)
	{
		msg.assign("Parse error: ");
		msg+=omsg;
		msg+=" while parsing: \"";
		msg.append(buf,strlen(buf));
		msg+="\" ; exact position: \"";
		msg.append(ptr,strlen(ptr));
		msg+="\"";
	}
	*/
	parse_error(const std::string& omsg, const std::string& buf, unsigned int pos)
	{
		msg.assign("Parse error: ");
		msg+=omsg;
		msg+=" while parsing: \"";
		msg.append(buf);
		msg+="\" ; exact position: \"";
		msg.append(buf.substr(pos));
		msg+="\"";
	}
};

class parser_base {
public:
	/*
	static const unsigned int buffer_size=32768;
	*/
private:
	/*
	char	buf[buffer_size];
	*/
	const std::string& buf;
	unsigned int pos;
	std::stack<unsigned int> pos_stack;
protected:
	inline std::string string_version() const { return buf.substr(pos); }
	void discard(unsigned int x) {
		pos += x;
		if (pos > buf.size())
			pos = buf.size();
	}
	void save() { pos_stack.push(pos); }
	void restore() { pos = pos_stack.top(); pos_stack.pop(); }
	bool at_end() const { return pos >= buf.size(); }
	void require(bool t) const {
		if (!t) throw parse_error("glourg",buf,pos);
	}
	void skipws() {
		for(; pos < buf.size() && isspace(buf[pos]) ; pos++);
	}
	unsigned int numlen() {
		unsigned int w;
		skipws();
		for(w = 0 ; (pos+w) < buf.size() && isdigit(buf[pos+w]) ; w++);
		return w;
	}
	unsigned int numlen_sign() {
		unsigned int w;
		skipws();
		w=0;
		if (buf[pos+w] == '-') w++;
		for( ; (pos+w) < buf.size() && isdigit(buf[pos+w]) ; w++);
		return w;
	}
	unsigned int num_in_cset(const std::string& s) {
		std::string::size_type w;
		skipws();
		w = buf.find_first_not_of(s, pos);
		if (w == buf.npos)
			w = buf.size();
		return w - pos;
	}
	std::string eat_cset(const std::string& s) {
		std::string r;
		skipws();
		unsigned int w = num_in_cset(s);
		r = buf.substr(pos, pos + w);
		pos += w;
		return r;
	}
	void error(const std::string& msg) {
		throw parse_error(msg,buf,pos);
	}
	bool ptry(const char * s) {
		skipws();
		if (buf.substr(pos, strlen(s)) == s) {
			pos += strlen(s);
			return true;
		}
		return false;
	}
	bool ptry(char c) {
		skipws();
		if (pos < buf.size() && buf[pos]==c) { pos++; return true; }
		return false;
	}
	std::string want_cset(const std::string& s) {
		std::string r = eat_cset(s);
		if (s.empty()) {
			std::string msg("unexpected literal (want from ");
			(msg+=s)+=")";
			error(msg);
		}
		return r;
	}
	void want(const char * s) {
		if (!ptry(s)) {
			std::string msg("unexpected literal (want ");
			(msg+=s)+=")";
			error(msg);
		}
	}
	void want(char c) {
		if (!ptry(c)) {
			std::string msg("unexpected literal (want ");
			(msg+=c)+=")";
			error(msg);
		}
	}
#ifdef	WITH_GMP
	mpz_class want_mpz() {
		int sign = 1;
		if (ptry('-')) sign=-1;
		int w = numlen();
		if (w == 0) error("want a number");
		std::string s;
		s.append(buf.substr(pos, w));
		pos += w;
		mpz_class res;
		res.set_str(s,10);
		if (sign == -1) res = -res;
		return res;
	}
#endif
	int want_si() {
		int sign = 1; 
		if (ptry('-')) sign = -1;
		int w = numlen();
		if (w == 0) error("want a number");
		std::string s;
		s.append(buf.substr(pos, w));
		pos += w;
		int res;
		std::istringstream(s) >> res;
		return sign * res;
	}
	unsigned int want_ui() {
		int w = numlen();
		if (w == 0) error("want a number");
		std::string s;
		s.append(buf.substr(pos, w));
		pos += w;
		unsigned int res;
		std::istringstream(s) >> res;
		return res;
	}
	/*
	 * FIXME : the callers must be fixed to match what cab_rel_parser
	 * does.
	parser_base(const char * obuf, unsigned int osize) 
		: pos(0), size(osize), pos_stack()
	{
		memcpy(buf,obuf,size);
		for(;size>=1 && buf[size-1]==0;size--);
	}
	*/
	parser_base(const std::string& b) : buf(b), pos(0), pos_stack() {}
};

#endif	/* PARSER_BASE_HPP_ */
