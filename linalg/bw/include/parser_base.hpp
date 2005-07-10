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
};

class parser_base {
public:
	static const unsigned int buffer_size=32768;
private:
	char	buf[buffer_size];
	unsigned int pos;
	unsigned int size;
	std::stack<unsigned int> pos_stack;
protected:
	std::string string_version() const {
		return std::string(buf,pos,size);
	}
	void discard(unsigned int x) { pos=pos+x; if (pos>size) pos=size; }
	void save() { pos_stack.push(pos); }
	void restore() { pos=pos_stack.top(); pos_stack.pop(); }
	bool at_end() const { return pos >= size; }
	void require(bool t) const {
		if (!t) throw parse_error("glourg",buf,buf+pos);
	}
	void skipws() {
		for(;pos<size && isspace(buf[pos]);pos++);
	}
	unsigned int numlen() {
		unsigned int w;
		skipws();
		for(w=0;(pos+w)<size && isdigit(buf[pos+w]);w++);
		return w;
	}
	unsigned int numlen_sign() {
		unsigned int w;
		skipws();
		w=0;
		if (buf[pos+w]=='-') w++;
		for(;(pos+w)<size && isdigit(buf[pos+w]);w++);
		return w;
	}
	unsigned int num_in_cset(const std::string& s) {
		unsigned int w;
		skipws();
		for(w=0;(pos+w)<size && s.find(buf[pos+w]) != s.npos ;w++);
		return w;
	}
	std::string eat_cset(const std::string& s) {
		std::string r;
		skipws();
		for(;pos < size && s.find(buf[pos]) != s.npos ; pos++)
			r.append(1, buf[pos]);
		return r;
	}
	void error(const std::string& msg) {
		throw parse_error(msg,buf,buf+pos);
	}
	bool ptry(const char * s) {
		skipws();
		if (strncmp(buf+pos,s,MIN(size-pos,strlen(s)))==0) {
			pos+=strlen(s);
			return true;
		}
		return false;
	}
	bool ptry(char c) {
		skipws();
		if (pos < size && buf[pos]==c) { pos++; return true; }
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
		int sign=1;
		if (ptry('-')) sign=-1;
		int w=numlen();
		if (w==0) error("want a number");
		std::string s;
		s.append(buf+pos,w);
		pos+=w;
		mpz_class res;
		res.set_str(s,10);
		if (sign==-1) res=-res;
		return res;
	}
#endif
	int want_si() {
		int sign=1; 
		if (ptry('-')) sign=-1;
		int w=numlen();
		if (w==0) error("want a number");
		std::string s;
		s.append(buf+pos,w);
		pos+=w;
		int res;
		std::istringstream(s) >> res;
		return sign*res;
	}
	unsigned int want_ui() {
		int w=numlen();
		if (w==0) error("want a number");
		std::string s;
		s.append(buf+pos,w);
		pos+=w;
		unsigned int res;
		std::istringstream(s) >> res;
		return res;
	}
	parser_base(const char * obuf, unsigned int osize) 
		: pos(0), size(osize), pos_stack()
	{
		memcpy(buf,obuf,size);
		for(;size>=1 && buf[size-1]==0;size--);
	}
};

#endif	/* PARSER_BASE_HPP_ */
