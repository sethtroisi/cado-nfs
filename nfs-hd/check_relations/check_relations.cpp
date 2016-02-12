//Usage: ./check_relations <file_relation> [<output>, <output> <poly.poly> <lpb>
//<error>]

#include <stdio.h>
#include <cassert>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <NTL/ZZXFactoring.h>

using namespace std;
using namespace NTL;

// Read a polynomial in the correct form.
void read_poly(ZZX& a, string poly)
{
  stringstream sline;
  replace(poly.begin(), poly.end(), ',', ' ');
  sline.str(string());
  sline.clear();
  sline << "[" << poly << "]";
  sline >> a;
}

// Read a polynomial inside a relation.
void read_poly_relation(ZZX& a, string line)
{
  string split;
  stringstream sline;
  sline.str(line);
  getline(sline, split, ':');
  read_poly(a, split);
}

// Read the file *.poly.
void read_poly_file(Vec < ZZX >& f, char * path_poly)
{
  ifstream file_poly (path_poly);
  string line;
  Pair < ZZX, unsigned int > p_tmp;
  Vec< Pair< ZZX, unsigned int > > v_tmp;
  ZZX a;
  while (getline(file_poly, line)) {
    if (line[0] == 'p' && line[1] == 'o' && line[2] == 'l' && line[3] == 'y') {
      p_tmp.b = (int) (line[4] - '0');
      char * split;
      char * linec = (char * ) malloc(sizeof(char) * line.size() + 1);
      copy(line.begin(), line.end(), linec);
      linec[line.size()] = '\0';
      split = strtok(linec, ":");
      split = strtok(NULL, ":");
      line.assign(split);
      line.erase(0, 1);
      read_poly(a, line);
      p_tmp.a = a;
      v_tmp.append(p_tmp);
      free(linec);
    }
  }
  file_poly.close();

  // Stupid sort
  long n = v_tmp.length();
  do {
    long n_tmp = 0;
    for (long i = 1; i < n; i++) {
      if (v_tmp[i - 1].b > v_tmp[i].b) {
        p_tmp = v_tmp[i - 1];
        v_tmp[i - 1] = v_tmp[i];
        v_tmp[i] = p_tmp;
        n_tmp = i;
      }
    }
    n = n_tmp;
  } while (n != 0);

  for (long i = 0; i < v_tmp.length(); i++) {
    f.append(v_tmp[i].a);
  }
}

//Read the lpb.
void read_lpb(Vec <ZZ>& lpb, char * lpbs, long length)
{
  char * split;
  for (long i = 0; i < length; i++) {
    split = strtok(lpbs, ",");
    lpb.append(power2_ZZ(atol(split)));
  }
}

// Convert char of hex to dec.
long hexDec(char c)
{
#ifndef NDEBUG
  assert(c == '0' ||
      c == '1' ||
      c == '2' ||
      c == '3' ||
      c == '4' ||
      c == '5' ||
      c == '6' ||
      c == '7' ||
      c == '8' ||
      c == '9' ||
      c == 'a' ||
      c == 'b' ||
      c == 'c' ||
      c == 'd' ||
      c == 'e' ||
      c == 'f');
#endif // NDEBUG

  if (c == 'a'){
    return 10;
  } else if (c == 'b'){
    return 11;
  } else if (c == 'c') {
    return 12;
  } else if (c == 'd') {
    return 13;
  } else if (c == 'e') {
    return 14;
  } else if (c == 'f') {
    return 15;
  } else {
    return (long) (c - '0');
  }
}

// Convert a hex string in ZZ.
ZZ hexZZ(string hex)
{
  ZZ z;
  z = 0;
  ZZ base;
  base = 1;
  for (size_t i = hex.size(); i-- > 0 ;) {
    z = z + hexDec(hex[i]) * base;
    base = base * 16;
  }
  return z;
}

// TODO: Maybe nb_facto not useful.
// Extract factorizations from relation.
void extract_factorisation(Vec < Vec < ZZ > >& factorizations, string line,
    long nb_facto)
{
  string facto;
  Vec < string > factors;
  char * split;
  char * linec = (char * ) malloc(sizeof(char) * line.size() + 1);
  copy(line.begin(), line.end(), linec);
  linec[line.size()] = '\0';
  split = strtok(linec, ":");
  for (long i = 0; i < nb_facto; i++) {
    split = strtok (NULL, ":");
    facto.assign(split);
    factors.append(facto);
  }
  free(linec);

  for (long i = 0; i < factors.length(); i++) {
    Vec < ZZ > v_tmp;
    linec = (char * ) malloc(sizeof(char) * factors[i].size() + 1);
    copy(factors[i].begin(), factors[i].end(), linec);
    linec[factors[i].size()] = '\0';
    split = strtok(linec, ",");
    facto.assign(split);
    v_tmp.append(hexZZ(facto));
    while (split != NULL) {
      split = strtok(NULL, ",");
      if (split != NULL) {
        facto.assign(split);
        v_tmp.append(hexZZ(facto));
      }
    }
    factorizations.append(v_tmp);
    free(linec);
  }
}

void compute_factorization(ZZ& factorization, Vec <ZZ> factor)
{
  factorization = 1;
  for (long i = 0; i < factor.length(); i++) {
    factorization *= factor[i];
  }
}

unsigned int is_smooth(Vec <ZZ> factor, ZZ lpb) {
  for (long i = 0; i < factor.length(); i++) {
    if (factor[i] > lpb) {
      return 0;
    }
  }
  return 1;
}

int main(int argc, char ** argv)
{
  ZZX a;

  //Output and error files.
  ofstream filew;
  ofstream fileerr;

  //Input file.
  ifstream filer (argv[1]);

  //Polynomials f read from file_poly.
  Vec <ZZX> f;

  //Large prime bounds.
  Vec <ZZ> lpb;

  if (argc == 2) {
    char * filew_name = (char *) malloc(sizeof(char) * strlen(argv[1]) + 7);
    strcpy(filew_name, argv[1]);
    strcat(filew_name, ".output");
    filew.open(filew_name);
    free(filew_name);
  } else if (argc == 3) {
    filew.open(argv[2]);
  } else if (argc == 6) {
    filew.open(argv[2]);
    read_poly_file(f, argv[3]);
    read_lpb(lpb, argv[4], f.length());
    fileerr.open(argv[5]);
  } else {
    cerr << "Usage error\n";
    cerr << "./irreducible <file_relation> [<output>, <output> <poly.poly> "
      "<lpb> <error>\n";
  }

  //Number of relations.
  ZZ nb_rel;
  nb_rel = 0;
  ZZ nb_not_rel;
  nb_not_rel = 0;
  ZZ nb_err_rel;
  nb_err_rel = 0;

  string line;

  //Needed to factor a ZZX.
  Vec< Pair< ZZX, long > > factors;
  ZZ c;

  while (getline(filer, line)) {
    if (line[0] != '#') {
      read_poly_relation(a, line);
      assert(content(a) == 1);
      assert(LeadCoeff(a) > 0);
      assert(deg(a) > 0);
      //TODO: if there exists a irreducible test, use it.
      factor(c, factors, a);
#ifdef NDEBUG
      assert(c == 1);
#endif // NDEBUG

      if (factors.length() == 1) {
        if (argc == 2 || argc == 3) {
          filew << line << "\n";
          nb_rel++;
        } else if (argc == 6) {
          unsigned int error = 0;
          Vec < Vec < ZZ > > factorization;
          extract_factorisation(factorization, line, f.length());
          for (long j = 0; j < f.length(); j++) {
            ZZ norm;
            norm = abs(resultant(a, f[j]));
            ZZ test;
            compute_factorization(test, factorization[j]);

            if (norm != test || !is_smooth(factorization[j], lpb[j])) {
              fileerr << line << "\n";
              error = 1;
              nb_err_rel++;
              break;
            }
          }
          if (!error) {
            filew << line << "\n";
            nb_rel++;
          }
        }
      } else {
        nb_not_rel++;
      }
    }
  }

  filer.close();
  filew.close();
  if (argc == 6) {
    fileerr.close();
  }

  cout << "Number of relations: " << nb_rel << ".\n";
  cout << "Number of false relations: " << nb_not_rel << ".\n";
  if (nb_err_rel != 0) {
    cout << "Number of relations with error: " << nb_err_rel << "\n";
  }

  return 0;
}
