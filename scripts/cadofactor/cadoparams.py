import sys

"""
Parameter file format

Parameters for tasks and programs follow a hierarchical namespace, a tree 
similar to a directory structure, but with segments separated by the period 
character: "."
E.g., the parameters of the task "foo" live under
tasks.foo
Any parameters specified in the path to a node are considered when looking 
for a parameter in a node; parameters late in the path take precedence. 
Hence with:

threads = 2
tasks.sieving.threads = 1

the sieving task will use the value 1 for the threads parameter, while other 
tasks use the value 2 (unless they specify a different value in their subtree).

Tasks run programs, and those have their own node in the parameter tree.
For example, with
threads = 2
tasks.sieving.las.threads = 1
the theads=1 parameter would apply only to the las program, but not to any 
other programs run during the sieving tasks, if there are any. The name of 
the node of a program is equal to the name of the binary executable.
"""

class Parameters(dict):
    def __init__(self, *args, **kwargs):
        super(Parameters, self).__init__(*args, **kwargs)
        self._separator = '.'
    
    def set_sep(self, separator):
        if self:
            raise Exception("Can't change separator on non-empty Parameters "
                            "instance")
        self._separator = separator
    
    def get_sep(self):
        return self._separator
    
    def myparams(self, keys, path):
        ''' From the hierarchical dictionary params, generate a flat 
        dictionary with those parameters which are listed in keys and that 
        are found along path. 
        path is specified as a string with path segments separated by the 
        separator
        
        >>> d = {'a':1,'b':2,'c':3,'foo':{'a':3},'bar':{'a':4,'baz':{'a':5}}}
        >>> Parameters(d).myparams(keys=('a', 'b'), path = 'foo')
        {'a': 3, 'b': 2}
        
        >>> Parameters(d).myparams(keys=('a', 'b'), path = 'bar.baz')
        {'a': 5, 'b': 2}
        '''
        result = {}
        source = self
        idx = 0
        # path can be an array of partial paths, i.e., each entry can contain
        # one or more path segments separated by the separator. First join
        # them all, then split them again
        if isinstance(path, str):
            joinpath = path
        else:
            joinpath = self.get_sep().join(path)
        splitpath = joinpath.split(self.get_sep())
        while source:
          tomerge = {key:source[key] for key in keys 
                     if key in source and not isinstance(source[key], dict)}
          result.update(tomerge)
          if idx < len(splitpath) and splitpath[idx] in source:
            source = source[splitpath[idx]]
            idx = idx + 1
          else:
            source = None
        return result
    
    def _insertkey(self, path, value):
        ''' path is a path with segments delimited by the separator or an 
        array of pieces of the path, 
        value is inserted in the hierarchical parameters dictionary at the 
        location specified by path
        Keys overwrite previously existing keys, but a conflict between a key
        and a sub-dictionary causes a KeyError (similar to open('filepath','w')
        when filepath exists as a subdirectory).
        
        >>> p=Parameters({'a':1, 'b':2, 'foo':{'a':3}, 'bar':{'a':4}})
        >>> p._insertkey('c', 3)
        >>> p == {'a': 1, 'b': 2, 'c': 3, 'bar': {'a': 4}, 'foo': {'a': 3}}
        True
        >>> p._insertkey('bar.c', 5)
        >>> p == {'a': 1, 'b': 2, 'c': 3, 'bar': {'a': 4, 'c': 5}, 'foo': {'a': 3}}
        True
        >>> p._insertkey('bar.baz.c', 6)
        >>> p == {'a': 1, 'b': 2, 'c': 3, 'bar': {'a': 4, 'c': 5, 'baz': {'c': 6}}, 'foo': {'a': 3}}
        True
        >>> p._insertkey('bar.baz.c.d', 6)
        Traceback (most recent call last):
        KeyError: 'Subdirectory c already exists as key'
        >>> p._insertkey('bar', 6)
        Traceback (most recent call last):
        KeyError: 'Key bar already exists as subdictionary'
        '''
        
        if isinstance(path, str):
            joinpath = path
        else:
            joinpath = self.get_sep().join(path)
        splitpath = joinpath.split(self.get_sep())
        key = splitpath.pop()
        dest = self
        for segment in splitpath:
            if segment in dest.keys():
                if not isinstance(dest[segment], dict):
                    raise KeyError('Subdirectory %s already exists as key' 
                                   % segment)
                dest = dest[segment]
            else:
                dest[segment] = {}
                dest = dest[segment]
        if key in dest.keys() and isinstance(dest[key], dict):
            raise KeyError('Key %s already exists as subdictionary' % key)
        dest[key] = value
    
    def _readfile(self, infile):
        """ 
        Read configuration file lines from infile, which must be an iterable.
        An open file handle, or an array of strings, work.

        >>> p = Parameters()
        >>> p._readfile(DEFAULTS)
        >>> p["tasks"]["parallel"]
        '0'
        >>> p.myparams(["degree", "incr", "parallel"], "tasks.polyselect") == \
        {'parallel': '0', 'incr': '60', 'degree': '5'}
        True
        """
        for line in infile:
            line2 = line.split('#', 1)[0].strip()
            if not line2:
               continue
            if not '=' in line2:
                raise Exception('Invalid line, missing "=": %s' % line)
            # Which one is worse?
            # (key, value) = re.match(r'(\S+)\s*=\s*(\S+)', line).groups()
            (key, value) = (s.strip() for s in line2.split('=', 1))
            self._insertkey(key, value)
    
    @staticmethod
    def __str_internal__(dic, sep, path):
        ''' Returns all entries of the dictionary dic as key=sep strings
        in an array
        '''
        result = []
        for key in dic:
            if not isinstance(dic[key], dict):
                result.append(path + key + " = " + dic[key])
        for key in dic:
            if isinstance(dic[key], dict):
                result = result + Parameters.__str_internal__(dic[key], sep, 
                                                              path + key + sep)
        return result
    
    def __str__(self, sep = None):
        r = Parameters.__str_internal__(self, self.get_sep(), "")
        return "\n".join(r)
    
DEFAULTS = (
    "logfile = cado.log",
    "tasks.parallel = 0",
    "tasks.niceness = 10",
    
    "tasks.polyselect.degree = 5",
    "tasks.polyselect.lq = 1",
    "tasks.polyselect.nq = 1000",
    "tasks.polyselect.incr = 60",
    "tasks.polyselect.admin = 0",
    "tasks.polyselect.adrange = 1e7",
    "tasks.polyselect.delay  = 120",
    "tasks.polyselect.maxnorm = 1e9",
    
    "tasks.sieve.rlim = 8000000",
    "tasks.sieve.alim = 8000000",
    "tasks.sieve.lpbr = 29",
    "tasks.sieve.lpba = 29",
    "tasks.sieve.mfbr = 58",
    "tasks.sieve.mfba = 58",
    "tasks.sieve.rlambda = 2.3",
    "tasks.sieve.alambda = 2.3",
    "tasks.sieve.I = 13",
    "tasks.sieve.qmin = 12000000",
    "tasks.sieve.qrange = 1000000",
    "tasks.sieve.checkrange = 1",
    "tasks.sieve.firstcheck = 1",
    "tasks.sieve.delay = 120",
    "tasks.sieve.niceness = 19",
    "tasks.sieve.keeprelfiles = 0",
    "tasks.sieve.sieve_max_threads = 2",
    "tasks.sieve.poly_max_threads = 1",
    "tasks.sieve.ratq = 0",

    # filtering
    "tasks.purge.skip = -1", # should be about bwc_mn - 32
    "tasks.purge.keep = 208", # should be 160 + #ideals <= FINAL_BOUND 
                              # (cf purge.c)
    "tasks.purge.nslices_log = 1",
    "tasks.purge.filterlastrels = 1",
    
    "tasks.merge.skip = -1", # should be about bwc_mn - 32
    "tasks.merge.forbw = 3",
    "tasks.merge.coverNmax = 100",
    "tasks.merge.ratio = 1.5",
    "tasks.merge.keep = -1", # should be 128 + skip
    "tasks.merge.maxlevel = 15",

    # linalg
    "tasks.linalg.algo = bwc",
    "tasks.linalg.threads = 2",
    "tasks.linalg.mpi = 0",
    "tasks.linalg.hosts = """,
    "tasks.linalg.bwc.interval = 1000",
    "tasks.linalg.bwc.mm_impl = 'bucket'",
    "tasks.linalg.bwc.interleaving = 0",
    # bwc_mn should be 64 or 128
    "tasks.linalg.bwc.mn = 64",
    # shuffled product is expected to be better in most cases, at least
    # when we use MPI. Since it is the preferred communication algorithm
    # for large runs, we prefer to force its use also for mid-range
    # examples.
    "tasks.linalg.bwc.shuffled_product = 1",

    # characters
    "tasks.sqrt.nkermax = 30",
    "tasks.sqrt.nchar = 50",
    "tasks.sqrt.nthchar = 2"
    )



if __name__ == "__main__":
    import doctest
    doctest.testmod()
