import sys
import os
import re
import abc

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

class Parameters(object):
    # TODO: Add a self.used dictionary with the same structure and keys as
    # self.data, but booleans as the values which indicate whether a parameter
    # has ever been access with myparams(). Add a function that tests that
    # all parameters have been accessed, so that a warning can be printed
    # about parameters in the parameter file that are not used by anything,
    # which might indicate a misspelling, etc.
    def __init__(self, *args, **kwargs):
        self.data = dict(*args, **kwargs)
    
    def myparams(self, keys, path):
        ''' From the hierarchical dictionary params, generate a flat 
        dictionary with those parameters which are listed in keys and that 
        are found along path. 
        path is specified as a string with path segments separated '.'
        
        >>> d = {'a':1,'b':2,'c':3,'foo':{'a':3},'bar':{'a':4,'baz':{'a':5}}}
        >>> Parameters(d).myparams(keys=('a', 'b'), path = 'foo')
        {'a': 3, 'b': 2}
        
        >>> Parameters(d).myparams(keys=('a', 'b'), path = 'bar.baz')
        {'a': 5, 'b': 2}
        '''
        # path can be an array of partial paths, i.e., each entry can contain
        # one or more path segments separated by '.'. First join
        # them all, then split them again
        if isinstance(path, str):
            joinpath = path
        else:
            joinpath = '.'.join(path)
        splitpath = joinpath.split('.')
        
        source = self.data
        result = self._extract_by_keys(source, keys)
        for segment in splitpath:
            if not segment in source:
                break
            source = source[segment]
            result.update(self._extract_by_keys(source, keys))
        return result
    
    @staticmethod
    def _extract_by_keys(source, keys):
        return {key:source[key] for key in keys 
                if key in source and not isinstance(source[key], dict)}
    
    def __iter__(self):
        return self._recurse_iter(self.data, [])
    
    @staticmethod
    def _recurse_iter(source, path):
        for key in source:
            if isinstance(source[key], dict):
                for y in Parameters._recurse_iter(source[key], path + [key]):
                    yield y
            else:
                yield (path, key, source[key])
    
    def _get_subdict(self, path):
        source = self.data
        for d in path:
            assert d in source
            assert isinstance(source[d], dict)
            source = source[d]
        return source
    
    def find(self, path, regex):
        source = self._get_subdict(path)
        result = []
        pattern = re.compile(regex)
        for l in self._recurse_iter(source, path):
            if pattern.search(l[1]):
                result.append([l[0], l[1]])
        return result
    
    def _insertkey(self, path, value):
        ''' path is a path with segments delimited by '.' or an 
        array of pieces of the path, 
        value is inserted in the hierarchical parameters dictionary at the 
        location specified by path
        Keys overwrite previously existing keys, but a conflict between a key
        and a sub-dictionary causes a KeyError (similar to open('filepath','w')
        when filepath exists as a subdirectory).
        
        >>> p=Parameters({'a':1, 'b':2, 'foo':{'a':3}, 'bar':{'a':4}})
        >>> p._insertkey('c', 3)
        >>> p.data == {'a': 1, 'b': 2, 'c': 3, 'bar': {'a': 4}, 'foo': {'a': 3}}
        True
        >>> p._insertkey('bar.c', 5)
        >>> p.data == {'a': 1, 'b': 2, 'c': 3, 'bar': {'a': 4, 'c': 5}, 'foo': {'a': 3}}
        True
        >>> p._insertkey('bar.baz.c', 6)
        >>> p.data == {'a': 1, 'b': 2, 'c': 3, 'bar': {'a': 4, 'c': 5, 'baz': {'c': 6}}, 'foo': {'a': 3}}
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
            joinpath = '.'.join(path)
        splitpath = joinpath.split('.')
        key = splitpath.pop()
        dest = self.data
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
    
    @staticmethod
    def subst_env_var(fqn, value):
        """ Substitute strings like '${HOME}' in value by the corresponding
        shell environment variables
        """
        while True:
            match = re.search("^(.*)\$\{(.*)\}(.*)$", value)
            if not match:
                break
            (prefix, varname, postfix) = match.groups()
            if not varname in os.environ:
                raise KeyError('Shell environment variable ${%s} referenced '
                               'in key %s is not defined' % (varname, fqn))
            value = prefix + os.environ[varname] + postfix
        return value
    
    def _subst_reference(self, path, key, value):
        """ Substitute strings like '$(somekey)' in a value by the value of
        "somekey" found along the current path, e.g.,
        foo.bar.k = $(m)
        foo.m = 5
        k = $(m)
        m = 3
        results in foo.bar.k = 5 and k = 3
        """
        while True:
            match = re.search("^(.*)\$\((.*)\)(.*)$", value)
            if not match:
                break
            (prefix, varname, postfix) = match.groups()
            if key == varname:
                raise KeyError("Self-referential substitution $(%s) in key %s")
            result = self.myparams([varname], path)
            if not result:
                raise KeyError('Key $(%s) referenced in key %s is not defined'
                               % (varname, '.'.join(path + [key])))
            value = prefix + result[varname] + postfix
        return value
    
    def _subst_references(self, dic, path):
        for key in dic:
            if isinstance(dic[key], dict):
                self._subst_references(dic[key], path + [key])
            else:
                dic[key] = self._subst_reference(path, key, dic[key])
    
    def _readfile(self, infile):
        """ 
        Read configuration file lines from infile, which must be an iterable.
        An open file handle, or an array of strings, work.
        
        >>> p = Parameters()
        >>> p._readfile(DEFAULTS)
        >>> p.data["tasks"]["parallel"]
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
            value = self.subst_env_var(key, value)
            self._insertkey(key, value)
        self._subst_references(self.data, [])
    
    @staticmethod
    def __str_internal__(dic, path):
        ''' Returns all entries of the dictionary dic as key=sep strings
        in an array
        '''
        result = []
        for key in dic:
            if not isinstance(dic[key], dict):
                result.append(path + key + " = " + dic[key])
        for key in dic:
            if isinstance(dic[key], dict):
                result = result + Parameters.__str_internal__(dic[key], 
                                                              path + key + '.')
        return result
    
    def __str__(self):
        r = Parameters.__str_internal__(self.data, "")
        return "\n".join(r)

class UseParameters(object, metaclass=abc.ABCMeta):
    @abc.abstractproperty
    def name(self):
        pass
    
    def __init__(self, parameters, path_prefix):
        self.parameters = parameters
        self.path_prefix = path_prefix
        self.parampath = path_prefix + [self.name]
    
    def myparams(self, keys, extrapath = None):
        path = self.parampath
        if not extrapath is None:
            if isinstance(extrapath, str):
                path += [extrapath]
            else:
                path += extrapath
        return self.parameters.myparams(keys, path)
    
    def get_param_prefix(self):
        return self.path_prefix
    
    def get_param_path(self):
        return self.parampath
    
    def get_parameters(self):
        return self.parameters
    

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
    "tasks.linalg.hosts = foo",
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
