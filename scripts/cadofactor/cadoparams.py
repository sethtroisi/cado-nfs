#!/usr/bin/env python3

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
tasks.sieve.threads = 1

the sieving task will use the value 1 for the threads parameter, while other 
tasks use the value 2 (unless they specify a different value in their subtree).

Tasks run programs, and those have their own node in the parameter tree.
For example, with
threads = 2
tasks.sieve.las.threads = 1
the theads=1 parameter would apply only to the las program, but not to any 
other programs run during the sieving tasks, if there were any. The name of 
the node of a program is equal to the name of the binary executable.
"""

import os
import re
import abc
import logging
import cadologger

logger = logging.getLogger("Parameters")

def BoolParam(value):
    """
    >>> BoolParam(True)
    True
    >>> BoolParam(False)
    False
    >>> BoolParam("yes")
    True
    >>> BoolParam("no")
    False
    """
    if value is True or isinstance(value, str) and \
            value.lower() in ["yes", "true", "on", "1"]:
        return True
    elif value is False or isinstance(value, str) and \
            value.lower() in ["no", "false", "off", "0"]:
        return False
    else:
        raise ValueError("Could not parse '%s' as truth value" % value)

class Parameters(object):
    """ Class that stores parameters for cadofactor in hierarchical dictionaries
    """
    # TODO: Add a self.used dictionary with the same structure and keys as
    # self.data, but booleans as the values which indicate whether a parameter
    # has ever been accessed with myparams(). Add a function that tests that
    # all parameters have been accessed, so that a warning can be printed
    # about parameters in the parameter file that are not used by anything,
    # which might indicate a misspelling, etc.
    
    def __init__(self, *args, **kwargs):
        self.data = dict(*args, **kwargs)
        self._have_read_defaults = False
        self.key_types = {}
    
    def myparams(self, keys, path):
        ''' From the hierarchical dictionary params, generate a flat 
        dictionary with those parameters which are listed in keys and that 
        are found along path. 
        path is specified as a string with path segments separated '.',
        or as a list of path segments.
        
        If keys is a dictionary, then the dictionary key will be used as the
        parameter key, and its value will be used as the default value. The
        parameter value is converted to the same type as the dictionary value
        is. If the dictionary value is a class (such as int, str, or bool),
        then we assume that there is no default value and the key is mandatory;
        an error will be raised if it is not found in the parameter hierarchy.
        
        >>> d = {'a':'1','b':'2','c':'3','foo':{'a':'3'},'bar':{'a':'4','baz':{'a':'5'}}}
        >>> Parameters(d).myparams(keys=('a', 'b'), path = 'foo') == {'a': '3', 'b': '2'}
        True
        
        >>> Parameters(d).myparams(keys=('a', 'b'), path = 'bar.baz') == {'a': '5', 'b': '2'}
        True

        Test returning the default value of a parameter not provided in the
        parameter file
        >>> Parameters(d).myparams(keys={'d': 1}, path=[])
        {'d': 1}

        Test converting to the same type as the default value
        >>> Parameters(d).myparams(keys={'a': 1}, path='foo')
        {'a': 3}
        
        Test converting to an explicit type
        >>> Parameters(d).myparams(keys={'a': int}, path='foo')
        {'a': 3}

        Test converting to a non-mandatory explicit type
        >>> Parameters(d).myparams(keys={'a': [int], 'x': [int]}, path='foo')
        {'a': 3}

        Test converting if default value is bool
        >>> Parameters({'foo': 'yes'}).myparams(keys={'foo': False}, path=[])
        {'foo': True}
        
        Test converting if explicit type is bool
        >>> Parameters({'foo': 'yes'}).myparams(keys={'foo': bool}, path=[])
        {'foo': True}
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
        if isinstance(keys, dict):
            for key in keys:
                # A default of None means no default or type given: if the
                # parameter is in the parameter file, then store it in result
                # as a string, and if not, then don't store it in result
                if keys[key] is None:
                    continue
                # If only the type without default value is specified, then
                # the value must exist in the parameter file, and is converted
                # to the specified type
                if type(keys[key]) is type:
                    target_type = keys[key]
                    if not key in result:
                        logger.critical("Parameter %s not found under path %s",
                            key, joinpath)
                        raise(KeyError)
                elif type(keys[key]) is list:
                    # If a list is given as the default value, its first
                    # element must be a type and we interpret it as a non-
                    # mandatory type: if the parameter exists in the parameter
                    # file, it is cast to that type.
                    target_type = keys[key][0]
                    if not key in result:
                        continue
                else:
                    target_type = type(keys[key])
                    result.setdefault(key, keys[key])
                
                # Convert type
                result[key] = self._convert_one_type(splitpath, key,
                                                     result[key], target_type)
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
            if not d in source:
                return None
            assert isinstance(source[d], dict)
            source = source[d]
        return source
    
    def find(self, path, regex):
        source = self._get_subdict(path)
        if not source:
            source = {}
        pattern = re.compile(regex)
        result = [[l[0], l[1]] for l in self._recurse_iter(source, path)
            if pattern.search(l[1])]
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
            match = re.search(r"^(.*)\$\{(.*)\}(.*)$", value)
            if not match:
                break
            (prefix, varname, postfix) = match.groups()
            if not varname in os.environ:
                raise KeyError('Shell environment variable ${%s} referenced '
                               'in key %s is not defined (maybe not exported?)'
                               % (varname, fqn))
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
        while isinstance(value, str):
            match = re.search(r"^(.*)\$\((.*)\)(.*)$", value)
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
    
    @staticmethod
    def _cast_to_int(value):
        """ Return value cast to int
        
        If we can't cast to int directly, try going through float first to
        parse scientific notation
        
        >>> Parameters._cast_to_int("1")
        1
        >>> Parameters._cast_to_int("1x")
        Traceback (most recent call last):
        ValueError: invalid literal for int() with base 10: '1x'
        >>> Parameters._cast_to_int("1.5e4")
        15000
        >>> Parameters._cast_to_int("1.05e1")
        Traceback (most recent call last):
        ValueError: Value 1.05e1 cannot be converted to int without loss
        """
        try:
            return int(value)
        except ValueError as e:
            try:
                floatvalue = float(value)
                if float(int(floatvalue)) == floatvalue:
                    return int(floatvalue)
            except ValueError:
                # If that didn't work either, raise the original cast-to-int
                # exception
                raise e
        raise ValueError("Value %s cannot be converted to int without loss"
                         % value)
    
    def _convert_one_type(self, path, key, orig_value, datatype,
                          fatal_keytype=False):
        """ Convert orig_value to type datatype
        
        For datatype=None, return orig_value unchanged.
        For datatype=bool, use BoolParam for the conversion.
        For datatype=int, try converting to float first if necessary to parse
        scientific notation.
        Add datatype to key_types, and check for conflict.
        
        >>> p = Parameters({})
        
        >>> p._convert_one_type([], "foo", "1", int)
        1
        
        >>> p._convert_one_type([], "bar", "1", bool)
        True
        
        >>> p._convert_one_type([], "foo", "1x", int)
        Traceback (most recent call last):
        ValueError: Cannot convert value 1x for parameter foo to type int: invalid literal for int() with base 10: '1x'
        
        >>> p._convert_one_type([], "foo", "1", bool, fatal_keytype=True)
        Traceback (most recent call last):
        Exception: Conflicting type request for parameter foo: previously used with type int, now bool
        """
        param = ".".join(path + [key])
        # print ("Trying to convert param=%s, value=%s" % (param, orig_value))
        if datatype is None:
            return value

        if datatype is int:
            castfunction = self._cast_to_int
        # bool Type is special, we use BoolParam for the conversion
        elif datatype is bool:
            castfunction = BoolParam
        else:
            castfunction = datatype
        
        try:
            value = castfunction(orig_value)
            # print ("Converted param=%s, value=%r to %r, type %s" %
            #       (param, orig_value, value, datatype.__name__))
        except ValueError as err:
            raise ValueError("Cannot convert value %s for parameter %s to type "
                             "%s: %s" % (orig_value, param, datatype.__name__,
                                         err))

        if key in self.key_types:
            if not datatype is self.key_types[key]:
                if fatal_keytype:
                    raise Exception("Conflicting type request for parameter "
                             "%s: previously used with type %s, now %s" %
                             (param, self.key_types[key].__name__,
                             datatype.__name__))
                logger.error("Conflicting type request for parameter "
                             "%s: previously used with type %s, now %s",
                             param, self.key_types[key].__name__,
                             datatype.__name__)
        else:
            self.key_types[key] = datatype

        return value


    def parseline(self, line):
        line2 = line.split('#', 1)[0] if '#' in line else line
        line2 = line2.strip()
        if not line2:
            return (None, None)
        if not '=' in line2:
            raise Exception('Invalid line, missing "=": %s' % line)
        # Which one is worse?
        # (key, value) = re.match(r'(\S+)\s*=\s*(\S+)', line).groups()
        (key, value) = (s.strip() for s in line2.split('=', 1))
        return (key, value)

    def readparams(self, infile):
        """ 
        Read configuration file lines from infile, which must be an iterable.
        An open file handle, or an array of strings, work.
        
        >>> p = Parameters()
        >>> p.readparams(["tasks.sieve.rels_wanted = 1", \
                          "tasks.polyselect.degree=5", \
                          "tasks.polyselect.incr =60"])
        >>> p.data["tasks"]["sieve"]["rels_wanted"]
        '1'
        >>> p.myparams(["degree", "incr"], "tasks.polyselect") == \
        {'incr': '60', 'degree': '5'}
        True
        """
        for line in infile:
            line = line.strip('\n')
            (key ,value) = self.parseline(line)
            # print ("%s = %s # %s", (key ,value))
            if key is None:
                continue
            value = self.subst_env_var(key, value)
            self._insertkey(key, value)
        self._subst_references(self.data, [])

    def readfile(self, filename):
        """ Read parameters from a file """
        logger.debug("Reading parameter file %s", filename)
        with open(filename, "r") as handle:
            self.readparams(handle)

    def __str_internal__(self):
        ''' Returns all entries of the dictionary dic as key=sep strings
        in an array
        '''
        return ("%s = %s" % (".".join(path + [key]), value) for
                (path, key, value) in self)
    
    def __str__(self):
        r = self.__str_internal__()
        return "\n".join(r)


class UseParameters(metaclass=abc.ABCMeta):
    """ Mix-in class for objects that take parameters from the parameter file
    """
    @abc.abstractproperty
    def param_nodename(self):
        """ The name of this object's node in the parameter tree """
        pass
    
    @abc.abstractproperty
    def paramnames(self):
        # A list of parameter keywords which this task uses.
        # This is used for extracting relevant parameters from the parameter
        # hierarchical dictionary.
        # Sub-classes need to define a property 'paramnames' which returns a
        # list of parameters they accept, plus super()'s paramnames list
        pass

    @staticmethod
    def list_to_dict(a):
        if a is None:
            return {}
        elif isinstance(a, dict):
            return a.copy()
        else:
            return {k:None for k in a}

    @staticmethod
    def join_params(a, b):
        """ Join two dictionaries
        
        The values from the second take precedence in case of collision.
        Lists are converted to dictionaries whose keys map to None.

        >>> UseParameters.join_params(None, [2])
        {2: None}
        >>> UseParameters.join_params([1], None)
        {1: None}
        >>> UseParameters.join_params([1], [2]) == {1:None, 2:None}
        True
        >>> UseParameters.join_params({1:"a"}, [2]) == {1:"a", 2:None}
        True
        >>> UseParameters.join_params([1], {2:"a"}) == {1:None, 2:"a"}
        True
        """
        c = UseParameters.list_to_dict(a)
        c.update(UseParameters.list_to_dict(b))
        return c

    class MyParameters():
        """ Class that encapsules info on this node's parameters
        
        It stores a reference to the parameter dictionaries and info on this
        node's path and node name in the parameter tree.
        It can be initialised from a MyParameters or a Parameters instance;
        in the latter case the path prefix defaults to empty. If path_prefix
        is specified, it is always used as the path prefix.
        """
        def __init__(self, parent, name, path_prefix = None):
            self.name = name
            if isinstance(parent, UseParameters.MyParameters):
                self.parameters = parent.parameters
                self.path_prefix = parent.get_param_path()
            else:
                # Assume parent is a Parameters instance
                self.parameters = parent
                self.path_prefix = []
            if not path_prefix is None:
                self.path_prefix = path_prefix
        
        def get_param_path(self):
            if self.name is None:
                return self.path_prefix[:]
            else:
                return self.path_prefix + [self.name]
        
        def get_parameters(self):
            return self.parameters
        
        def myparams(self, keys, extrapath = None):
            path = self.get_param_path()
            if not extrapath is None:
                if isinstance(extrapath, str):
                    path += [extrapath]
                else:
                    path += extrapath
            return self.parameters.myparams(keys, path)
        
        def __str__(self):
            return "MyParameters: name = %s, prefix = %s" % (
                self.name, self.path_prefix)
        
    
    def __init__(self, *args, parameters, path_prefix = None,
                 **kwargs):
        self.parameters = UseParameters.MyParameters(
            parameters, self.param_nodename, path_prefix)
        super().__init__(*args, **kwargs)


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        import doctest
        doctest.testmod()
