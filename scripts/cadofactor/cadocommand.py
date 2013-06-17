import subprocess
import logging
import re
import cadoprograms
import cadologger

logger = logging.getLogger("Command")

class Command(object):
    ''' Represents a running subprocess
    
    The subprocess is started when the instance is initialised
    '''
    @staticmethod
    def _shellquote(s):
        ''' Quote a command line argument
        
        Currently does it the hard way: encloses the argument in single
        quotes, and escapes any single quotes that are part of the argument
        '''
        # If only characters that are known to be shell-safe occur, don't quote
        if re.match("^[a-zA-Z0-9+-_.:@/]*$", s):
            return s
        return "'" + s.replace("'", "'\\''") + "'"
    
    @staticmethod
    def is_same(a,b):
        """ Return true if a and b are the same object and not None """
        return (not a is None and a is b)
    
    @staticmethod
    def _open(fn, is_output):
        none = [None, subprocess.PIPE]
        mode = ["r", "w"]
        if fn is None:
            return none[int(is_output)]
        return open(fn, mode[int(is_output)])
    
    @staticmethod
    def _close_or_not(fd, ref):
        if not fd in [None, subprocess.PIPE, subprocess.STDOUT, ref]:
            fd.close()
    
    def __init__(self, program, *args, **kwargs):
        self.program = program
        progargs = self.program.make_command_array()
        
        # Each of stdin, stdout, stderr is one of these three:
        # - subprocess.PIPE, if program.get_std*() returned None
        # - subprocess.STDOUT for stderr if program.get_stdout() and
        #   program.get_stdin() return the same object
        # - a newly opened file handle which needs to be closed later
        self.stdin = self._open(self.program.get_stdin(), False)
        self.stdout = self._open(self.program.get_stdout(), True)
        if self.is_same(self.program.get_stdout(), self.program.get_stderr()):
            self.stderr = subprocess.STDOUT
        else:
            self.stderr = self._open(self.program.get_stderr(), True)
        
        self.child = subprocess.Popen(progargs, *args, stdin=self.stdin,
            stdout=self.stdout, stderr=self.stderr, **kwargs)
        cmdline = self.program.make_command_line()
        logger.cmd(cmdline, self.child.pid)
    
    def wait(self):
        ''' Wait for command to finish executing, capturing stdout and stderr 
        in output tuple '''
        (stdout, stderr) = self.child.communicate()
        if self.child.returncode != 0:
            logger.warning("Process with PID %d finished with return code %d",
                           self.child.pid, self.child.returncode)
        else:
            logger.debug("Process with PID %d finished successfully",
                         self.child.pid)
        if stdout:
            logger.debug("Process with PID %d stdout: %s", 
                         self.child.pid, stdout)
        if stderr:
            logger.debug("Process with PID %d stderr: %s", 
                         self.child.pid, stderr)
        self.returncode = self.child.returncode
        
        self._close_or_not(self.stdin, self.program.get_stdin())
        self._close_or_not(self.stdout, self.program.get_stdout())
        self._close_or_not(self.stderr, self.program.get_stderr())
        
        return (self.returncode, stdout, stderr)

class RemoteCommand(Command):
    def __init__(self, program, host, parameters, path_prefix, *args, **kwargs):
        # We use a make_command_line() instead of make_command_array() so that,
        # e.g., stdio redirection to files specified in program can be added
        # to the command line with and the redirection happens on the remote
        # host
        cmdline = program.make_command_line()
        self.prog = cadoprograms.SSH
        progparams = parameters.myparams(self.prog.get_config_keys(),
                                         path_prefix + [self.prog.name])
        ssh = self.prog([host, "env", "sh", "-c", Command._shellquote(cmdline)],
            progparams)
        super().__init__(ssh, *args, **kwargs)

class SendFile(Command):
    rsync_options = []
    def __init__(self, localfile, hostname, hostpath, parameters, *args,
                 **kwargs):
        if hostname == "localhost":
            target = hostpath
        else:
            target = hostname + ":" + hostpath
        rsync = cadoprograms.RSync([localfile, target], parameters)
        super().__init__(rsync, *args, **kwargs)

if __name__ == '__main__':
    import cadoparams

    parameters = {"long": True}
    program = cadoprograms.Ls("/", parameters)
    c = Command(program)
    (rc, out, err) = c.wait()
    if out:
        print("Stdout: " + str(out, encoding="utf-8"))
    if err:
        print("Stderr: " + str(err, encoding="utf-8"))
    del(c)

    program = cadoprograms.Ls("/", parameters, stdout = "ls.out")
    ssh_parameters = cadoparams.Parameters({"verbose": False})
    c = RemoteCommand(program, "localhost", ssh_parameters, [])
    (rc, out, err) = c.wait()
    print("Stdout: " + str(out, encoding="utf-8"))
    print("Stderr: " + str(err, encoding="utf-8"))
    del(c)

    c = SendFile("/users/caramel/kruppaal/ls.out", "quiche", "/tmp/foo", {})
    (rc, out, err) = c.wait()
    print("Stdout: " + str(out, encoding="utf-8"))
    print("Stderr: " + str(err, encoding="utf-8"))
    del(c)
