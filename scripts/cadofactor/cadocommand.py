import subprocess
import logging
import re
import cadoprograms

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

    def __init__(self, program, *args, **kwargs):
        self.logger = logging.getLogger("Command")
        progargs = program.make_command_array()
        self.child = subprocess.Popen(progargs, *args, stdin=program.stdin,
            stdout=program.stdout, stderr=program.stderr, **kwargs)
        cmdline = program.make_command_line()
        self.logger.cmd(cmdline, self.child.pid)
    
    def wait(self):
        ''' Wait for command to finish executing, capturing stdout and stderr 
        in output tuple '''
        (stdout, stderr) = self.child.communicate()
        if self.child.returncode != 0:
            self.logger.error("Process with PID %d finished with return "
                              "code %d",
                              self.child.pid, self.child.returncode)
        else:
            self.logger.debug("Process with PID %d finished sucessfully",
                              self.child.pid)
        if stdout:
            self.logger.debug("Process with PID %d stdout: %s", 
                              self.child.pid, stdout)
        if stderr:
            self.logger.debug("Process with PID %d stderr: %s", 
                              self.child.pid, stderr)
        self.returncode = self.child.returncode
        return (self.returncode, stdout, stderr)

class RemoteCommand(Command):
    def __init__(self, program, host, parameters, stdin = None, 
                 stdout = subprocess.PIPE, stderr = subprocess.PIPE, **kwargs):
        # The remote command line. Need not be quoted as it is given to ssh
        # as an array element. 
        # We use a shell command line instead of an array so that, e.g., stdio
        # redirection to files specified in program can be added to the command
        # line with and the redirection happens on the remote host
        cmdline = program.make_command_line()
        # Hostname should likewise be quoted, even though we probably won't
        # use user or host names with shell meta characters
        ssh = cadoprograms.SSH([self._shellquote(host), cmdline], parameters)
        super().__init__(ssh, **kwargs)

# Things to consider: allow "user" parameter? Expect it as part of hostname?
# Expect port as part of hostname? If we have a user, copying to localhost
# may not be possible directly as the user running this script may not have 
# write privilege.
# We should allow switching between rsync and scp, possibly even (s)ftp?
class SendFile(Command):
    rsync="/usr/bin/rsync"
    rsync_options = []
    def __init__(self, localfile, hostname, hostpath, port = None, 
                 rsync_options = None, **kwargs):
        if hostname != "localhost":
            target = hostname + ":"
            if not port is None:
                target += str(port)
        target += hostpath
        copy_command = [self.rsync] + self.rsync_options \
            + [localfile, target]
        super().__init__(copy_command, **kwargs)

if __name__ == '__main__':
    parameters = {"long": True}
    program = cadoprograms.Ls("/", parameters)
    c = Command(program)
    (rc, out, err) = c.wait()
    print("Stdout: " + str(out, encoding="utf-8"))
    print("Stderr: " + str(err, encoding="utf-8"))
    del(c)

    program = cadoprograms.Ls("/", parameters, stdout = "ls.out")
    ssh_parameters = {"verbose": False}
    c = RemoteCommand(program, "localhost", ssh_parameters)
    (rc, out, err) = c.wait()
    print("Stdout: " + str(out, encoding="utf-8"))
    print("Stderr: " + str(err, encoding="utf-8"))
    del(c)

#    c = SendFile("WU", "quiche", "/tmp/foo")
#    (rc, out, err) = c.wait()
#    print("Stdout: " + str(out, encoding="utf-8"))
#    print("Stderr: " + str(err, encoding="utf-8"))
#    del(c)
