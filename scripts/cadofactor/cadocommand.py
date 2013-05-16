import subprocess
import cadologger

class Command(object):
    ''' Represents a running subprocess
    
    The subprocess is started when the instance is initialised
    '''
    def __init__(self, args, stdin = None, stdout = subprocess.PIPE, 
                 stderr = subprocess.PIPE, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.logger = cadologger.Logger()
        # Convert args array to a string for printing if necessary
        if isinstance(self.args, str):
            self.cmdline = self.args
        else:
            self.cmdline = " ".join(self.args)
        self.child = subprocess.Popen(
            self.args, stdin=stdin, stdout=stdout, stderr=stderr, 
            **self.kwargs
            )
        self.logger.cmd(self.cmdline, self.child.pid)
    
    def wait(self):
        ''' Wait for command to finish executing, capturing stdout and stderr 
        in output tuple '''
        (self.stdout, self.stderr) = self.child.communicate()
        if self.child.returncode != 0:
            self.logger.error("Process with PID %d finished with return "
                              "code %d",
                              self.child.pid, self.child.returncode)
        if self.stdout:
            self.logger.debug("Process with PID %d stdout: %s", 
                              self.child.pid, self.stdout)
        if self.stderr:
            self.logger.debug("Process with PID %d stderr: %s", 
                              self.child.pid, self.stderr)
        self.returncode = self.child.returncode
        return (self.returncode, self.stdout, self.stderr)

class RemoteCommand(Command):
    ssh="/usr/bin/ssh"
    # Options can be overridden with ssh_options to __init__().
    # An option set to None is skipped which allows un-setting these
    ssh_options = {
        "ConnectTimeout": 30,
        "ServerAliveInterval": 10,
        "PasswordAuthentication": "no"
    }
    def __init__(self, command, host, port = None, ssh_options = None, 
                 **kwargs):
        # ssh_command is list of strings, like argv
        ssh_command = [self.__class__.ssh]
        options = self.__class__.ssh_options.copy()
        if not ssh_options is None:
            options.update(ssh_options)
        if not port is None:
            ssh_command += ["-p", str(port)];
        for (opt, val) in options.items():
            if not val is None:
                ssh_command += ["-o", opt + "=" + str(val)]
        ssh_command.append(host)
        if isinstance(command, str):
            ssh_command.append(command)
        else:
            ssh_command += command
        super().__init__(ssh_command, **kwargs)

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
        copy_command = [self.__class__.rsync] + self.__class__.rsync_options \
            + [localfile, target]
        super().__init__(copy_command, **kwargs)

if __name__ == '__main__':
    c = Command(["ls", "/"])
    (rc, out, err) = c.wait()
    print("Stdout: " + str(out, encoding="utf-8"))
    print("Stderr: " + str(err, encoding="utf-8"))
    del(c)

    c = RemoteCommand(["ls", "/"], "localhost")
    (rc, out, err) = c.wait()
    print("Stdout: " + str(out, encoding="utf-8"))
    print("Stderr: " + str(err, encoding="utf-8"))
    del(c)

    c = SendFile("WU", "quiche", "/tmp/foo")
    (rc, out, err) = c.wait()
    print("Stdout: " + str(c.stdout, encoding="utf-8"))
    print("Stderr: " + str(c.stderr, encoding="utf-8"))
    del(c)
