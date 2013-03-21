import subprocess
import cadologger

class Command(object):
    def __init__(self, args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.logger = Logger()
        # Convert args array to a string for printing if necessary
        if isinstance(self.args, str):
            self.cmdline = self.args
        else:
            self.cmdline = " ".join(self.args)

        self.child = subprocess.Popen(self.args, stdout = subprocess.PIPE, 
            stderr = subprocess.PIPE, **self.kwargs)
        
        self.logger.info("Running command: " + self.cmdline)
        self.logger.cmd(self.cmdline, extra={"pid": self.child.pid})

    def wait(self):
        # Wait for command to finish executing, capturing stdout and stderr 
        # in output tuple
        (self.stdout, self.stderr) = self.child.communicate()

        if self.child.returncode == 0:
            logger.info("Process with PID " + str(self.child.pid) + " finished successfully")
        else:
            logger.error("Process with PID " + str(self.child.pid) + " finished with return code " + str(self.child.returncode))
        self.returncode = self.child.returncode
        return self.returncode

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
    c = cadocommand.Command(["ls", "/"])
    rc = c.wait()
    print("Stdout: " + str(c.stdout, encoding="utf-8"))
    print("Stderr: " + str(c.stderr, encoding="utf-8"))
    del(c)

    c = cadocommand.RemoteCommand(["ls", "/"], "localhost")
    rc = c.wait()
    print("Stdout: " + str(c.stdout, encoding="utf-8"))
    print("Stderr: " + str(c.stderr, encoding="utf-8"))
    del(c)

    c = cadocommand.SendFile("WU", "quiche", "/tmp/foo")
    rc = c.wait()
    print("Stdout: " + str(c.stdout, encoding="utf-8"))
    print("Stderr: " + str(c.stderr, encoding="utf-8"))
    del(c)
