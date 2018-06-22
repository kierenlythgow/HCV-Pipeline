#!/usr/bin/env python

import itertools
import os
import psutil
import re
import shlex
import signal
import subprocess as sp
import sys

class snorklib(object):

    def __init__(self):

      # constants

      # err

        self.err = {
        'SuccessReturn'         : [  0, 'Info', 'Successfully completed' ],
        'SuccessHelp'           : [  1, 'Info', 'Successfully printed help message' ],
        'ErrorArgsMissing'      : [  2, 'Erro', 'Failed to retrieve expected command-line argument(s)' ],
        'ErrorArgsInvalid'      : [  3, 'Erro', 'Failed due to invalid command-line argument(s)' ],
        'ErrorArgsParseFailure' : [  4, 'Erro', 'Failed to parse command-line argument(s)' ],
        'ErrorSagaConnect'      : [  5, 'Erro', 'Failed to connect to SCAMPI db' ],
        'ErrorSagaDisconnect'   : [  6, 'Erro', 'Failed to disconnect from SCAMPI db' ],
        'ErrorFileOpen'         : [  7, 'Erro', 'Failed to open file' ],
        'ErrorSysCall'          : [  8, 'Erro', 'Failed system call' ],
        'ErrorFork1Failure'     : [  9, 'Erro', 'Failed to do first process fork' ],
        'ErrorFork2Failure'     : [ 10, 'Erro', 'Failed to do second process fork' ],
        'SuccessUsage'          : [ 11, 'Info', 'Successfully printed usage message' ],
        'SuccessVersion'        : [ 12, 'Info', 'Successfully printed version message' ],
        'ErrorDirCreate'        : [ 13, 'Erro', 'Failed to create directory' ],
        'ErrorExistingLogFile'  : [ 14, 'Erro', 'Failed because log file already exists' ],
        'ErrorDirDoesNotExist'  : [ 15, 'Erro', 'Failed because output directory does not exist' ],
        'SuccessNothingToDo'    : [ 16, 'Info', 'Successfully terminated because nothing to do' ],
        'ErrorIniFileMissing'   : [ 17, 'Erro', 'Failed to find config.ini file' ],
        'ErrorInvalidData'      : [ 18, 'Erro', 'Failed due to invalid input data' ],
        'ErrorExternalSoftware' : [ 19, 'Erro', 'Failed in call to external software' ],
        'SuccessDryRun'         : [ 20, 'Info', 'Successfully terminated in dryrun mode' ],
        'SuccessKilled'         : [ 21, 'Info', 'Successfully terminated by KILL signal' ],
        'ErrorKilled'           : [ 22, 'Erro', 'Failed because terminated by KILL signal' ],
        'ErrorFileMissing'      : [ 23, 'Erro', 'Failed because file is missing' ],
        'ErrorChmod'            : [ 24, 'Erro', 'Failed to change mode of file or directory' ],
        'ErrorScriptCompute'    : [ 25, 'Erro', 'Failed in spawned script' ],
        'ErrorDeletingFile'     : [ 26, 'Erro', 'Failed to delete file' ],
        'ErrorEmptyFile'        : [ 27, 'Erro', 'Failed due to file empty' ],
        'ErrorOutputIncomplete' : [ 28, 'Erro', 'Failed due to incomplete output data' ],
        'ErrorInvalidRetKey'    : [ 29, 'Erro', 'Failed due to invalid return code key in program' ],
        'ErrorSagaDataGet'      : [ 30, 'Erro', 'Failed to retrieve data from SCAMPI' ],
        'ErrorOutfileExists'    : [ 31, 'Erro', 'Failed because output file already exists' ],
        'ErrorDirRename'        : [ 32, 'Erro', 'Failed to rename dir' ],
        'ErrorFileRename'       : [ 33, 'Erro', 'Failed to rename file' ],
        'ErrorReadingData'      : [ 34, 'Erro', 'Failed to read data' ],
        'ErrorFileCopy'         : [ 35, 'Erro', 'Failed to copy file' ],
        'ErrorFileNotLink'      : [ 36, 'Erro', 'Have a file not a symlink' ],
        'ErrorOutputNotComplete': [ 37, 'Erro', 'Output is incomplete' ],
        'ErrorFileCreate'       : [ 38, 'Erro', 'Failed to create file' ]
        }
        self.section = []
        self.option = {}        # self.option[section][var] = val
        self.lasterror = ''

    # cmdfile

        self.cmdfile_logmsgstart = [
            'echo "# ========== #"',
            'echo "SGE job ID: "$JOB_ID',
            'echo "SGE task ID: "$SGE_TASK_ID',
            'echo "SGE job name: "$JOB_NAME',
            'echo "Run on host: "$HOSTNAME',
            'echo "Operating system: "`uname -a`',
            'echo "Username: "$LOGNAME',
            'echo "Started at: "`date +"%Y-%m-%d %H:%M:%S"`',
            'echo "# ========== #',
            'echo "Number of hosts: "$NHOSTS',
            'echo "Number of queues: "$NQUEUES',
            'echo "Number of slots: "$NSLOTS',
            'echo "Job submitted from host: "$SGE_O_HOST',
            'echo "SGE root dir: "$SGE_ROOT',
            'echo "SGE working dir: "$SGE_O_WORKDIR',
            'echo "PATH at job submission: "$SGE_O_PATH',
            'echo "# ========== #',
            'echo'
        ]
        self.cmdfile_logmsgfinish = [
            'echo "Finished at: "`date +"%Y-%m-%d %H:%M:%S"`',
            'echo "# ========== #"'
        ]


    # err

    def err_dump(self):
        'Return a string containing all the error codes.'
        L = [[self.err[key][1], self.err[key][0], key, self.err[key][2]] for key in self.err.keys()]
        sorted(L, key = lambda x: int(x[0]))
        return '\n'.join(['\t'.join([str(x) for x in elt]) for elt in L])

    def err_code(self, retkey):
        'Return the numerical return code value corresponding to the return code descriptor.'
        return self.err[retkey][0]

    def err_retkey(self, retcode):
        'Return the return code descriptor string corresponding to the numerical return code.'
        return [key for key in self.err.keys() if self.err[key][0] == retcode]

    def err_text(self, retkey):
        'Return a string reporting the return code and the standardised message for that return code.'
        if self.err.has_key(retkey):
            return '{type}: Program finished with marcoporo retcode={code} - {msg}\n'.format(
                type=self.err[retkey][1], code=self.err[retkey][0], msg=self.err[retkey][2])
        return ''

    def err_exit(self, retkey, exitprogram=False):
        'Function to report exit code and exit program.'
        if exitprogram:
            sys.exit(self.err[retkey][0])

    def err_last(self):
        'Return the last error message.'
        return self.lasterror

    # str

    def str_2bool(self, s):
        'Return the Python boolean value True or False depending on what the string is.'
        return (s is not None and s.lower() in ['1', 't', 'true', 'y', 'yes'])

    # conf

    def config_read(self, inipath):
        'Read all the data from the inipath into section and option variables; set lasterror if anything goes wrong.'
        line_section = "UNDEFINED"
        optionflat = {}
        with (open(inipath, "r")) as in_fp:
            for line in in_fp:
              # Strip all leading or trailing spaces from the line, replace all tabs by spaces.
              # At this point, there still could be lots of spaces in the 'var' and/or 'val' part.
                line = line.strip().replace("\t", "")
                if (line.startswith(';') or line.startswith('#') or not len(line)):
                    continue
                if (line.startswith('[')):
                    line_section = line.replace("[", "").replace("]", "")
                    self.section.append(line_section)
                    self.option[line_section] = {}
                if ("=" in line):
                    if (line_section == "UNDEFINED"):
                        self.lasterror = "Erro: ini file formatting error: var=val pairs outside of a [section]."
                        in_fp.close()
                        return False
                    L = line.split("=")
                    var = L[0].replace(" ", "")         # Contains no spaces
                    val = "=".join(L[1:]).strip()       # Could contain spaces
                    self.option[line_section][var] = val
                    optionflat[var] = val
            in_fp.close()
      # The section that replaces any instances of {variable} with the value of self.option[section][variable]
      # For every variable, replace {variable} with value in every other value in the dictionary.
        for section, D in self.option.iteritems():
            for var, val in D.iteritems():
                try:
                    localenvvar = re.search("[\w\.]*{(\w*:\w*)}[\w\.]*", val).group(1)
                    localsection, localvar = localenvvar.split(':')
                    oldstring = '{0}{1}{2}'.format('{', localenvvar, '}')
                    newstring = self.option[localsection][localvar]
                    self.option[section][var] = self.option[section][var].replace(oldstring, newstring)
                except:
                    pass
        return True

  # sys

    def sys_pids(self, machineid=None, procname=None, procuser=None, includecpid=False):
        '''
        Return the list of process ids for all processes called procname
        owned by procuser, possibly including the current pid. Does not use
        list comprehensions to find the pid lists because in between getting
        the pid list and testing properties of the pids, the processes may
        terminate and then it is not possible to get the information associated
        with those pids any more.
        '''
        pidlist = []
        pidcurr = psutil.get_pid_list()
        for pid in pidcurr:
            try:
                cmdline = ' '.join(psutil.Process(pid).cmdline())
                if procname is None \
                or (len(psutil.Process(pid).cmdline()) >= 2 and psutil.Process(pid).cmdline()[0].lower().startswith('python') and procname in cmdline):
                    if procuser is None or psutil.Process(pid).username() == procuser:
                        if machineid is None or machineid in cmdline:
                            pidlist.append(pid)
            except:
                pass
        cpid = os.getpid()
        if not includecpid and cpid in pidlist:
            pidlist.remove(cpid)
        return pidlist

    def sys_kill(self, pidlist):
        'Try to kill all process ids in pidlist.'
        faillist = []
        for pid in pidlist:
            try:
                psutil.Process(pid).kill()
            except:
                faillist.append(pid)
                return False
        if len(faillist):
            self.lasterr = 'Failed to kill pids: {0}'.format(','.join([str(x) for x in faillist]))
            return False
        return True

    def sys_exec(self, cmd):
        '''
        Execute a command using the subprocess module to trap the
        return value, the stdout string and stderr string.
        The command may contain pipes.
        '''
        cmdL = [x.strip() for x in cmd.split('|')]
        if not len(cmdL):
            return [0, '', '']
        try:
            proc_handleL = []
            ph = sp.Popen(shlex.split(cmdL[0]), stdout=sp.PIPE, stderr=sp.PIPE)
            proc_handleL.append(ph)
            for cmd in cmdL[1:]:
                ph = sp.Popen(shlex.split(cmd), stdin=proc_handleL[-1].stdout, stdout=sp.PIPE, stderr=sp.PIPE)
                proc_handleL.append(ph)
            proc_stdout, proc_stderr = proc_handleL[-1].communicate()
            proc_returncode = proc_handleL[-1].returncode
        except OSError as e:
            proc_returncode = e.errno
            proc_stdout = 'Error executing: cmd=*{0}*'.format(cmd)
            proc_stderr = e.strerror
        return [proc_returncode, proc_stdout, proc_stderr]
