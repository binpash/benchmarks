import os
import subprocess
import logging
import time


## TODO: Figure out how logging here plays out together with the log() in PaSh

## Setup TRACE level for logging
level = logging.TRACE = logging.DEBUG - 5 

def log_logger(self, message, *args, **kwargs):
    if self.isEnabledFor(level):
        self._log(level, message, args, **kwargs)
logging.getLoggerClass().trace = log_logger

def log_root(msg, *args, **kwargs):
    logging.log(level, msg, *args, **kwargs)
logging.addLevelName(level, 'TRACE')
logging.trace = log_root


GIT_TOP_CMD = [ 'git', 'rev-parse', '--show-toplevel', '--show-superproject-working-tree']
if 'PASH_SPEC_TOP' in os.environ:
    PASH_SPEC_TOP = os.environ['PASH_SPEC_TOP']
else:
    PASH_SPEC_TOP = subprocess.run(GIT_TOP_CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).stdout.rstrip()


## Ensure that PASH_TMP_PREFIX is set by pa.sh
PASH_SPEC_TMP_PREFIX = os.getenv('PASH_SPEC_TMP_PREFIX')

SOCKET_BUF_SIZE = 8192

SCHEDULER_SOCKET = os.getenv('PASH_SPEC_SCHEDULER_SOCKET')

INSIGNIFICANT_VARS = {'PWD', 'OLDPWD', 'SHLVL', 'PASH_SPEC_TMP_PREFIX', 'PASH_SPEC_SCHEDULER_SOCKET', 'PASH_SPEC_TOP',
                      'PASH_TOP', 'PASH_TOP_LEVEL','RANDOM', 'LOGNAME', 'MACHTYPE', 'MOTD_SHOWN', 'OPTERR', 'OPTIND',
                      'PPID', 'PROMPT_COMMAND', 'PS4', 'SHELL', 'SHELLOPTS', 'SHLVL', 'TERM', 'UID', 'USER', 'XDG_SESSION_ID'}

SIGNIFICANT_VARS = {'foo', 'bar', 'baz', 'file1', 'file2', 'file3', 'file4', 'file5', 'LC_ALL', 'nchars', 'filename'}

START_TIME = time.time()

NAMED_TIMESTAMPS = {}

SANDBOX_KILLING = False
SPECULATE_IMMEDIATELY = False