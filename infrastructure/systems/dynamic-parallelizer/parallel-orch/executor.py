import config
import logging
import subprocess
import util
import os

from dataclasses import dataclass

@dataclass
class ExecCtxt:
    process: subprocess.Popen
    trace_file: str
    outfds: str
    stderr: str
    pre_env_file: str
    post_env_file: str
    sandbox_dir: str

@dataclass
class ExecResult:
    exit_code: int
    proc_id: int

@dataclass
class ExecArgs:
    command: str
    concrete_node_id: "ConcreteNodeId"
    execution_id: int
    pre_execution_env_file: str
    speculate_mode: bool
    lower_sandboxes: list[str]

# This module executes a sequence of commands
# and traces them with Riker.
# All commands are run inside an overlay sandbox.

def set_pgid():
    os.setpgid(0, 0)

def run_assignment_and_return_env_file(assignment: str, pre_execution_env_file: str):
    post_execution_env_file = util.ptempfile(prefix='hs_assignment_post_env')
    logging.debug(f'Running assignment: {assignment} | pre_execution_env_file: {pre_execution_env_file} | post_execution_env_file: {post_execution_env_file}')
    run_script = f'{config.PASH_SPEC_TOP}/parallel-orch/run_assignment.sh'
    args = ["/bin/bash", run_script, assignment, pre_execution_env_file, post_execution_env_file]
    process = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    util.copy(pre_execution_env_file + '.fds', post_execution_env_file + '.fds')
    return post_execution_env_file

def run_trace_sandboxed(args: ExecArgs):
    run_script = f'{config.PASH_SPEC_TOP}/parallel-orch/run_command.sh'

    trace_file  = util.ptempfile(prefix='hs_trace')
    outfiles_dir = util.ptempdir(prefix='hs_outfiles')
    stderr_file = util.ptempfile(prefix='hs_stderr')
    logging.debug(f'Scheduler: Trace file for: {args.concrete_node_id}: {trace_file}')
    logging.debug(f'Scheduler: Stdout file for: {args.concrete_node_id} is: {outfiles_dir}')
    logging.debug(f'Scheduler: Stderr file for: {args.concrete_node_id} is: {stderr_file}')

    sandbox_dir, tmp_dir = util.create_sandbox()
    post_execution_env_file = util.ptempfile(prefix='hs_post_env')
    lower_dirs_str = ':'.join(args.lower_sandboxes)
    speculate_mode = "speculate" if args.speculate_mode else "standard"

    cmd = ["/bin/bash", run_script, args.command, trace_file, outfiles_dir, args.pre_execution_env_file, sandbox_dir, tmp_dir, speculate_mode, str(args.concrete_node_id), post_execution_env_file, str(args.execution_id), lower_dirs_str ]
    logging.debug(cmd)
    process = subprocess.Popen(cmd, stdout=None, stderr=None, preexec_fn=set_pgid)
    # For debugging
    # process = subprocess.Popen(cmd)

    return ExecCtxt(process, trace_file, outfiles_dir, stderr_file, args.pre_execution_env_file, post_execution_env_file, sandbox_dir)

def commit_workspace(workspace_path):
    ## Call commit-sandbox.sh to commit the uncommitted sandbox to the main workspace
    run_script = f'{config.PASH_SPEC_TOP}/deps/try/try'
    args = ["/bin/bash", run_script, "-i", "/run/mount", "commit", workspace_path]
    process = subprocess.check_output(args)
    return process

## Read trace and capture each command
def read_trace(sandbox_dir, trace_file):
    if sandbox_dir == "":
        path = trace_file
    else:
        path = f"{sandbox_dir}/upperdir/{trace_file}"
    logging.debug(f'Reading trace from: {path}')
    with open(path) as f:
        return f.read().split('\n')[:-1]

def read_env_file(env_file, sandbox_dir=None):
    if sandbox_dir is None:
        path = env_file
    else:
        path = f"{sandbox_dir}/upperdir/{env_file}"
    logging.debug(f'Reading env from: {path}')
    out = subprocess.check_output([f"{os.getenv('PASH_TOP')}/compiler/orchestrator_runtime/pash_filter_vars.sh", path])
    return out.decode("utf-8")
