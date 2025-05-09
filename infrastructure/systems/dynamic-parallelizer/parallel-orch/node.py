from enum import Enum, auto
import logging
import os
import re
import executor
from executor import ExecCtxt, ExecResult, ExecArgs
import trace_v2
import util
import signal
from dataclasses import dataclass
from typing import Tuple
from enum import Enum, auto
from pathlib import Path
import util
import analysis

STATE_LOG = '[STATE_LOG] '

def state_log(s):
    # logging.info(STATE_LOG + s)
    pass

class NodeState(Enum):
    INIT = auto()
    READY = auto()
    COMMITTED = auto()
    STOP = auto()
    SPECULATED = auto()
    EXECUTING = auto()
    SPEC_EXECUTING = auto()
    UNSAFE = auto()
    COMMITTED_UNSAFE = auto()

def state_pstr(state: NodeState):
    same_length_state_str = {
        NodeState.INIT:           '  INIT',
        NodeState.READY:          ' READY',
        NodeState.COMMITTED:      'COMMIT',
        NodeState.STOP:           '  STOP',
        NodeState.SPECULATED:     'SPEC_F',
        NodeState.EXECUTING:      '   EXE',
        NodeState.SPEC_EXECUTING: 'SPEC_E',
        NodeState.UNSAFE:         'UNSAFE'
    }
    return same_length_state_str[state]

class RWSet:

    def __init__(self, read_set: set, write_set: set):
        self.read_set = read_set
        self.write_set = write_set

    def add_to_read_set(self, item: str):
        self.read_set.add(item)

    def add_to_write_set(self, item: str):
        self.write_set.add(item)

    def get_read_set(self) -> set:
        return self.read_set

    def get_write_set(self) -> set:
        return self.write_set

    def has_conflict(self, other: 'RWSet') -> bool:
        if (self.write_set.intersection(other.read_set) or
            self.read_set.intersection(other.write_set) or
            self.write_set.intersection(other.write_set)):
            return True
        else:
            return False

    def get_conflict(self, other: 'RWSet') -> set:
        return self.write_set.intersection(other.read_set).union(
            self.read_set.intersection(other.write_set)).union(
                self.write_set.intersection(other.write_set))

    def __str__(self):
        return f"RW(R:{self.get_read_set()}, W:{self.get_write_set()})"


class NodeId:
    def __init__(self, id_: int):
        self.id_ = id_

    def get_non_iter_id(self):
        return NodeId(self.id_)

    def __repr__(self):
        ## TODO: Represent it using n.
        output = f'{self.id_}'
        return output

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        # return self.loop_iters == other.loop_iters and self.id == other.id
        return self.id_ == other.id_

    def __ne__(self, other):
        return not(self == other)

    def __lt__(self, obj):
        return (str(self) < str(obj))

    def __gt__(self, obj):
        return (str(self) > str(obj))

    @staticmethod
    def parse_node_id(node_id_str: str):
        return NodeId(int(node_id_str))


class LoopStack:
    def __init__(self, loop_contexts_or_iters=None):
        if loop_contexts_or_iters is None:
            self.loops = []
        else:
            self.loops = loop_contexts_or_iters

    def __repr__(self):
        ## TODO: Represent it using 'it', 'it0', 'it1', etc
        ##       or -(iters)- in front of it.
        output = "-".join([str(it) for it in self.loops])
        return output
    def __eq__(self, other):
        return self.loops == other.loops


class HSLoopListContext:
    def __init__(self, loop_list_context=None):
        if loop_list_context is None:
            loop_list_context = []
        self.loop_list_context = loop_list_context

    def push(self, loop_list):
        loop_list_context = self.loop_list_context[:]
        loop_list_context.append(loop_list)
        return HSLoopListContext(loop_list_context)

    def get_ith(self, i):
        pass

    def get_top(self):
        return self.loop_list_context[-1][:]

    def pop(self):
        loop_list_context = self.loop_list_context[:]
        loop_list_context.pop()
        return HSLoopListContext(loop_list_context)

def get_loop_list_from_env(env):
    with open(env) as f:
        d = util.parse_env_string_to_dict(f.read())
    if 'IFS' in d:
        ifs = d['IFS']
    else:
        ifs = ' \t\n'
    ifs_ws = ' ' if ' ' in ifs else ''
    ifs_ws += '\t' if '\t' in ifs else ''
    ifs_ws += '\n' if '\n' in ifs else ''
    new_loop_list = re.split(f'[{ifs_ws}]*[{ifs}][{ifs_ws}]*', d['HS_LOOP_LIST'].strip(ifs_ws))
    return new_loop_list

@dataclass
class Node:
    id_: NodeId
    cmd: str
    asts: "list[AstNode]"
    basic_block_id: int
    assignment: bool
    loop_list_change: bool

    def __init__(self, id_, cmd, asts, basic_block_id, var_assignment, loop_list_change):
        self.id_ = id_
        self.cmd = cmd
        self.asts = asts
        self.basic_block_id = basic_block_id
        self.assignment = var_assignment
        self.loop_list_change = loop_list_change

    def is_assignment(self):
        return self.assignment

    def is_loop_list_change(self):
        return self.loop_list_change

    def is_loop_list_push(self):
        return self.loop_list_change and self.cmd.startswith('HS_LOOP_LIST=')

    def is_loop_list_pop(self):
        return self.loop_list_change and self.cmd.startswith('unset')

    def pretty_format(self):
        v = 'q' if self.assignment else ''
        l = 'l' if self.loop_list_change else ''
        return self.cmd.strip() + f'  --- {v}{l} {self.id_}@'

    def simulate_env(self, env):
        return executor.run_assignment_and_return_env_file(self.cmd, env)

    def simulate_loop_list(self, env, loop_list_context: 'HSLoopListContext'):
        assert self.loop_list_change
        if self.cmd == 'unset HS_LOOP_LIST':
            return loop_list_context.pop()
        else:
            new_env = self.simulate_env(env)
            new_loop_list = get_loop_list_from_env(new_env)
            return loop_list_context.push(new_loop_list)

def loop_iters_do_action(loop_iters, edge_type: 'CFGEdgeType'):
    loop_iters_list = list(loop_iters)
    if edge_type == CFGEdgeType.LOOP_BACK:
        loop_iters_list[0] += 1
    elif edge_type == CFGEdgeType.LOOP_SKIP:
        loop_iters_list.pop(0)
    elif edge_type == CFGEdgeType.LOOP_BEGIN:
        loop_iters_list.insert(0, 1)
    elif edge_type == CFGEdgeType.LOOP_END:
        loop_iters_list.pop(0)
    return loop_iters_list

class ConcreteNodeId:
    def __init__(self, node_id: NodeId, loop_iters = list()):
        self.node_id = node_id
        self.loop_iters = tuple(loop_iters)

    def __repr__(self):
        return f'cnid({self.node_id.id_})'

    def __hash__(self):
        return hash((self.node_id, self.loop_iters))

    def __eq__(self, other):
        return self.node_id == other.node_id and self.loop_iters == other.loop_iters

    def __str__(self):
        return f'{self.node_id}@' + ''.join(['-' + str(n) for n in self.loop_iters])

    def do_action(self, next_abstract_id: NodeId, edge_type: 'CFGEdgeType'):
        loop_iters_list = list(self.loop_iters)
        if edge_type == CFGEdgeType.LOOP_BACK:
            loop_iters_list[0] += 1
        elif edge_type == CFGEdgeType.LOOP_SKIP:
            loop_iters_list.pop(0)
        elif edge_type == CFGEdgeType.LOOP_BEGIN:
            loop_iters_list.insert(0, 1)
        elif edge_type == CFGEdgeType.LOOP_END:
            loop_iters_list.pop(0)
        return ConcreteNodeId(next_abstract_id, loop_iters_list)

    @staticmethod
    def parse(input_str):
        node_id_str, loop_iters_str = input_str.split('@')
        return ConcreteNodeId(NodeId(int(node_id_str)), [int(cnt) for cnt in loop_iters_str.split('-')[1:]])

class ConcreteNode:
    cnid: ConcreteNodeId
    abstract_node: Node
    state: NodeState
    # exists for EXEC or SPEC_E or subsequent states, erased for READY
    exec_id: int
    # Nodes to check for fs dependencies before this node can be committed
    # for this particular execution of the main sandbox.
    # No need to do the same for the background sandbox since it will always get committed.
    to_be_resolved_snapshot: "set[NodeId]"
    # Read and write sets for this node
    rwset: RWSet
    # The wait trace file for this node
    wait_env_file: str
    # This can only be set while in the frontier and the background node execution is enabled
    # TODO: For now ignore this. Maybe there is a better way to do this.
    # background_sandbox: Sandbox

    # Exists when the node is in EXE or SPEC_EXE or after those states
    exec_ctxt: ExecCtxt

    # Exists when the node is in COMMITED or SPEC_F
    exec_result: ExecResult

    # Updated when the node is loop changing and the node is transitioning
    # into COMMITTED or SPEC_F
    loop_list_context: HSLoopListContext

    # read-only, the value of initial loop_list_context
    # used when reset_to_ready
    init_loop_list_context: HSLoopListContext
    
    spec_pre_env: str

    # Exists when node is in READY
    assignments: "list[NodeId]"

    # Exists when node is in EXE or SPEC_EXE, it acts as a cache for
    # the trace file content
    trace_lines: list
    # Exists when node is in EXE or SPEC_EXE, it it an opened file
    # or none when such file doesn't exist
    trace_fd=None
    trace_ctx=None

    def __init__(self, cnid: ConcreteNodeId, node: Node, loop_list_context: HSLoopListContext,
                 spec_pre_env=None):
        self.cnid = cnid
        self.abstract_node = node
        self.state = NodeState.INIT
        self.tracefile = None
        self.rwset = None
        self.wait_env_file = None
        self.exec_ctxt = None
        self.exec_id = None
        self.spec_pre_env = spec_pre_env
        self.loop_list_context = loop_list_context
        self.init_loop_list_context = loop_list_context
        self.trace_fd = None
        self.init_trace_lines()

    def __str__(self):
        return f'Node(id:{self.id_}, cmd:{self.cmd}, state:{self.state}, wait_env_file:{self.wait_env_file}, exec_ctxt:{self.exec_ctxt})'

    def __repr__(self):
        return str(self)

    @property
    def id_(self):
        return self.abstract_node.id_

    @property
    def cmd(self):
        return self.abstract_node.cmd

    @property
    def asts(self):
        return self.abstract_node.asts

    def pretty_state_repr(self):
        return f'{state_pstr(self.state)} {self.cmd} --- {self.cnid}'

    def is_initialized(self):
        return self.state == NodeState.INIT

    def is_ready(self):
        return self.state == NodeState.READY

    def is_committed(self):
        return self.state == NodeState.COMMITTED

    def is_stopped(self):
        return self.state == NodeState.STOP

    def is_speculated(self):
        return self.state == NodeState.SPECULATED

    def is_executing(self):
        return self.state == NodeState.EXECUTING

    def is_spec_executing(self):
        return self.state == NodeState.SPEC_EXECUTING

    def is_unsafe(self):
        return self.state == NodeState.UNSAFE


    # This function fixes up fds redirected by the runtime
    # which facilitates fd edit merging optimization.
    def fixup_fds(self):
        replace_map = {}
        with open(self.exec_ctxt.pre_env_file + '.fds', 'r') as f:
            "line format: fd mode offset path"
            lines = f.read().split('\n')[:-1]
            for line in lines:
                fd, mode, offset, path = line.split(' ', maxsplit=3)
                replace_map[f'{self.exec_ctxt.outfds}/{fd}'] = path
        new_lines = []
        if self.state == NodeState.SPEC_EXECUTING:
            post_path = util.sandboxed_path(self.exec_ctxt.sandbox_dir,
                                            self.exec_ctxt.post_env_file + '.fds')
        else:
            post_path = self.exec_ctxt.post_env_file + '.fds'
        with open(post_path, 'r') as f:
            "line format: fd mode offset path"
            lines = f.read().split('\n')[:-1]
            for line in lines:
                fd, mode, offset, path = line.split(' ', maxsplit=3)
                offset = int(offset)
                if path in replace_map:
                    path = replace_map[path]
                    if not path.startswith('pipe:['):
                        offset += os.path.getsize(path)
                new_lines.append((fd, mode, str(offset), path))
        with open(post_path, 'w') as f:
            for line in new_lines:
                f.write(' '.join(line))
                f.write('\n')

    def start_command(self, env_file: str, speculate=False, speculated_nodes=None):
        # TODO: implement speculate
        # TODO: built-in commands
        execute_func = executor.run_trace_sandboxed
        if speculated_nodes is None:
            lower_sandboxes = []
        else:
            lower_sandboxes = [node.exec_ctxt.sandbox_dir for node in reversed(speculated_nodes)]
        # Set the execution id
        self.exec_id = util.generate_id()
        self.exec_ctxt = execute_func(ExecArgs(command=self.cmd, concrete_node_id=self.cnid, execution_id=self.exec_id, pre_execution_env_file=env_file, speculate_mode=speculate, lower_sandboxes=lower_sandboxes))
        util.debug_log(f'Node {self.cnid} executing with pid {self.exec_ctxt.process.pid}')

    def execution_outcome(self) -> Tuple[int, str, str]:
        assert self.exec_result is not None
        return self.exec_result.exit_code, self.exec_ctxt.post_env_file, self.exec_ctxt.outfds

    def command_unsafe(self):
        if len(self.asts) == 0:
            return True
        return not analysis.safe_to_execute(self.asts, {})

    def update_loop_list_context(self):
        if self.abstract_node.is_loop_list_push():
            real_env_path = util.sandboxed_path(self.exec_ctxt.sandbox_dir,
                                                self.exec_ctxt.post_env_file)
            new_loop_list = get_loop_list_from_env(real_env_path)
            self.loop_list_context = self.loop_list_context.push(new_loop_list)
        elif self.abstract_node.is_loop_list_pop():
            self.loop_list_context = self.loop_list_context.pop()

    def guess_post_env(self):
        if self.command_unsafe() or self.is_unsafe():
            env_file = self.spec_pre_env
        elif self.is_committed():
            env_file = self.exec_ctxt.post_env_file
        elif self.is_speculated():
            env_file = util.sandboxed_path(self.exec_ctxt.sandbox_dir,
                                               self.exec_ctxt.post_env_file)
        elif self.is_executing():
            env_file = self.exec_ctxt.pre_env_file
        else:
            env_file = self.spec_pre_env
        assert env_file is not None
        return env_file

    def trace_state(self):
        state_log(f'{self.cnid}: {state_pstr(self.state)}')
        if self.state in [NodeState.EXECUTING, NodeState.SPEC_EXECUTING, NodeState.COMMITTED, NodeState.SPECULATED]:
            state_log("id: {}, pre_env: {}, sandbox: {}, out_fd_dir: {}, trace: {}".format(
                self.cnid,
                self.exec_ctxt.pre_env_file,
                self.exec_ctxt.sandbox_dir,
                self.exec_ctxt.outfds,
                self.exec_ctxt.trace_file
            ))

    def kill(self):
        assert self.state in [NodeState.EXECUTING, NodeState.SPEC_EXECUTING]
        self.exec_ctxt.process.kill()

    def init_trace_lines(self):
        self.trace_lines = ['']

    def gather_fs_actions(self) -> RWSet:
        assert self.state in [NodeState.EXECUTING, NodeState.SPEC_EXECUTING]
        sandbox_dir = self.exec_ctxt.sandbox_dir
        trace_file = self.exec_ctxt.trace_file
        if self.trace_fd is None:
            try:
                self.trace_fd = open(util.sandboxed_path(sandbox_dir, trace_file))
                self.trace_ctx = trace_v2.Context()
                self.trace_ctx.set_dir(os.getcwd())
            except FileNotFoundError:
                return
        new_trace = self.trace_fd.read()
        new_lines = new_trace.split('\n')
        start_parse = len(self.trace_lines)-1
        self.trace_lines[-1] = self.trace_lines[-1] + new_lines[0]
        self.trace_lines.extend(new_lines[1:])
        stop_parse = len(self.trace_lines)-1
        read_set, write_set = trace_v2.parse_and_gather_cmd_rw_sets(
            self.trace_lines[start_parse:stop_parse], self.trace_ctx)
        self.update_rw_set(read_set, write_set)

    def get_rw_set(self):
        # if self.state in [NodeState.EXECUTING, NodeState.SPEC_EXECUTING]:
        #     self.gather_fs_actions()
        return self.rwset

    fd_line = re.compile(r'(\d+) ([rwd]) (\d+) (.+)')
    def has_env_conflict_with(self, other_env) -> bool:
        # Early return if paths are the same
        if self.exec_ctxt.pre_env_file == other_env:
            return False

        ignore_vars = set([
            "_", 'RANDOM', "msg", "pash_runtime_final_status", "pash_previous_set_status",
            "pash_runtime_shell_variables_file", "from_set", "output_variable_file",
            "pash_loop_iter_counters", "daemon_response", "vars_file",
            "pash_speculative_command_id", "prev_env", "PREVIOUS_SET_STATUS",
            "BASH_LINENO", "response_args", "stdout_file", "pash_spec_command_id",
            "cmd_exit_code", "pash_set_to_add", "POST_EXEC_ENV",
            "HISTCMD", "EXEC_MODE", "TRACE_FILE", "CMD_STRING",
            "CMD_ID", "STDOUT_FILE", "DIRSTACK", "SECONDS", "TMPDIR",
            "UPDATED_DIRS_AND_MOUNTS", "EPOCHSECONDS", "LATEST_ENV_FILE",
            "TRY_COMMAND", "SRANDOM", "speculate_flag", "EXECUTION_ID",
            "EPOCHREALTIME", "OLDPWD", "exit_code", "BASHPID", "BASH_COMMAND", "BASH_ARGV0",
            "cmd", "BASH_ARGC", "BASH_ARGV", "BASH_SUBSHELL", "LINENO", "GROUPS", "BASH_SOURCE",
            "PREVIOUS_SHELL_EC", "pash_previous_exit_status", "filter_vars_file", "pash_spec_loop_id",
            "pash_loop_iters", "LINES", "COLUMNS",
        ])

        ignore_prefix = "pash_loop_"

        re_scalar_string = re.compile(r'declare (?:-x|--)? (\w+)="([^"]*)"')
        re_scalar_int = re.compile(r'declare -i (\w+)="(\d+)"')
        re_array = re.compile(r'declare -a (\w+)=(\([^)]+\))')
        re_fn = re.compile(r'declare -fx (\w+)=(\([^)]+\))')

        def parse_env(content):
            env_vars = {}
            for line in content.splitlines():
                if line.startswith('#') or not line.strip():
                    continue
                for regex in [re_scalar_string, re_scalar_int, re_array]:
                    match = regex.match(line)
                    if match:
                        key, value = match.groups()
                        if key not in ignore_vars and not key.startswith(ignore_prefix):
                            env_vars[key] = value
                        break
            inside_function = False
            current_function = ''
            function_body_lines = []
            for line in content.splitlines():
                if line.startswith('#') or not line.strip():
                    continue
                if not inside_function and not line.startswith('declare') and line.endswith('() '):
                    inside_function = True
                    current_function = line[:-len(' () ')]
                elif inside_function:
                    function_body_lines.append(line)
                    if line == '}':
                        inside_function = False
                        if not current_function in ignore_vars:
                            env_vars[current_function] = '\n'.join(function_body_lines)
                        function_body_lines = []
            return env_vars

        with open(self.exec_ctxt.pre_env_file, 'r') as file:
            node_env_vars = parse_env(file.read())

        with open(other_env, 'r') as file:
            other_env_vars = parse_env(file.read())

        util.env_log(f"Comparing env files {self.exec_ctxt.pre_env_file} and {other_env}")

        conflict_exists = False
        for key in set(node_env_vars.keys()).union(other_env_vars.keys()):
            if key not in node_env_vars:
                util.env_log(f"Variable {key} missing in node environment")
                conflict_exists = True
            elif key not in other_env_vars:
                util.env_log(f"Variable {key} missing in other environment")
                conflict_exists = True
            elif node_env_vars[key] != other_env_vars[key]:
                util.env_log(f"Variable {key} differs: node environment has {node_env_vars[key]}, other has {other_env_vars[key]}")
                conflict_exists = True

        with open(self.exec_ctxt.pre_env_file + '.fds', 'r') as file1, open(other_env + '.fds', 'r') as file2:
            # Since we assume stdin, stdout, and stderr don't change during the script
            # We are omitting them from the comparison
            s1 = file1.read().strip().split('\n')
            s1 = [ConcreteNode.fd_line.match(line).groups() for line in s1]
            s2 = file2.read().strip().split('\n')
            s2 = [ConcreteNode.fd_line.match(line).groups() for line in s2]
            if len(s1) != len(s2):
                util.env_log(f"fds diff: \n{s1}\n{s2}")
                conflict_exists = True
            else:
                for a, b in zip(s1, s2):
                    if (a[0] != b[0] or a[1] != b[1] or (a[1] == 'r' and a[2] != b[2])
                        or a[3] != b[3]):
                        util.env_log(f"fds diff: \n{s1}\n{s2}")
                        conflict_exists = True
                        break

        return conflict_exists

    def kill_children(self):
        util.overhead_log(f"KILL|{self.cnid}")
        assert self.state in [NodeState.EXECUTING, NodeState.SPEC_EXECUTING]
        process = self.exec_ctxt.process
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()
        util.overhead_log(f"KILL_END|{self.cnid}")

    def commit_fd_writes(self):
        with open(self.exec_ctxt.pre_env_file + '.fds', 'r') as f:
            "line format: fd mode offset path"
            lines = f.read().split('\n')[:-1]
            for line in lines:
                fd, mode, offset, path = line.split(' ', maxsplit=3)
                if mode == 'w':
                    util.append(self.exec_ctxt.outfds + '/' + fd, path)

    def runtime_finished(self):
        # TODO: update this when exec doesn't use sandbox anymore
        if self.state in [NodeState.SPEC_EXECUTING, NodeState.EXECUTING]:
            post_path = util.sandboxed_path(self.exec_ctxt.sandbox_dir,
                                            self.exec_ctxt.post_env_file + '.fds')
        # elif self.state == NodeState.EXECUTING:
        #     post_path = self.exec_ctxt.post_env_file + '.fds'
        else:
            assert False
        return Path(post_path).exists()

    ##                                      ##
    ##          Transition Functions        ##
    ##                                      ##

    def transition_from_init_to_ready(self, spec_pre_env):
        assert self.state == NodeState.INIT
        self.state = NodeState.READY
        self.rwset = RWSet(set(), set())
        self.spec_pre_env = spec_pre_env
        # self.spec_pre_env = ConcreteAssignmentNode.execute_assignments_and_get_most_recent_spec_pre_env(assignments)
        # Also, probably unroll here?
        self.trace_state()

    def transition_from_ready_to_unsafe(self):
        assert self.state == NodeState.READY
        self.state = NodeState.UNSAFE
        self.trace_state()

    def try_reset_to_ready(self, spec_pre_env: str=None, loop_list_context=None):
        if self.state in [NodeState.READY, NodeState.UNSAFE]:
            return
        else:
            self.reset_to_ready(spec_pre_env, loop_list_context)

    def reset_to_ready(self, spec_pre_env: str = None, loop_list_context: HSLoopListContext = None):
        assert self.state in [NodeState.EXECUTING, NodeState.SPEC_EXECUTING,
                              NodeState.SPECULATED]

        state_log(f"Resetting node {self.id_} to ready {self.exec_id}")
        # We reset the exec id so if we receive a message
        # due to a race condition, we will ignore it.
        self.exec_id = None

        # TODO: make this more sophisticated
        if self.state in [NodeState.EXECUTING, NodeState.SPEC_EXECUTING]:
            self.kill_children()
        util.overhead_log(f"DELETE_SANDBOX|{self.cnid}")
        util.delete_sandbox(self.exec_ctxt.sandbox_dir)
        util.overhead_log(f"DELETE_SANDBOX_END|{self.cnid}")
        self.exec_ctxt = None
        self.exec_result = None
        if loop_list_context is not None:
            self.init_loop_list_context = loop_list_context
        self.loop_list_context = self.init_loop_list_context
        if spec_pre_env is not None:
            self.spec_pre_env = spec_pre_env
        self.init_trace_lines()
        if self.trace_fd is not None:
            self.trace_fd.close()
            self.trace_fd = None
            self.trace_ctx = None
        self.state = NodeState.READY
        self.trace_state()

    def start_executing(self, env_file):
        assert self.state == NodeState.READY

        self.start_command(env_file)
        self.state = NodeState.EXECUTING
        self.init_trace_lines()
        self.rwset = RWSet(set(), set())
        self.trace_state()

    def start_spec_executing(self, env_file, speculated_nodes):
        # raise NotImplementedError
        assert self.state == NodeState.READY
        self.start_command(env_file, speculate=True, speculated_nodes=speculated_nodes)
        self.init_trace_lines()
        self.rwset = RWSet(set(), set())
        self.state = NodeState.SPEC_EXECUTING
        self.trace_state()

    def collect_result(self):
        assert self.state in [NodeState.EXECUTING, NodeState.SPEC_EXECUTING]
        self.exec_ctxt.process.wait()
        self.exec_result = ExecResult(self.exec_ctxt.process.returncode, self.exec_ctxt.process.pid)
        return self.exec_result.exit_code == 137, self.runtime_finished()

    def commit_frontier_execution(self):
        assert self.state == NodeState.EXECUTING
        self.gather_fs_actions()
        self.init_trace_lines()
        self.kill_children()
        if self.trace_fd is not None:
            self.trace_fd.close()
            self.trace_fd = None
            self.trace_ctx = None
        self.update_loop_list_context()
        util.overhead_log(f"COMMIT|{self.cnid}")
        executor.commit_workspace(self.exec_ctxt.sandbox_dir)
        # self.commit_fd_writes()
        util.overhead_log(f"COMMIT_END|{self.cnid}")
        # util.delete_sandbox(self.exec_ctxt.sandbox_dir)
        self.fixup_fds()
        self.state = NodeState.COMMITTED
        self.trace_state()

    def finish_spec_execution(self):
        assert self.state == NodeState.SPEC_EXECUTING
        self.update_loop_list_context()
        self.gather_fs_actions()
        self.init_trace_lines()
        self.kill_children()
        if self.trace_fd is not None:
            self.trace_fd.close()
            self.trace_fd = None
            self.trace_ctx = None
        self.fixup_fds()
        self.state = NodeState.SPECULATED
        self.trace_state()

    def commit_speculated(self):
        assert self.state == NodeState.SPECULATED
        util.overhead_log(f"COMMIT|{self.cnid}")
        executor.commit_workspace(self.exec_ctxt.sandbox_dir)
        # self.commit_fd_writes()
        util.overhead_log(f"COMMIT_END|{self.cnid}")
        # util.delete_sandbox(self.exec_ctxt.sandbox_dir)
        self.state = NodeState.COMMITTED
        self.trace_state()

    def transition_from_stopped_to_executing(self, env_file=None):
        assert self.state == NodeState.READY
        self.state = NodeState.EXECUTING
        self._attempt_start_command(env_file)

    def transition_from_spec_executing_to_speculated(self):
        pass

    def commit_unsafe_node(self):
        assert self.state == NodeState.UNSAFE
        self.state = NodeState.COMMITTED_UNSAFE

    def update_rw_set(self, r_set, w_set):
        for rfile in r_set:
            self.rwset.add_to_read_set(rfile)
        for wfile in w_set:
            self.rwset.add_to_write_set(wfile)



class CFGEdgeType(Enum):
    IF_TAKEN = auto()
    ELSE_TAKEN = auto()
    LOOP_TAKEN = auto()
    LOOP_SKIP = auto()
    LOOP_BACK = auto()
    LOOP_BEGIN = auto()
    LOOP_END = auto()
    OTHER = auto()

class HSBasicBlock:
    def __init__(self, bb_id: int, nodes: list[Node]):
        self.bb_id = bb_id
        self.nodes = nodes

    def __str__(self):
        return ''.join([node.pretty_format() + '\n' for node in self.nodes])

    @property
    def loop_context(self):
        return self.nodes[0].loop_context

    @property
    def node_ids(self):
        return [node.id_ for node in self.nodes]

    def get_node(self, node_id: NodeId) -> Node:
        nodes = [node for node in self.nodes if node.id_ == node_id]
        assert len(nodes) == 1
        return nodes[0]

class HSProg:
    basic_blocks: list[HSBasicBlock] = []
    block_adjacency: "dict[int, dict[int, CFGEdgeType]]"

    def __init__(self, basic_blocks: list, block_edges: list[tuple]):
        self.basic_blocks = [HSBasicBlock(i, []) for i, bb in enumerate(basic_blocks)]
        self.block_adjacency = {}
        for bb_id in range(len(basic_blocks)):
            self.block_adjacency[bb_id] = {}

        for from_bb, to_bb, edge_type, aux_info in block_edges:
            self.block_adjacency[from_bb][to_bb] = (CFGEdgeType[edge_type], aux_info)

    def is_start_of_block(self, node_id: NodeId):
        for bb in self.basic_blocks:
            bb : HSBasicBlock
            if len(bb.nodes) and bb.nodes[0].id_ == node_id:
                return True
        return False

    def append_node_to(self, bb_id, node: Node):
        self.basic_blocks[bb_id].nodes.append(node)

    def find_basic_block(self, node_id: NodeId):
        for bb in self.basic_blocks:
            bb : HSBasicBlock
            for node in bb.nodes:
                if node.id_ == node_id:
                    return bb
        raise ValueError('no such node_id')

    def is_last_block(self, bb: HSBasicBlock):
        bb_id = self.basic_blocks.index(bb)
        if len(self.block_adjacency[bb_id]) == 0:
            return True
        else:
            return False

    def guess_next_block(self, bb: HSBasicBlock, loop_iters: list,
                         loop_list_context: HSLoopListContext):
        bb_id = self.basic_blocks.index(bb)
        pick_dict = {}
        for next_bb_id, (edge_type, aux_info) in self.block_adjacency[bb_id].items():
            pick_dict[edge_type] = (next_bb_id, aux_info)
        if CFGEdgeType.LOOP_BEGIN in pick_dict:
            assert len(pick_dict) == 1
            return (CFGEdgeType.LOOP_BEGIN, self.basic_blocks[pick_dict[edge_type][0]],
                    pick_dict[edge_type][1])
        elif CFGEdgeType.LOOP_SKIP in pick_dict:
            assert CFGEdgeType.LOOP_TAKEN in pick_dict
            if len(loop_list_context.get_top()) < loop_iters[0]:
                return (CFGEdgeType.LOOP_SKIP,
                        self.basic_blocks[pick_dict[CFGEdgeType.LOOP_SKIP][0]],
                        pick_dict[edge_type][1])
            else:
                return (CFGEdgeType.LOOP_TAKEN,
                        self.basic_blocks[pick_dict[CFGEdgeType.LOOP_TAKEN][0]],
                        pick_dict[edge_type][1])
        for edge_type in [CFGEdgeType.LOOP_END, CFGEdgeType.LOOP_BACK,
                          CFGEdgeType.IF_TAKEN, CFGEdgeType.ELSE_TAKEN,
                          CFGEdgeType.OTHER]:
            if edge_type in pick_dict:
                return edge_type, self.basic_blocks[pick_dict[edge_type][0]], pick_dict[edge_type][1]
        assert False

    def find_node(self, node_id):
        for bb in self.basic_blocks:
            for node in bb.nodes:
                if node.id_ == node_id:
                    return node
        raise ValueError('no such node_id')

    def __str__(self):
        return 'prog:\n' + '\n'.join(
            [f'block {i}:\n' + str(bb) + f'goto block {self.block_adjacency[i]}\n' for i, bb in enumerate(self.basic_blocks)])
