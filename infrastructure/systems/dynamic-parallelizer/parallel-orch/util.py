import config
import logging
import os
import socket
import subprocess
import tempfile
import time
import re
import psutil
import signal
import analysis
import shutil
from node import Node, NodeId, LoopStack, HSProg, HSBasicBlock
from partial_program_order import PartialProgramOrder
from config import PASH_SPEC_TMP_PREFIX

DEBUG_LOG = '[DEBUG_LOG] '
ENV_LOG = '[ENV_LOG] '
GOOD_LOG = '[GOOD_LOG] '

def debug_log(s):
    logging.info(DEBUG_LOG + s)

def env_log(s):
    logging.info(ENV_LOG + s)

def overhead_log(s):
    # logging.debug(DEBUG_LOG + s)
    pass

def good_log(s):
    logging.info(GOOD_LOG + s)

def ptempfile(prefix=''):
    fd, name = tempfile.mkstemp(dir=config.PASH_SPEC_TMP_PREFIX, prefix=prefix+'_')
    ## TODO: Get a name without opening the fd too if possible
    os.close(fd)
    return name

def ptempdir(prefix=''):
    name = tempfile.mkdtemp(dir=config.PASH_SPEC_TMP_PREFIX, prefix=prefix+'_')
    return name

def copy(path_from, path_to):
    shutil.copy(path_from, path_to)

def append(path_from, path_to):
    with open(path_from, 'r') as source, open(path_to, 'a') as destination:
        shutil.copyfileobj(source, destination)

def cp_to_ptmpfile(source, prefix=''):
    fd, name = tempfile.mkstemp(dir=config.PASH_SPEC_TMP_PREFIX, prefix=prefix+'_')
    os.close(fd)
    shutil.copy(source, name)
    return name

def create_sandbox():
    os.makedirs(f"{config.PASH_SPEC_TMP_PREFIX}/tmp/pash_spec/a", exist_ok=True)
    os.makedirs(f"{config.PASH_SPEC_TMP_PREFIX}/tmp/pash_spec/b", exist_ok=True)
    sdir = tempfile.mkdtemp(dir=f"{config.PASH_SPEC_TMP_PREFIX}/tmp/pash_spec/a", prefix="sandbox_")
    tdir = tempfile.mkdtemp(dir=f"{config.PASH_SPEC_TMP_PREFIX}/tmp/pash_spec/b", prefix="sandbox_")
    return sdir, tdir

def delete_sandbox(sandbox):
    if not sandbox.startswith(f"{config.PASH_SPEC_TMP_PREFIX}/tmp/pash_spec/a"):
        breakpoint()
    assert sandbox.startswith(f"{config.PASH_SPEC_TMP_PREFIX}/tmp/pash_spec/a")
    shutil.rmtree(os.path.join(sandbox, 'upperdir'), ignore_errors=True)

def sandboxed_path(sandbox_dir, path):
    if len(sandbox_dir):
        return f"{sandbox_dir}/upperdir/{path}"
    else:
        return path

def init_unix_socket(socket_file: str) -> socket.socket:
    server_address = socket_file

    # Make sure the socket does not already exist
    ## TODO: Is this necessary?
    try:
        os.unlink(server_address)
    except OSError:
        if os.path.exists(server_address):
            raise
    logging.debug("SocketManager: Made sure that socket does not exist")

    # Create a UDS socket
    sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    logging.debug("SocketManager: Created socket")

    sock.bind(server_address)
    logging.debug("SocketManager: Successfully bound to socket")

    ## TODO: Check if we need to configure the backlog
    sock.listen()
    logging.debug("SocketManager: Listenting on socket")

    return sock

def socket_get_next_cmd(sock: socket.socket) -> "tuple[socket.socket, str]" :
    connection, client_address = sock.accept()
    data = connection.recv(config.SOCKET_BUF_SIZE)

    ## TODO: This could be avoided for efficiency
    str_data = data.decode('utf-8')
    logging.debug(f'Received data: {str_data}')
    ## TODO: Lift this requirement if needed
    ##
    ## We need to ensure that we read a command at once or the command was empty (only relevant in the first invocation)
    assert(str_data.endswith("\n") or str_data == "")

    return (connection, str_data)

def socket_respond(connection: socket.socket, message: str):
    bytes_message = message.encode('utf-8')
    connection.sendall(bytes_message)
    connection.close()

def parse_env_string_to_dict(content):
    # Parse scalar string vars
    scalar_vars_string = re.findall(r'declare (?:-x|--)? (\w+)="([^"]*)"', content, re.DOTALL)

    # regex magic, capturing string with quotation
    escape_vars_string = re.findall(r"declare (?:-x|--)? (\w+)=\$'((?:\\.|[^'])*)'", content, re.DOTALL)

    # Parse scalar integer vars
    scalar_vars_int = re.findall(r'declare -i (\w+)="(\d+)"', content)

    # Parse array declarations
    array_vars = re.findall(r'declare -a (\w+)=(\([^)]+\))', content)

    # Merge all parsed variables
    result = {key: value for key, value in scalar_vars_string}
    result.update({key: int(value) for key, value in scalar_vars_int})
    result.update({key: value for key, value in array_vars})
    result.update({key: value.encode('ascii').decode('unicode_escape')
                   for key, value in escape_vars_string})

    return result

def compare_dicts(dict1, dict2):
    only_in_first = {}
    only_in_second = {}
    different_in_both = {}
    # Check for keys in dict1 but not in dict2 and for different values
    for key, value in dict1.items():
        if key not in dict2:
            only_in_first[key] = value
        elif dict1[key] != dict2[key]:
            different_in_both[key] = (dict1[key], dict2[key])
    # Check for keys in dict2 but not in dict1
    for key, value in dict2.items():
        if key not in dict1:
            only_in_second[key] = value
    return only_in_first, only_in_second, different_in_both

def compare_env_strings(file1_content, file2_content):
    dict1 = parse_env_string_to_dict(file1_content)
    dict2 = parse_env_string_to_dict(file2_content)
    return compare_dicts(dict1, dict2)

def log_time_delta_from_start(module: str, action: str, node=None):
    logging.info(f">|{module}|{action}{',' + str(node) if node is not None else ''}|Time From start:{to_milliseconds_str(time.time() - config.START_TIME)}")

def set_named_timestamp(action: str, node=None, key=None):
    if key is None:
        key = f"{action}{',' + str(node) if node is not None else ''}"
    config.NAMED_TIMESTAMPS[key] = time.time()

def invalidate_named_timestamp(action: str, node=None, key=None):
    if key is None:
        key = f"{action}{',' + str(node) if node is not None else ''}"
    del config.NAMED_TIMESTAMPS[key]

def log_time_delta_from_start_and_set_named_timestamp(module: str, action: str, node=None, key=None):
    try:
        set_named_timestamp(action, node, key)
        logging.info(f">|{module}|{action}{',' + str(node) if node is not None else ''}|Time from start:{to_milliseconds_str(time.time() - config.START_TIME)}")
    except KeyError:
        logging.error(f"Named timestamp {key} already exists")

def log_time_delta_from_named_timestamp(module: str, action: str, node=None, key=None, invalidate=True):
    try:
        if key is None:
            key = f"{action}{',' + str(node) if node is not None else ''}"
        logging.info(f">|{module}|{action}{',' + str(node) if node is not None else ''}|Time from start:{to_milliseconds_str(time.time() - config.START_TIME)}|Step time:{to_milliseconds_str(time.time() - config.NAMED_TIMESTAMPS[key])}")
        if invalidate:
            invalidate_named_timestamp(action, node, key)
    except KeyError:
        logging.error(f"Named timestamp {key} does not exist")

def to_milliseconds_str(seconds: float) -> str:
    return f"{seconds * 1000:.3f}ms"



def get_all_child_processes(pid):
    try:
        parent = psutil.Process(pid)
    except psutil.NoSuchProcess:
        return []

    children = parent.children(recursive=True)
    parent_of_parent = parent.parent()
    logging.critical("PARENT_PROCESS: " + str(parent_of_parent))
    logging.critical("MAIN_PROCESS: " + str(parent))
    all_processes = [parent] + children
    for process in all_processes:
        logging.critical("PROCESS: " + str(process))
    return all_processes


def kill_process_tree(pid, sig=signal.SIGTERM):
    processes = get_all_child_processes(pid)
    for proc in processes:
        try:
            os.kill(proc.pid, sig)
        except (psutil.NoSuchProcess):
            pass
        except (PermissionError):
            logging.critical("NO PERMISSION")
        except (ProcessLookupError):
            logging.critical("PROCESS LOOKUP ERROR")

    # Check if processes are still alive
    alive_processes = []
    for proc in processes:
        try:
            if proc.is_running():
                alive_processes.append(f"{proc}-({proc.status()})")
        except:
            pass
    return alive_processes


## TODO: Try to move those to PaSh and import them here
def parse_cmd_from_file(file_path: str) -> "tuple[str,list[AstNode]]":
    logging.debug(f'Parsing: {file_path}')
    with open(file_path) as f:
        cmd = f.read()
    asts = analysis.parse_shell_to_asts(file_path)
    return cmd, asts

def parse_edge_line(line: str) -> "tuple[int, int]":
    from_str, to_str = line.split(" -> ")
    return (int(from_str), int(to_str))

def parse_loop_context_line(line: str) -> "tuple[int, list[int]]":
    node_id, loop_contexts_raw = line.split("-loop_ctx-")
    if loop_contexts_raw != "":
        loop_contexts_str = loop_contexts_raw.split(",")
        loop_contexts = [int(loop_ctx) for loop_ctx in loop_contexts_str]
    else:
        loop_contexts = []
    return int(node_id), loop_contexts

def parse_loop_contexts(lines):
    loop_contexts = {}
    for line in lines:
        node_id, loop_ctx = parse_loop_context_line(line)
        loop_contexts[node_id] = loop_ctx
    return loop_contexts

def parse_var_assignment_lines(lines: "list[str]") -> list[int]:
    return {int(line.split("-var")[0]) for line in lines}

def parse_partial_program_order_from_file(file_path: str):
    with open(file_path) as f:
        raw_lines = f.readlines()

    ## Filter comments and remove new lines
    lines = [line.rstrip() for line in raw_lines
            if not line.startswith("#")]

    ## The directory in which cmd_files are
    cmds_directory = str(lines[0])
    logging.debug(f'Cmds are stored in: {cmds_directory}')

    ## The initial env file
    initial_env_file = str(lines[1])

    ## The number of nodes
    number_of_nodes = int(lines[2])
    logging.debug(f'Number of po cmds: {number_of_nodes}')

    basic_blocks_start = lines.index('Basic blocks:') + 1
    basic_blocks_end = lines.index('Basic block edges:')
    basic_block_edges_start = basic_blocks_end + 1
    basic_block_edges_end = lines.index('Loop context:')

    block_num_max = 0
    basic_block_edges = []
    for line in lines[basic_block_edges_start:basic_block_edges_end]:
        from_block, remain = line.split(' -> ')
        to_block, edge_type, aux_info = remain.split(':')
        from_block, to_block = int(from_block), int(to_block)
        edge_type = edge_type.strip()
        if from_block > block_num_max:
            block_num_max = from_block
        if to_block > block_num_max:
            block_num_max = to_block
        basic_block_edges.append((from_block, to_block, edge_type, aux_info))
    hs_prog = HSProg(list(range(block_num_max+1)), basic_block_edges)

    # TODO: rewrite loop_context functions to bb_id functions

    ## The loop context for each node
    loop_context_start = basic_block_edges_end + 1
    loop_context_end = number_of_nodes + loop_context_start
    loop_context_lines = lines[loop_context_start:loop_context_end]
    loop_contexts = parse_loop_contexts(loop_context_lines)
    debug_log(f'Loop contexts: {loop_contexts}')

    var_assignment_lines = int(lines[loop_context_end])
    var_assignment_start = loop_context_end + 1
    var_assignment_end = var_assignment_start + var_assignment_lines
    var_assignment_lines = lines[var_assignment_start:var_assignment_end]
    var_assignments = parse_var_assignment_lines(var_assignment_lines)
    debug_log(f'Var assignments: {var_assignments}')

    ## The rest of the lines are edge_lines
    edge_lines = lines[var_assignment_end:]
    debug_log(f'Edges: {edge_lines}')

    ab_nodes = {}
    for i in range(number_of_nodes):
        file_path = f'{cmds_directory}/{i}'
        cmd, asts = parse_cmd_from_file(file_path)
        loop_ctx = loop_contexts[i]
        var_assignment = (i in var_assignments)
        is_loop_list_change = cmd.startswith('HS_LOOP_LIST=') or cmd.strip() == 'unset HS_LOOP_LIST'
        ab_nodes[NodeId(i)] = Node(NodeId(i), cmd.strip(),
                                   asts, loop_ctx[0],
                                   var_assignment,
                                   is_loop_list_change)
        hs_prog.append_node_to(loop_ctx[0], ab_nodes[NodeId(i)])

    debug_log(str(hs_prog))
    edges = {NodeId(i) : [] for i in range(number_of_nodes)}
    for edge_line in edge_lines:
        from_id, to_id = parse_edge_line(edge_line)
        edges[NodeId(from_id)].append(NodeId(to_id))

    debug_log(f"Nodes|{','.join([str(node) for node in ab_nodes])}")
    debug_log(f"Edges|{edges}")
    return PartialProgramOrder(ab_nodes, edges, hs_prog)

def generate_id() -> int:
    return int(time.time() * 1000000)

# nodes is iterable of node
# edges is dict[node, list[node]]
def invert_graph(nodes, edges):
    graph = {n: [] for n in nodes}
    for from_id, to_ids in edges.items():
        for to_id in to_ids:
            graph[to_id].append(from_id)
    return graph
