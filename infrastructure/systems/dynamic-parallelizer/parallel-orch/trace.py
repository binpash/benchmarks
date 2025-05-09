import re
import sys
import os
from typing import Tuple
from enum import Enum
import logging
from copy import deepcopy


class Ref(Enum):

    STDIN = sys.stdin
    STDOUT = os.path.abspath(os.sep)
    STDERR = sys.stderr
    ROOT = os.path.abspath(os.sep)
    # Not sure this is always correct
    # but it doesn't affect the results
    CWD = os.getcwd()
    LAUNCH_EXE = ""


class PathRef:

    def __init__(self, ref, path, permissions, no_follow, env):
        self.ref = PathRefKey(env, ref)
        self.path = path
        self.is_read, self.is_write, self.is_exec = self.resolve_permissions(
            permissions)
        self.is_nofollow = no_follow

    def resolve_permissions(self, permissions: str):
        if "r" in permissions:
            is_read = True
        else:
            is_read = False
        if "w" in permissions:
            is_write = True
        else:
            is_write = False
        if "x" in permissions:
            is_exec = True
        else:
            is_exec = False
        return is_read, is_write, is_exec

    def __str__(self):
        return f"PathRef({self.ref}, {self.path}, {'r' if self.is_read else '-'}{'w' if self.is_write else '-'}{'x' if self.is_exec else '-'} {'no follow' if self.is_nofollow else ''})"
    
    def __repr__(self) -> str:
        return self.__str__()

    def get_resolved_path(self):
        
        if isinstance(self.ref, PathRef):
            self.ref = self.ref.get_resolved_path()
        
        # Remove dupliate prefixes
        if not self.path.startswith("/"):
            modified_path = "/" + self.path
        else:
            modified_path = self.path
        commonprefix = os.path.commonprefix([self.ref, modified_path])
        ref_without_prefix = self.ref.replace(commonprefix, "", 1)
        path_without_prefix = modified_path.replace(commonprefix, "", 1)

        if path_without_prefix.startswith("/"):
            path_without_prefix = path_without_prefix.replace("/", "", 1)

        return os.path.join(commonprefix, ref_without_prefix, path_without_prefix).replace("/./", "/")

    

class PathRefKey:

    def __init__(self, env, lhs_ref) -> None:
        self.env = env
        self.lhs_ref = lhs_ref

    def __eq__(self, other):
        return (self.env, self.lhs_ref) == (other.env, other.lhs_ref)

    def __ne__(self, other) -> bool:
        return not (self == other)

    def __hash__(self):
        return hash((self.env, self.lhs_ref))

    def __str__(self):
        return f"Key({self.lhs_ref}@{self.env})"
    
    def __repr__(self) -> str:
        return self.__str__()


class ExpectResult():

    def __init__(self, ref, result):
        self.ref = ref
        self.result = result

    def __str__(self, ref, result):
        return f"ExpectResult({self.ref}, {self.result})"


class PipeRef:
    
    def __init__(self, lhs_ref, env):
        self.ref = PathRefKey(env, lhs_ref)
        
    def __str__(self):
        return f"PipeRef({self.lhs_ref})"


def log_resolved_trace_items(resolved_dict):
    for k, v in resolved_dict.items():
        try:
            logging.debug(f" {k}: {v}")
        except:
            logging.debug(f'{k}: {v}')


def remove_command_redir(cmd):
    return cmd.split(">")[0].rstrip()


def remove_command_prefix(line) -> str:
    return line.split(f"]: ")[1].rstrip()


def get_command_prefix(line):
    return line.split(f"]: ")[0].lstrip("[")


def is_no_command_prefix(line):
    return "No Command" in get_command_prefix(line)


def is_new_path_ref(trace_item):
    return "PathRef" in trace_item

def is_pipe_ref(trace_item):
    return "PipeRef" in trace_item


def get_path_ref_id(trace_item):
    return trace_item.split("=")[0].strip()


def get_path_ref_open_config(trace_item):
    assert (is_new_path_ref(trace_item))
    # WARNING: HACK
    open_config_suffix = trace_item.split(", ")[2]

    open_config = re.split('\(|\)', open_config_suffix)[0]
    # WARNING: HACK: hard-coded replacement
    open_config = open_config.replace("truncate create", "").rstrip()
    return open_config


def get_path_ref_no_follow(trace_item):
    return "nofollow" in trace_item


def is_path_ref_read(trace_item: PathRef):
    return trace_item.is_read


def is_path_ref_write(trace_item: PathRef):
    return trace_item.is_write


def is_path_ref_execute(trace_item: PathRef):
    return trace_item.is_write


def is_path_ref_empty(trace_item: PathRef):
    return not trace_item.is_read and not trace_item.is_write and not trace_item.is_exec


def get_path_ref_name(trace_item):
    assert (is_new_path_ref(trace_item))
    return trace_item.split(", ")[1].replace('"', '')


def get_path_ref_ref(trace_item):
    assert (is_new_path_ref(trace_item))
    return trace_item.split(", ")[0].split("(")[1]


def is_no_command_prefix(line):
    if line.startswith(f"[No Command"):
        return True
    return False


def is_launch(line):
    return "Launch(" in line


def parse_launch_command(trace_item):
    assignment_prefix = trace_item.split("], ")[0].split(
        "([Command ")[1].rstrip("]").strip()
    assignment_suffix = ", ".join(trace_item.split("], ")[1:]).strip()
    assignment_string = assignment_suffix[1:-2].split(",")
    assignments = [(x.split("=")) for x in assignment_string]
    return assignment_prefix, assignments


def get_lauch_name(trace_item):
    assert (is_launch(trace_item))
    launch_name_dirty = trace_item.split("],")[0]
    launch_name = launch_name_dirty.split("Command ")[1]
    return launch_name


def get_no_command_ref_id(trace_item):
    return trace_item.split("=")[0].strip()


def get_no_command_ref_ref(trace_item):
    return trace_item.split("=")[1].strip()


def is_prefix_of_cmd(line, prefix):
    if prefix is not None and prefix in get_command_prefix(line):
        return True
    return False


def is_expect_result(trace_item):
    return "ExpectResult(" in trace_item


def parse_expect_result(trace_item):
    return trace_item.lstrip("ExpectResult(").split(")")[0].split(", ")

def parse_pipe_ref(trace_item):
    return trace_item.split("] = ")[0].lstrip("[").split(", ")

def parse_launch(refs_dict, keys_order, env, line) -> None:
    assignment_prefix, assignments = parse_launch_command(
        remove_command_prefix(line))
    for assignment in assignments:
        lhs_ref = PathRefKey(assignment_prefix, assignment[0].strip())
        rhs_ref = PathRefKey(env, assignment[1].strip())
        refs_dict[lhs_ref] = refs_dict[rhs_ref]
        keys_order.append(lhs_ref)

def add_ref_to_refs_dict(refs_dict, keys_order, lhs_ref, ref):
    refs_dict[lhs_ref] = ref
    keys_order.append(lhs_ref)


def parse_final_refs(refs_dict, keys_order, env, line) -> None:
    line = remove_command_prefix(line)
    path_ref_id = get_no_command_ref_id(line).strip()
    lhs_ref = PathRefKey(env, path_ref_id)
    rhs_ref = get_no_command_ref_ref(line)
    if rhs_ref == "CWD":
        add_ref_to_refs_dict(refs_dict, keys_order, lhs_ref, Ref.CWD)
    elif rhs_ref == "ROOT":
        add_ref_to_refs_dict(refs_dict, keys_order, lhs_ref, Ref.ROOT)
    elif rhs_ref == "STDERR":
        add_ref_to_refs_dict(refs_dict, keys_order, lhs_ref, Ref.STDERR)
    elif rhs_ref == "STDIN":
        add_ref_to_refs_dict(refs_dict, keys_order, lhs_ref, Ref.STDIN)
    elif rhs_ref == "STDOUT":
        add_ref_to_refs_dict(refs_dict, keys_order, lhs_ref, Ref.STDOUT)
    elif rhs_ref == "LAUNCH_EXE":
        add_ref_to_refs_dict(refs_dict, keys_order, lhs_ref, Ref.LAUNCH_EXE)


def parse_new_path_ref(refs_dict, keys_order, env, line):
    line = remove_command_prefix(line).strip()
    lhs_ref = PathRefKey(env, get_path_ref_id(line).strip())
    ref = get_path_ref_ref(line)
    name = get_path_ref_name(line)
    open_config = get_path_ref_open_config(line)
    no_follow = get_path_ref_no_follow(line)
    path_ref = PathRef(ref, name, open_config, no_follow, env)
    refs_dict[lhs_ref] = path_ref
    keys_order.append(lhs_ref)

def parse_expect_result_item(expect_result_dict, env, line):
    line = remove_command_prefix(line).strip()
    path_ref_id, result = parse_expect_result(line)
    lhs_ref = PathRefKey(env, path_ref_id)
    expect_result_dict[lhs_ref] = ExpectResult(lhs_ref, result)
    
def parse_pipe_ref_item(refs_dict, keys_order, env, line):
    line = remove_command_prefix(line).strip()
    # lhs_ref, rhs_ref = parse_pipe_ref(line)
    # Warning HACK: This is a hack to get the correct lhs_ref
    # we are probably ok with this because it.
    rhs_ref, lhs_ref = parse_pipe_ref(line)
    lhs_key = PathRefKey(env, lhs_ref)
    lhs_key_rev = PathRefKey(env, rhs_ref)
    pipe_ref = PipeRef(rhs_ref, env)
    pipe_ref_rev = PipeRef(lhs_ref, env)
    refs_dict[lhs_key] = pipe_ref.ref
    refs_dict[lhs_key_rev] = pipe_ref_rev.ref
    keys_order.append(lhs_key)

def parse_rw_sets(trace_object) -> None:
    # logging.trace("".join(trace_object))
    refs_dict = {}
    expect_result_dict = {}
    keys_order = []
    # In the first iteration, we get the refs
    for line in trace_object:
        # This branch will always execute first
        env = get_command_prefix(line).lstrip("Command").strip()
        if is_no_command_prefix(line):
            if is_launch(line):
                parse_launch(refs_dict, keys_order, env, line)
            elif " = " in line:
                parse_final_refs(refs_dict, keys_order, env, line)
        # Parses Launch(...)
        elif is_launch(line):
            parse_launch(refs_dict, keys_order, env, line)
        # Parses PathRef(...)
        elif is_new_path_ref(line):
            parse_new_path_ref(refs_dict, keys_order, env, line)
        # Parses PipeRef
        elif is_pipe_ref(line):
            parse_pipe_ref_item(refs_dict, keys_order, env, line)
        # Parses ExpectResult(...)
        elif is_expect_result(line):
            parse_expect_result_item(expect_result_dict, env, line)
    return refs_dict, expect_result_dict, keys_order


def traverse_path_ref(refs_dict: dict, ref: PathRef):
    if isinstance(ref, PathRef) and not ref.is_nofollow and isinstance(refs_dict[ref.ref], PathRef):
        return traverse_path_ref(refs_dict, refs_dict[ref.ref])
    else:
        return ref.ref


def resolve_rw_set_refs(refs_dict):
    for ref_item, ref in refs_dict.items():
        if isinstance(ref, PathRef):
            refs_dict[ref_item].ref = traverse_path_ref(refs_dict, ref)
    return refs_dict


def replace_path_ref_terminal_nodes(refs_dict: dict):
    refs_dict_new = {}
    for i, ref in refs_dict.items():
        if isinstance(ref, PathRef):
            
            # HACK: This is hard-coded stdout
            if ref.path == "" and ref.is_nofollow:
                continue
            else:
                # If ref of ref is string, it means that we reached a terminal node.
                if isinstance(ref.ref, str):
                    pass
                elif ref.ref not in refs_dict:
                    key = PathRefKey("No Command", "r1")
                    ref.ref = refs_dict[key].value
                else:
                    
                    if isinstance(refs_dict[ref.ref], Ref):
                        ref.ref = deepcopy(str(refs_dict[ref.ref].value))
                    else:
                        ref.ref = deepcopy(refs_dict[ref.ref])
                assert(i not in refs_dict_new)
                refs_dict_new[i] = deepcopy(ref)
    return refs_dict_new


def resolve_dir_rw_paths(read_set, write_set, dir_set):
    prefix = os.path.commonprefix(dir_set)
    suffixes = [dir.replace(prefix, "") for dir in dir_set]
    # Warning: HACK
    dir_string = prefix
    for dir in suffixes:
        dir_string = os.path.join(dir_string, dir)
        to_add = os.path.join(prefix, dir_string)
        if to_add.endswith("/"):
            write_set.append(to_add)
        else:
            write_set.append(to_add + "/")


def resolve_dir_accesses_from_parsed_items(resolved_dict_replaced, expect_result_dict,
                                           key, previous_key, resolved_trace_object,
                                           write_set, dir_set):
    if previous_key in resolved_dict_replaced:
        previous_resolved_trace_object = resolved_dict_replaced[previous_key]
        relevant_current_expect_result = expect_result_dict.get(previous_key)
        relevant_previous_expect_result = expect_result_dict.get(key)
        if relevant_current_expect_result is not None and relevant_previous_expect_result is not None:
            if isinstance(previous_resolved_trace_object, PathRef) and \
                    is_path_ref_write(previous_resolved_trace_object) and \
                    relevant_current_expect_result.result == "SUCCESS" and \
                    relevant_previous_expect_result.result == "SUCCESS":
                dir_set.append(resolved_trace_object.get_resolved_path())
                write_set.pop()


def resolve_rw_sets_from_parsed_items(resolved_dict_replaced, expect_result_dict, keys_order):
    read_set = set()
    write_set = []
    dir_set = []
    for i, key in enumerate(keys_order):
        if key not in resolved_dict_replaced:
            continue
        resolved_trace_object = resolved_dict_replaced[key]
        if isinstance(resolved_trace_object, Ref):
            continue
        # WARNING: HACK: We need to make sure this condition does not lead to missed dependencies
        # KK 2023-05-03: I don't see where this is useful
        # if resolved_trace_object.get_resolved_path().startswith(os.path.abspath('/tmp/pash_spec')) or 
        #     continue

        ## Each separate node has a different /dev/tty (even though they all seem to write to it)
        ##  so we never want to keep /dev/tty in the read-write sets of any node.
        ## We take care of writes to stdout and stderr elsewhere in the code.
        ## TODO: Generalize this to other special files too (make a global list of such files)
        ## TODO: We actually want to add these to the read-write sets, but then don't take them
        ##       into account when doing the resolution. Trace should not have any scheduling
        ##       logic, it should just parse the trace.
        if resolved_trace_object.get_resolved_path().startswith(os.path.abspath('/dev/tty')):
            continue
        if is_path_ref_read(resolved_trace_object):
            read_set.add(resolved_trace_object.get_resolved_path())
        if is_path_ref_write(resolved_trace_object):
            write_set.append(resolved_trace_object.get_resolved_path())
        # This is a sign that a directory declaration might exist
        if is_path_ref_empty(resolved_trace_object) and i > 0:
            pass
    resolve_dir_rw_paths(read_set, write_set, dir_set)
    return read_set, write_set

# Parse the trace object and gather rw sets for this command
# TODO: PathRefs now also contain environments.
#       Figure out a way to resolve ref_id+env key combinations.


def parse_and_gather_cmd_rw_sets(trace_object) -> Tuple[set, set]:
    refs_dict, expect_result_dict, keys_order = parse_rw_sets(trace_object)
    resolved_dict = resolve_rw_set_refs(refs_dict)
    resolved_dict_replaced = replace_path_ref_terminal_nodes(resolved_dict)
    read_set, write_set = resolve_rw_sets_from_parsed_items(
        resolved_dict_replaced, expect_result_dict, keys_order)
    return read_set, set(write_set)


def parse_exit_code(trace_object) -> int:
    for line in reversed(trace_object):
        if "Exit(" in line:
            return int(line.split("Exit(")[1].rstrip(")\n"))

# Trace can be called as a script with the trace file to analyze as an argument
def main():
    logging.basicConfig(level=logging.DEBUG)
    trace_file = sys.argv[1]
    with open(trace_file, "r") as f:
        trace_object = f.readlines()
    read_set, write_set = parse_and_gather_cmd_rw_sets(trace_object)
    print("Read set:")
    for r in read_set:
        print(r)
    print("Write set:")
    for w in write_set:
        print(w)
    print("Exit code:")
    print(parse_exit_code(trace_object))
    
if __name__ == "__main__":
    main()