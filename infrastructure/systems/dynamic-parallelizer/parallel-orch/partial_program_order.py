from enum import Enum
from node import NodeId, Node, CFGEdgeType, ConcreteNodeId, ConcreteNode, HSProg, HSBasicBlock, HSLoopListContext, loop_iters_do_action, get_loop_list_from_env
import logging
import util
from pathlib import Path
from collections import deque
from executor import run_assignment_and_return_env_file

PROG_LOG =  '[PROG_LOG] '
EVENT_LOG = '[EVENT_LOG] '

def event_log(s):
    logging.info(EVENT_LOG + s)
    pass

def progress_log(s):
    # logging.info(PROG_LOG + s)
    pass

def simulate_loop_iter_env(env, var, loop_list_context, loop_iters):
    loop_list = loop_list_context.get_top()
    try:
        val = loop_list[loop_iters[0]-1]
    except IndexError:
        return env
    cmd = f'{var}="{val}"'
    return run_assignment_and_return_env_file(cmd, env)

class PartialProgramOrder:
    frontier: set  # Set of nodes at the frontier
    # Di: I'm going to ignore this for now and implement the feature without a local data structure
    # Later we can add this back as a caching mechanism to avoid doing RWSet
    # intersections of files all the time
    # run_after: "dict[NodeId, list[Node]]"  # Nodes that should run after certain conditions

    # Mapping of concrete nodes to lists of uncommitted concrete nodes the precedes them.
    # It is the snapshot of the reachable uncommited concrete nodes from prev_concrete_node graph
    # at the time the concrete node enters execution. So if there is fs conflict in them,
    # it needs to be rerun
    to_be_resolved: "dict[NodeId, list[Node]]"
    concrete_nodes: "dict[NodeId, Node]"
    temp_new_env: "tuple[NodeId, str]"

    # this keeps track of the loop list context of the last
    # committed node or loop_list assignment node
    current_loop_list: HSLoopListContext
    
    def __init__(self, abstract_nodes: "dict[NodeId, Node]", edges: "dict[NodeId, list[NodeId]]",
                 hs_prog: HSProg):
        self.hsprog = hs_prog
        self.concrete_nodes: dict[ConcreteNodeId, ConcreteNode] = {}
        self.frontier = set()
        # self.run_after = {}
        # Nodes that we have received "wait" for
        self.canon_exec_order: list[ConcreteNodeId] = list()
        # Nodes that we think should happen, and haven't received "wait" for
        self.spec_exec_order: list[ConcreteNodeId] = list()
        self.to_be_resolved: dict[ConcreteNodeId, list[ConcreteNodeId]] = {}
        self.temp_new_env = None
        self.current_loop_list = HSLoopListContext()

    def commit_node(self, node):
        # Logic to handle committing a node
        node.transition_to_committed()
        # Maybe update dependencies here
        # etc.

    def get_concrete_node(self, concrete_node_id: ConcreteNodeId) -> ConcreteNode:
        return self.concrete_nodes[concrete_node_id]

    def get_all_nodes(self):
        return [node for node in self.concrete_nodes.values()]

    def get_committed_nodes(self):
        return [node for node in self.concrete_nodes.values() if node.is_committed()]

    # def get_ready_nodes(self):
    #     return [(cnid, n) for cnid, n in self.concrete_nodes.items() if n.is_ready()]

    def get_executing_nodes(self):
        return [node for node in self.concrete_nodes.values() if node.is_executing()]

    def get_spec_executing_nodes(self):
        return [node for node in self.concrete_nodes.values() if node.is_spec_executing()]

    def get_executing_normal_and_spec_nodes(self):
        return [node for node in self.concrete_nodes.values() if node.is_executing() or node.is_spec_executing()]

    def get_speculated_nodes(self):
        return [node for node in self.concrete_nodes.values() if node.is_speculated()]

    def get_uncommitted_nodes(self):
        return [node for node in self.concrete_nodes.values() if not node.is_committed()]

    def get_frontier(self):
        return self.frontier

    def log_info(self):
        # logging.info(f"Nodes: {self.concrete_nodes}")
        # logging.info(f"Adjacency: {self.adjacency}")
        # logging.info(f"Inverse adjacency: {self.inverse_adjacency}")
        self.log_state()

    def log_state(self):
        for node in self.concrete_nodes.values():
            progress_log(node.pretty_state_repr())
        progress_log('')
        progress_log('canon: ' + ' '.join([str(cnid) for cnid in self.canon_exec_order]))
        progress_log('spec:  ' + ' '.join([str(cnid) for cnid in self.spec_exec_order]))

    # Determines the next basic block to transition to when the end of the current block is reached.
    def find_next_basic_block(self, basic_block):
        walked_edges = []
        seen_block_ids = set()
        while True:
            if self.hsprog.is_last_block(basic_block) or basic_block.bb_id in seen_block_ids:
                return None, walked_edges
            seen_block_ids.add(basic_block.bb_id)
            edge_type, next_bb = self.hsprog.guess_next_block(basic_block)
            walked_edges.append(edge_type)
            if len(next_bb.nodes) > 0:
                return next_bb, walked_edges
            basic_block = next_bb


    def find_next_non_assignment_node_in_block(self, basic_block, current_node_id: ConcreteNodeId):
        assignment_nodes = []
        next_non_assignment_node = None

        start_index_found = False
        for node in basic_block.nodes:
            if node.id_ == current_node_id.node_id:
                start_index_found = True
                continue  # Start searching from the node after current_node_id
            if start_index_found:
                if node.assignment:
                    assignment_nodes.append(node.id_)
                else:
                    next_non_assignment_node = node
                    break

        return next_non_assignment_node, assignment_nodes


    def find_next_concrete_node_and_gather_var_assignments(self, prev_node: ConcreteNodeId, basic_block: HSBasicBlock):
        next_non_assignment_node, assignment_nodes = self.find_next_non_assignment_node_in_block(basic_block, prev_node)

        # If a non-assignment node is found, prepare its ConcreteNodeId. Otherwise, return None.
        if next_non_assignment_node:
            next_concrete_id = ConcreteNodeId(next_non_assignment_node.id_, prev_node.loop_iters)
        else:
            next_concrete_id = None

        return next_concrete_id, assignment_nodes

    # def make_new_spec_node_with_assignments(self, prev_node: ConcreteNodeId):
    #     basic_block = self.hsprog.find_basic_block(prev_node.node_id)
    #     assert type(prev_node) == ConcreteNodeId, f"Node should be concrete. Given: {type(prev_node)}"
    #     next_concrete_id, assignment_nodes = self.find_next_concrete_node_and_gather_var_assignments(prev_node, basic_block)
    #     # We have to look for the next basic block if we didn't find a non-assignment node in the current block.
    #     while next_concrete_id is None:
    #         next_basic_block, walked_edges = self.find_next_basic_block(basic_block)
    #         if not next_basic_block:
    #             return None, assignment_nodes
    #         next_node = next_basic_block.nodes[0]
    #         next_concrete_id, new_assignments, _ = self.find_next_concrete_node_and_gather_var_assignments(prev_node, next_node, next_basic_block)
    #         assignment_nodes.extend(new_assignments)

    #     return next_concrete_id, assignment_nodes

    def can_speculate_pass(self, prev_node: ConcreteNodeId):
        abstract_node = self.hsprog.find_node(prev_node.node_id)
        concrete_node = self.concrete_nodes[prev_node]
        if not abstract_node.is_loop_list_change():
            return True
        if abstract_node.is_assignment():
            return True
        if concrete_node.is_speculated() or concrete_node.is_committed():
            return True
        return False
    
    def create_concrete_node(self, concrete_node_id: ConcreteNodeId, spec_pre_env: str,
                             loop_list_context: HSLoopListContext):
        if (concrete_node_id in self.concrete_nodes and
            not self.concrete_nodes[concrete_node_id].is_ready()):
            self.concrete_nodes[concrete_node_id].try_reset_to_ready(spec_pre_env)
        else:
            abstract_node = self.hsprog.find_node(concrete_node_id.node_id)
            new_concrete_node = ConcreteNode(concrete_node_id, abstract_node, loop_list_context)
            self.concrete_nodes[concrete_node_id] = new_concrete_node
            new_concrete_node.transition_from_init_to_ready(spec_pre_env)
            if new_concrete_node.command_unsafe():
                new_concrete_node.transition_from_ready_to_unsafe()
    
    # def make_new_spec_node(self, prev_node: ConcreteNodeId):
    #     bb = self.hsprog.find_basic_block(prev_node.node_id)
    #     next_concrete_id = None
    #     assignment_node_ids = []
    #     last_abstract_node_id = prev_node.node_id
    #     walked_edges = []
    #     seen_block_ids = set()
    #     while next_concrete_id is None:
    #         if bb.node_ids[-1] != last_abstract_node_id:
    #             i = bb.node_ids.index(last_abstract_node_id)
    #             next_node_id = bb.node_ids[i+1]
    #             next_node = bb.get_node(next_node_id)
    #             if next_node.is_assignment():
    #                 assignment_node_ids.append(next_node_id)
    #                 last_abstract_node_id = next_node_id
    #                 continue
    #             next_concrete_id = ConcreteNodeId(next_node_id, prev_node.loop_iters)
    #         else:
    #             while True:
    #                 if self.hsprog.is_last_block(bb) or bb.bb_id in seen_block_ids:
    #                     return None
    #                 seen_block_ids.add(bb.bb_id)
    #                 edge_type, next_bb = self.hsprog.guess_next_block(bb)
    #                 walked_edges.append(edge_type)
    #                 bb = next_bb
    #                 if len(next_bb.nodes) > 0:
    #                     break
    #             next_node = next_bb.nodes[0]
    #             next_node_id = next_node.id_
    #             if next_node.is_assignment():
    #                 assignment_node_ids.append(next_node_id)
    #                 last_abstract_node_id = next_node_id
    #                 continue
    #             for edge in walked_edges:
    #                 next_concrete_id = prev_node.do_action(next_node_id, edge)
    #                 prev_node = next_concrete_id
    #     if (next_concrete_id in self.concrete_nodes and
    #         not self.concrete_nodes[next_concrete_id].is_ready()):
    #         self.concrete_nodes[next_concrete_id].reset_to_ready(assignment_node_ids)
    #     else:
    #         new_concrete_node = ConcreteNode(next_concrete_id, next_node)
    #         self.concrete_nodes[next_concrete_id] = new_concrete_node
    #         new_concrete_node.transition_from_init_to_ready(assignment_node_ids)
    #         if new_concrete_node.command_unsafe():
    #             new_concrete_node.transition_from_ready_to_unsafe()
    #     return next_concrete_id

    def make_new_spec_node(self, prev_node_id: ConcreteNodeId):
        bb = self.hsprog.find_basic_block(prev_node_id.node_id)
        prev_node = self.concrete_nodes[prev_node_id]
        abstract_node = self.hsprog.find_node(prev_node_id.node_id)
        if (abstract_node.is_loop_list_change() and not abstract_node.is_assignment()
            and not prev_node.is_committed() and not prev_node.is_speculated()):
            return None
        if prev_node.is_unsafe():
            return None
        prev_loop_list_context = prev_node.loop_list_context
        loop_iters = prev_node_id.loop_iters
        next_concrete_id = None
        last_abstract_node_id = prev_node_id.node_id
        pre_env_file = prev_node.guess_post_env()
        try:
            if len(open(pre_env_file).read()) == 0:
                raise ValueError(pre_env_file)
        except FileNotFoundError:
            breakpoint()
        while True:
            if bb.node_ids[-1] != last_abstract_node_id:
                i = bb.node_ids.index(last_abstract_node_id)
                next_node_id = bb.node_ids[i+1]
                next_node = bb.get_node(next_node_id)
            else:
                while True:
                    if self.hsprog.is_last_block(bb):
                        return None
                    edge_type, next_bb, aux_info = self.hsprog.guess_next_block(
                        bb, loop_iters, prev_loop_list_context)
                    if edge_type == CFGEdgeType.LOOP_TAKEN:
                        pre_env_file = simulate_loop_iter_env(pre_env_file, aux_info,
                                                              prev_loop_list_context,
                                                              loop_iters)
                    loop_iters = loop_iters_do_action(loop_iters, edge_type)
                    bb = next_bb
                    if len(next_bb.nodes) > 0:
                        break
                next_node = bb.nodes[0]
                next_node_id = next_node.id_
            if next_node.is_assignment():
                pre_env_file = next_node.simulate_env(pre_env_file)
                if next_node.is_loop_list_push():
                    prev_loop_list_context = prev_loop_list_context.push(
                        get_loop_list_from_env(pre_env_file))
                elif next_node.is_loop_list_pop():
                    prev_loop_list_context = prev_loop_list_context.pop()
                last_abstract_node_id = next_node_id
                continue
            else:
                cnid = ConcreteNodeId(next_node_id, loop_iters)
                util.debug_log(f'pick {pre_env_file} as pre_env_file')
                new_pre_env_file = util.cp_to_ptmpfile(pre_env_file, 'hs_spec_pre_env')
                util.copy(pre_env_file + '.fds', new_pre_env_file + '.fds')
                util.debug_log(f'copied to {new_pre_env_file}')
                self.create_concrete_node(cnid, new_pre_env_file, prev_loop_list_context)
                return cnid
    
    def get_schedulable_spec_nodes(self) -> list[ConcreteNodeId]:
        schedulable = []
        self.adjust_to_be_resolved_dict()
        for cnid in self.spec_exec_order:
            if self.concrete_nodes[cnid].is_ready() and not self._has_fs_deps(cnid):
                schedulable.append(cnid)
        return schedulable
                        
    def try_schedule_spec_nodes(self, window=2) -> list[ConcreteNodeId]:
        if len(self.canon_exec_order) == 0:
            return
        elif len(self.spec_exec_order) == 0:
            prev_node = self.canon_exec_order[-1]
        else:
            prev_node = self.spec_exec_order[-1]
            prev_node: ConcreteNodeId

        for cnid in self.get_schedulable_spec_nodes():
            self.schedule_spec_work(cnid)

        while len(self.spec_exec_order) < window:
            next_concrete_id = self.make_new_spec_node(prev_node)
            if next_concrete_id is None:
                window = len(self.spec_exec_order)
                break
            self.spec_exec_order.append(next_concrete_id)
            prev_node = self.spec_exec_order[-1]
            if self.concrete_nodes[next_concrete_id].is_ready():
                self.schedule_spec_work(next_concrete_id)

    # def get_prev_nodes(self, concrete_node_id: ConcreteNodeId) -> "list[ConcreteNodeId]":
    #     return self.exec_order[concrete_node_id][:]

    # def get_all_next(self, current_node_id: ConcreteNodeId, visited=None) -> "set[NodeId]":
    #     all_next = set()
    #     def reachable_rec(cur, reachable):
    #         if cur in reachable:
    #             return
    #         reachable.add(cur)
    #         for n in self.get_next_nodes(cur):
    #             reachable_rec(n, reachable)
    #     for n in self.get_next_nodes(current_node_id):
    #         reachable_rec(n, all_next)
    #     return all_next


    def get_all_previous(self, current_node_id: ConcreteNodeId):
        if current_node_id in self.canon_exec_order:
            i = self.canon_exec_order.index(current_node_id)
            return self.canon_exec_order[:i]
        elif current_node_id in self.spec_exec_order:
            i = self.spec_exec_order.index(current_node_id)
            return self.canon_exec_order[:] + self.spec_exec_order[:i]
        else:
            assert False

    def get_all_hypothetical_previous(self, current_node_id: ConcreteNodeId):
        if current_node_id in self.canon_exec_order:
            i = self.canon_exec_order.index(current_node_id)
            return self.canon_exec_order[:i]
        elif current_node_id in self.spec_exec_order:
            i = self.spec_exec_order.index(current_node_id)
            return self.canon_exec_order[:] + self.spec_exec_order[:i]
        else:
            return self.canon_exec_order[:] + self.spec_exec_order[:]

    # TODO: fixme
    # def get_all_next_uncommitted(self, node_id: NodeId) -> "set[NodeId]":
    #     next = self.get_all_next(node_id)
    #     return set([node for node in next if not self.concrete_nodes[node].is_committed()])

    def get_all_previous_uncommitted(self, concrete_node_id: ConcreteNodeId) -> "set[ConcreteNodeId]":
        previous = self.get_all_previous(concrete_node_id)
        cnid_eid_tuples = []
        for cnid in previous:
            concrete_node = self.concrete_nodes[cnid]
            if concrete_node.is_speculated():
                cnid_eid_tuples.append((concrete_node.cnid, concrete_node.exec_id))
            elif not concrete_node.is_committed():
                cnid_eid_tuples.append(concrete_node.cnid)
        return set(cnid_eid_tuples)

    def adjust_to_be_resolved_dict_entry(self, concrete_node_id: ConcreteNodeId):
        node = self.concrete_nodes.get(concrete_node_id)
        if node.is_committed():
            self.to_be_resolved[concrete_node_id] = []
        elif node.is_ready() and (concrete_node_id in self.spec_exec_order or
                                  concrete_node_id in self.canon_exec_order):
            self.to_be_resolved[concrete_node_id] = self.get_all_previous_uncommitted(concrete_node_id)

    def init_to_be_resolved_dict(self):
        for node_id in self.concrete_nodes:
            self.adjust_to_be_resolved_dict_entry(node_id)

    def adjust_to_be_resolved_dict(self):
        # TODO: this design seems to require the function to be called
        # each time before a node entering EXECUTING or SPEC_EXECUTING
        # to be optimal (that is, it might keep more things in the list).
        # It's safe as is so I'm not touching it.
        for node_id in self.to_be_resolved.keys():
            self.adjust_to_be_resolved_dict_entry(node_id)

    #TODO: Add partial order invariant checks
    def valid(self):
        return True

    def fetch_fs_actions(self):
        for node in self.get_executing_normal_and_spec_nodes():
            util.overhead_log(f"TRACE_FETCHING|{node.cnid}")
            node.gather_fs_actions()
            util.overhead_log(f"TRACE_FETCHING_END|{node.cnid}")

    def _has_fs_deps(self, concrete_node_id: ConcreteNodeId):
        node_of_interest : ConcreteNode = self.get_concrete_node(concrete_node_id)
        for dep_entry in self.to_be_resolved[concrete_node_id]:
            if isinstance(dep_entry, tuple):
                nid, eid = dep_entry
                if eid == self.get_concrete_node(nid).exec_id:
                    continue
            else:
                nid = dep_entry    
            node: ConcreteNode = self.get_concrete_node(nid)
            if node.get_rw_set().has_conflict(node_of_interest.get_rw_set()):
                return True
        return False

    # TODO: It's currently designed this way to avoid reading trace file all the time
    # When we have complex caching code for this we can make this go away
    def has_fs_deps(self, concrete_node_id: ConcreteNodeId):
        self.fetch_fs_actions()
        self._has_fs_deps(concrete_node_id)

    def schedule_spec_work(self, concrete_node_id: ConcreteNodeId):
        concrete_node = self.get_concrete_node(concrete_node_id)
        self.adjust_to_be_resolved_dict_entry(concrete_node_id)
        speculated_nodes = [self.concrete_nodes[cnid] for cnid in self.spec_exec_order if
                            self.concrete_nodes[cnid].is_speculated()]
        self.get_concrete_node(concrete_node_id).start_spec_executing(concrete_node.spec_pre_env,
                                                                      speculated_nodes)

    def simulate_var_assignments(self, env, assignments: "list[NodeId]"):
        for assignment in assignments:
            assignment_node: Node = self.hsprog.find_node(assignment)
            env = run_assignment_and_return_env_file(assignment_node.cmd, env)
        return env

    def reset_speculation(self):
        event_log(f"reset speculation")
        for cnid in self.spec_exec_order:
            self.concrete_nodes[cnid].try_reset_to_ready()
        self.spec_exec_order = []
    
    ### external handler events ###

    def schedule_work(self, concrete_node_id: ConcreteNodeId, env_file: str):
        self.get_concrete_node(concrete_node_id).start_executing(env_file)

    def handle_complete(self, concrete_node_id: ConcreteNodeId, has_pending_wait: bool,
                        current_env: str):
        event_log(f"handle_complete {concrete_node_id}")
        node = self.get_concrete_node(concrete_node_id)
        util.debug_log(f"outfds: {node.exec_ctxt.outfds}")
        # TODO: make collect_result a state transition and make more states
        is_killed, runtime_finished = node.collect_result()
        if is_killed:
            node.reset_to_ready()
            if has_pending_wait:
                node.start_executing(current_env)
            return
        if not runtime_finished:
            node.reset_to_ready()
            node.transition_from_ready_to_unsafe()
            return
        if node.is_executing():
            node.commit_frontier_execution()
            self.current_loop_list = node.loop_list_context
            self.adjust_to_be_resolved_dict()
        elif node.is_spec_executing():
            if self.has_fs_deps(concrete_node_id):
                node.reset_to_ready()
                # otherwise it stays in ready state and waits to be scheduled by the scheduler
                if has_pending_wait:
                    node.start_executing(current_env)
            else:
                node.finish_spec_execution()
                if has_pending_wait:
                    self.current_loop_list = node.loop_list_context
                    node.commit_speculated()
                    util.good_log(f"{concrete_node_id} speculation committed")
                    self.adjust_to_be_resolved_dict()
        else:
            assert False

    def reset_succeeding_nodes(self, node_id: NodeId, env_file: str):
        # TODO: fixme
        pass
        # for uncommitted_node_id in self.get_all_next(node_id):
        #     uncommitted_node = self.get_concrete_node(uncommitted_node_id)
        #     if uncommitted_node.is_spec_executing():
        #         uncommitted_node.reset_to_ready()
        #     # uncommitted_node.start_spec_executing(env_file)

    def finish_wait_unsafe(self, concrete_node_id: ConcreteNodeId, env: str):
        node = self.concrete_nodes[concrete_node_id]
        node.spec_pre_env = env

    # Returns whether handle_wait should be called.
    # This function exists because handle_wait always guarantees the creation of
    # a concrete_node, but something needs to be handled without concrete_node
    # (e.g. loop list assignment)
    def pre_handle_wait(self, concrete_node_id: ConcreteNodeId, env_file: str):
        node = self.hsprog.find_node(concrete_node_id.node_id)
        if node.is_assignment():
            if node.is_loop_list_change():
                self.current_loop_list = node.simulate_loop_list(env_file, self.current_loop_list)
            return False
        else:
            return True
        
    def handle_wait(self, concrete_node_id: ConcreteNodeId, env_file: str):
        event_log(f"handle_wait {concrete_node_id}")
        util.env_log(f"wait for {concrete_node_id}, incoming env {env_file}")

        if len(self.spec_exec_order) and concrete_node_id == self.spec_exec_order[0]:
            self.spec_exec_order.pop(0)
        else:
            self.reset_speculation()

        if not concrete_node_id in self.concrete_nodes:
            abstract_node = self.hsprog.find_node(concrete_node_id.node_id)
            new_concrete_node = ConcreteNode(concrete_node_id, abstract_node, self.current_loop_list)
            new_concrete_node.transition_from_init_to_ready(env_file)
            if new_concrete_node.command_unsafe():
                new_concrete_node.transition_from_ready_to_unsafe()
            self.concrete_nodes[concrete_node_id] = new_concrete_node

        self.canon_exec_order.append(concrete_node_id)

        node = self.get_concrete_node(concrete_node_id)
        # Invalid state check
        if node.is_committed() or node.is_initialized():
            logging.error(f'Error: Node {concrete_node_id} is in an invalid state: {node.state}')
            raise Exception(f'Error: Node {concrete_node_id} is in an invalid state: {node.state}')

        # Set the temp_new_env in case we have to restart an eagerly killed node
        self.temp_new_env = (concrete_node_id, env_file)
        
        if node.is_ready():
            event_log(f"schedule {concrete_node_id}")
            node.start_executing(env_file)
        elif node.is_unsafe():
            pass
        elif node.is_stopped():
            if node in self.get_frontier():
                logging.info(f'Node {concrete_node_id} is stopped and in the frontier.')
                node.transition_from_stopped_to_executing(env_file)
            else:
                logging.info(f'Node {concrete_node_id} is stopped but not in the frontier.')
        elif node.is_speculated():
            # Check if env conflicts exist
            if node.has_env_conflict_with(env_file):
                util.debug_log(f'prev_env: {node.exec_ctxt.pre_env_file}, real: {env_file}')
                node.reset_to_ready(loop_list_context=self.current_loop_list)
                node.start_executing(env_file)
                self.reset_speculation()
            # Optimization: It would make sense to perform the checks independently,
            # and if fs conflict, then update the run after dict.
            elif self.has_fs_deps(concrete_node_id):
                node.reset_to_ready(loop_list_context=self.current_loop_list)
                node.start_executing(env_file)
            else:
                node.commit_speculated()
                util.good_log(f"{concrete_node_id} speculation committed")
                self.current_loop_list = node.loop_list_context
                self.adjust_to_be_resolved_dict()
        elif node.is_executing():
            if node.has_env_conflict_with(env_file):
                node.reset_to_ready(loop_list_context=self.current_loop_list)
                node.start_executing(env_file)
                self.reset_speculation()
        elif node.is_spec_executing():
            if node.has_env_conflict_with(env_file):
                node.reset_to_ready(loop_list_context=self.current_loop_list)
                node.start_executing(env_file)
                self.reset_speculation()
        else:
            logging.error(f'Error: Node {concrete_node_id} is in an invalid state: {node.state}')
            raise Exception(f'Error: Node {concrete_node_id} is in an invalid state: {node.state}')

    def eager_fs_killing(self):
        event_log("try to eagerly kill conflicted speculation")
        to_be_killed: "list[ConcreteNode]" = []
        self.fetch_fs_actions()
        for node in self.get_all_nodes():
            if ((node.is_speculated() or node.is_spec_executing())
                and self._has_fs_deps(node.cnid)):
                to_be_killed.append(node)
        for node in to_be_killed:
            node.reset_to_ready()
            # If we don't restart the node with pending wait here, the scheduler will hang
            if node.cnid==self.temp_new_env[0]:
                node.start_executing(self.temp_new_env[1])
