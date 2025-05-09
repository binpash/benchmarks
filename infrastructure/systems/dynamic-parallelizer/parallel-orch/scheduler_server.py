import argparse
import logging
import signal
import util
import config
import os
from partial_program_order import PartialProgramOrder, NodeId
from node import LoopStack, ConcreteNodeId

##
## A scheduler server
##

def handler(signum, frame):
    logging.debug(f'Signal: {signum} caught')
    shutdown()

signal.signal(signal.SIGTERM, handler)

def parse_args():
    parser = argparse.ArgumentParser(add_help=False)
    ## TODO: Import the arguments so that they are not duplicated here and in orch
    parser.add_argument("-d", "--debug-level",
                        type=int,
                        default=0,
                        help="Set debugging level")
    parser.add_argument("-f", "--log_file",
                        type=str,
                        default=None,
                        help="Set logging output file. Default: stdout")
    parser.add_argument("--sandbox-killing",
                        action="store_true",
                        default=False,
                    help="Kill any running overlay instances before commiting to the lower layer")
    parser.add_argument("--speculate-immediately",
                    action="store_true",
                    default=False,
                    help="Speculate immediately instead of waiting for the first Wait message.")
    parser.add_argument("--window",
                        type=int,
                        default=5,
                        help="Number of commands to speculate.")

    args, unknown_args = parser.parse_known_args()
    return args

def init():
    args = parse_args()
    # config.set_config_globals_from_pash_args(args)
    return args

def success_response(string):
    return f'OK: {string}\n'

def unsafe_response(string):
    return f'UNSAFE: {string}\n'

def error_response(string):
    return f'ERROR: {string}\n'


class Scheduler:
    """ Schedules a partial order of commands to run out-of-order
    Flow:
        input cmd ->
                    |   Daemon Start -> Receive whens tarting
                    |   Init -> Read the partial order from a file
                    |   CommandExecComplete -> A command completed its execution
                    |   Wait -> The JIT component waits for the results of a specific command
                    |   Done -> We are done
    """
    window: int  # Integer representing the window
    latest_env: str # This variable should be initialized by the first wait, and always have a value since

    def __init__(self, socket_file, window):
        self.window = window
        self.done = False
        self.socket = util.init_unix_socket(socket_file)
        ## A map containing connections for node_ids that are waiting for a response
        self.waiting_for_response = {}
        self.partial_program_order = None

    def handle_init(self, input_cmd: str):
        assert(input_cmd.startswith("Init"))
        partial_order_file = input_cmd.split(":")[1].rstrip()
        util.debug_log(f'Scheduler: Received partial_order_file: {partial_order_file}')
        self.partial_program_order = util.parse_partial_program_order_from_file(partial_order_file)
        util.debug_log(str(self.partial_program_order.hsprog))

    def handle_wait(self, input_cmd: str, connection):
        concrete_node_id, env_file = self.__parse_wait(input_cmd)
        self.waiting_for_response[concrete_node_id] = connection
        util.debug_log(f'Scheduler: Received wait message - {concrete_node_id}.')
        self.latest_env = env_file
        if self.partial_program_order.pre_handle_wait(concrete_node_id, env_file):
            self.partial_program_order.handle_wait(concrete_node_id, env_file)
            concrete_node = self.partial_program_order.get_concrete_node(concrete_node_id)
            if concrete_node.is_committed():
                self.respond_to_pending_wait(concrete_node_id)
            elif concrete_node.is_unsafe():
                util.debug_log(f'unsafe {concrete_node_id}')
                self.partial_program_order.finish_wait_unsafe(concrete_node_id, env_file)
                self.respond_to_wait_on_unsafe(concrete_node_id)
        else:
            util.debug_log(f'ignoring var assignment {concrete_node_id}')
            self.respond_to_wait_on_unsafe(concrete_node_id)

    def process_next_cmd(self):
        connection, input_cmd = util.socket_get_next_cmd(self.socket)

        if(input_cmd.startswith("Init")):
            connection.close()
            self.handle_init(input_cmd)
        elif (input_cmd.startswith("Daemon Start") or input_cmd == ""):
            util.debug_log(f'Scheduler: Received daemon start message.')
            connection.close()
        elif (input_cmd.startswith("CommandExecComplete:")):
            node_id, exec_id, sandbox_dir, trace_file = self.__parse_command_exec_x(input_cmd)
            connection.close()
            if self.partial_program_order.get_concrete_node(node_id).exec_id == exec_id:
                util.debug_log(f'Scheduler: Received command exec complete message - {node_id}.')
                self.partial_program_order.handle_complete(node_id, node_id in self.waiting_for_response, self.latest_env)

                if self.partial_program_order.get_concrete_node(node_id).is_committed():
                    self.respond_to_pending_wait(node_id)
                elif self.partial_program_order.get_concrete_node(node_id).is_unsafe() and \
                     node_id in self.waiting_for_response:
                    self.partial_program_order.finish_wait_unsafe(node_id, self.latest_env)
                    self.respond_to_wait_on_unsafe(node_id)
            else:
                util.debug_log(f'Scheduler: Received command exec complete message for a killed instance, ignoring - {node_id}.')
        elif (input_cmd.startswith("Wait")):
            self.handle_wait(input_cmd, connection)
        elif (input_cmd.startswith("Done")):
            util.socket_respond(connection, success_response("All finished!"))
            self.partial_program_order.log_info()
            self.done = True
        else:
            logging.error(error_response(f'Error: Unsupported command: {input_cmd}'))
            raise Exception(f'Error: Unsupported command: {input_cmd}')

    def respond_to_frontend_core(self, node_id: NodeId, response: str):
        assert(node_id in self.waiting_for_response)
        ## Get the connection that we need to respond to
        connection = self.waiting_for_response.pop(node_id)
        util.socket_respond(connection, response)
        connection.close()

    def respond_to_wait_on_unsafe(self, node_id: ConcreteNodeId):
        response = unsafe_response('')
        self.respond_to_frontend_core(node_id, response)

    def respond_to_pending_wait(self, node_id: ConcreteNodeId):
        logging.debug(f'Responding to pending wait for node: {node_id}')
        ## Get the completed node info
        node = self.partial_program_order.get_concrete_node(node_id)
        msg = '{} {} {}'.format(*node.execution_outcome())
        util.debug_log(f'outcome for node {node_id} is {node.execution_outcome()}')
        response = success_response(msg)

        ## Send the response
        self.respond_to_frontend_core(node_id, response)

    def __parse_wait(self, input_cmd: str) -> "tuple[ConcreteNodeId, str]":
        try:
            node_id_component, loop_iter_counter_component, pash_runtime_vars_file_component = input_cmd.rstrip().split("|")
            node_id = NodeId(int(node_id_component.split(":")[1].rstrip()))
            loop_counters_str = loop_iter_counter_component.split(":")[1].rstrip()
            pash_env_filename = pash_runtime_vars_file_component.split(":")[1].rstrip()
            if loop_counters_str == "None":
                return ConcreteNodeId(node_id), pash_env_filename
            else:
                loop_counters = [int(cnt) for cnt in loop_counters_str.split("-")]
                return ConcreteNodeId(node_id, loop_counters), pash_env_filename
        except:
            raise Exception(f'Parsing failure for line: {input_cmd}')

    def __parse_command_exec_x(self, input_cmd: str) -> "tuple[int, int]":
        try:
            components = input_cmd.rstrip().split("|")
            command_id = ConcreteNodeId.parse(components[0].split(":")[1])
            exec_id = int(components[1].split(":")[1])
            sandbox_dir = components[2].split(":")[1]
            trace_file = components[3].split(":")[1]
            return command_id, exec_id, sandbox_dir, trace_file
        except:
            raise Exception(f'Parsing failure for line: {input_cmd}')


    def schedule_work(self):
        self.partial_program_order.try_schedule_spec_nodes(self.window)

    def run(self):
        ## The first command should be the daemon start
        self.process_next_cmd()

        ## The second command should be the partial order init
        self.process_next_cmd()

        self.partial_program_order.log_state()
        while not self.done:
            self.process_next_cmd()
            self.partial_program_order.log_state()
            self.schedule_work()
            self.partial_program_order.log_state()
            self.partial_program_order.eager_fs_killing()
            self.partial_program_order.log_state()
        self.socket.close()
        self.shutdown()

    def shutdown(self):
        ## There may be races since this is called through the signal handling
        logging.debug("PaSh-Spec scheduler is shutting down...")
        logging.debug("PaSh-Spec scheduler shut down successfully...")
        self.terminate_pending_commands()

    def terminate_pending_commands(self):
        for node in self.partial_program_order.get_executing_normal_and_spec_nodes():
            node.reset_to_ready()
            # proc.terminate()

def main():
    args = init()

    # Format logging
    # ref: https://docs.python.org/3/library/logging.html#formatter-objects
    if args.log_file is None:
        logging.basicConfig(format="%(levelname)s|%(asctime)s|%(message)s")
    else:
        logging.basicConfig(format="%(levelname)s|%(asctime)s|%(message)s",
                            filename=f"{os.path.abspath(args.log_file)}",
                            filemode="w")

    # Set debug level
    if args.debug_level == 1:
        logging.getLogger().setLevel(logging.INFO)
    elif args.debug_level >= 2:
        logging.getLogger().setLevel(logging.DEBUG)

    # Set optimization options
    config.SANDBOX_KILLING = args.sandbox_killing
    config.SPECULATE_IMMEDIATELY = args.speculate_immediately
    scheduler = Scheduler(config.SCHEDULER_SOCKET, args.window)
    scheduler.run()


if __name__ == "__main__":
    main()
