#!/usr/bin/env python3

# Good points of reference:
# https://github.com/binpash/shasta/blob/main/shasta/ast_node.py
# https://github.com/binpash/Shseer/blob/8bb9e72f7fe1b4703fc963bfa5d5bd2837e80ab3/src/shseer/symb.py

import sys
import libdash
import functools as ft
from shasta.json_to_ast import to_ast_node
from shasta.ast_node import (
    AndNode,
    ArgChar,
    EArgChar,
    CArgChar,
    TArgChar,
    AArgChar, VArgChar, QArgChar, BArgChar,
    AstNode,
    CArgChar,
    CaseNode,
    CommandNode,
    DefunNode,
    DupRedirNode,
    FileRedirNode,
    ForNode,
    HeredocRedirNode,
    IfNode,
    NotNode,
    OrNode,
    PipeNode,
    RedirNode,
    SemiNode,
    SubshellNode,
    WhileNode,
    BackgroundNode,
    AssignNode,
)
from collections import namedtuple

sys.setrecursionlimit (9001)

# Node With Source
NodeWSource = namedtuple('NodeWSource', ['node', 'source_syntax', 'linum_before', 'linum'])

first_time = True
def parse_shell_script(path):
    global first_time
    raw_asts = libdash.parse(path, first_time)
    first_time = False
    return [NodeWSource(to_ast_node(raw_ast),
                        source,
                        linum_before,
                        linum)
               for (raw_ast, source, linum_before, linum) in raw_asts]

def is_pipeline(ast):
    return isinstance(ast.node, PipeNode)

def argchar_to_subnodes(argchar):
    match argchar:
        case CArgChar() | EArgChar() | TArgChar():
            return []
        case AArgChar(arg=more) | VArgChar(arg=more) | QArgChar(arg=more):
            return flatmap(argchar_to_subnodes, more)
        case BArgChar(node=node):
            return [node]

def flatmap(f, l):
    res = []
    for rl in map(f, l):
        res += rl
    return res

def merge_counts(d1, d2):
    h = d1.copy()
    for k, v in d2.items():
        if k in h:
            h[k] += v
        else:
            h[k] = v
    return h

def add_ast(d, a):
    name = type(a).__name__
    if name in d:
        d[name] += 1
    else:
        d[name] = 1
    return d

def count_nodes(asts):
    match asts:
        case NodeWSource(node=subnode):
            return count_nodes(subnode)
        case PipeNode(items=subnodes):
            return add_ast(count_nodes(subnodes), asts)
        case SubshellNode(body=subnode) \
             | NotNode(body=subnode) \
             | RedirNode(node=subnode) \
             | BackgroundNode(node=subnode) \
             | DefunNode(body=subnode):
            return add_ast(count_nodes(subnode), asts)
        case AndNode(left_operand=l, right_operand=r) \
             | OrNode(left_operand=l, right_operand=r) \
             | SemiNode(left_operand=l, right_operand=r) \
             | WhileNode(test=l, body=r):
            return add_ast(count_nodes([l, r]), asts)
        case IfNode(cond=t, then_b=thn, else_b=els):
            return add_ast(count_nodes([t, thn, els]), asts)
        case ForNode(argument=lolo_argchar, body=subnode):
            return add_ast(count_nodes([subnode] + \
                                       flatmap(lambda lo_argchar: flatmap(argchar_to_subnodes, lo_argchar),
                                               lolo_argchar)),
                           asts)
        case CaseNode(argument=argchars, cases=cases):
            return add_ast(count_nodes([c['cbody'] for c in cases] + \
                                       flatmap(argchar_to_subnodes, argchars)), asts)
        case AssignNode(val=argchars) \
             | FileRedirNode(arg=argchars) \
             | DupRedirNode(arg=argchars) \
             | HeredocRedirNode(arg=argchars):
            return add_ast(count_nodes(flatmap(argchar_to_subnodes, argchars)), asts)
        case CommandNode(arguments=lolo_argchar, redir_list=redir_list):
            return add_ast(count_nodes(flatmap(lambda lo_argchar: flatmap(argchar_to_subnodes, lo_argchar),
                                               lolo_argchar) + \
                                       redir_list),
                           asts)
        case [*subnodes]:
            return ft.reduce(merge_counts, map(count_nodes, subnodes), dict())
        case other:
            raise Exception(f"oops: {other} of type {type(other)}")

if __name__ == '__main__':
    p = "script.sh"
    asts = parse_shell_script(p)
    print(count_nodes(asts))
