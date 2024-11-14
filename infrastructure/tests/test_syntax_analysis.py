import unittest
import tempfile
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from syntax_analysis import *

def parse_literal(code):
    with tempfile.NamedTemporaryFile() as fp:
        fp.write(code.encode())
        fp.seek(0)
        return parse_shell_script(fp.name)

def count(code):
    return count_nodes(parse_literal(code))

class TestBasic(unittest.TestCase):
    def test_cmd(self):
        self.assertEqual(count("echo hi"), {'CommandNode': 1})
        self.assertEqual(count('printf "hello\n"'), {'CommandNode': 1})
        self.assertEqual(count("./myprogram x y z"), {'CommandNode': 1})
        self.assertEqual(count("./myprogram x y z &"), {'CommandNode': 1,
                                                        'BackgroundNode': 1})
        self.assertEqual(count("! ./myprogram x y z"), {'CommandNode': 1,
                                                        'NotNode': 1})

    def test_pipe(self):
        self.assertEqual(count("echo hi | echo bye"),
                         {'CommandNode': 2,
                          'PipeNode': 1})
        self.assertEqual(count("echo hi | tr h b | echo"),
                         {'CommandNode': 3,
                          'PipeNode': 1})

    def test_seq(self):
        self.assertEqual(count("echo hi; echo bye"),
                         {'CommandNode': 2,
                          'SemiNode': 1})
        self.assertEqual(count("echo hi && echo bye"),
                         {'CommandNode': 2,
                          'AndNode': 1})
        self.assertEqual(count("echo hi || echo bye"),
                         {'CommandNode': 2,
                          'OrNode': 1})
        self.assertEqual(count("echo hi && echo bye || exit 1"),
                         {'CommandNode': 3,
                          'AndNode': 1,
                          'OrNode': 1})

    def test_ctrl(self):
        self.assertEqual(count("if echo hi; then echo yes; else echo no; fi"),
                         {'CommandNode': 3,
                          'IfNode': 1})
        self.assertEqual(count("for x in a b c; do echo $x; done"),
                         {'CommandNode': 1,
                          'ForNode': 1})
        self.assertEqual(count("while ping -n 1 test.com; do echo $x; done"),
                         {'CommandNode': 2,
                          'WhileNode': 1})
        self.assertEqual(count("case $x in (hi) echo bye;; esac"),
                         {'CommandNode': 1,
                          'CaseNode': 1})

    def test_redir(self):
        self.assertEqual(count("echo hi > somefile.txt"),
                         {'CommandNode': 1,
                          'FileRedirNode': 1})


if __name__ == '__main__':
    unittest.main()
