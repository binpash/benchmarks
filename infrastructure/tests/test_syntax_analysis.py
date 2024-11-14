from collections import Counter
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
    count = Counter()
    count_nodes(parse_literal(code), count)
    return {str(n): c for n, c in count.items()}

class TestBasic(unittest.TestCase):
    def test_cmd(self):
        self.assertEqual(count("echo hi"), {'command(echo)': 1})
        self.assertEqual(count('printf "hello\n"'), {'command(printf)': 1, 'quoted_control': 1})
        self.assertEqual(count("./myprogram x y z"), {'command(./myprogram)': 1})
        self.assertEqual(count("./myprogram x y z &"), {'command(./myprogram)': 1,
                                                        'background': 1})
        self.assertEqual(count("! ./myprogram x y z"), {'command(./myprogram)': 1,
                                                        'negate_command': 1})

    def test_pipe(self):
        self.assertEqual(count("echo hi | echo bye"),
                         {'command(echo)': 2,
                          'pipeline': 1})
        self.assertEqual(count("echo hi | tr h b | echo"),
                         {'command(echo)': 2,
                         'command(tr)': 1,
                          'pipeline': 1})

    def test_seq(self):
        self.assertEqual(count("echo hi; echo bye"),
                         {'command(echo)': 2})
        self.assertEqual(count("echo hi && echo bye"),
                         {'command(echo)': 2,
                          'and_command': 1})
        self.assertEqual(count("echo hi || echo bye"),
                         {'command(echo)': 2,
                          'or_command': 1})
        self.assertEqual(count("echo hi && echo bye || exit 1"),
                         {'command(echo)': 2,
                         'command(exit)': 1,
                          'and_command': 1,
                          'or_command': 1})

    def test_ctrl(self):
        self.assertEqual(count("if echo hi; then echo yes; else echo no; fi"),
                         {'command(echo)': 3,
                          'if_command': 1})
        self.assertEqual(count("for x in a b c; do echo $x; done"),
                {'command(echo)': 1, 'for_command': 1, 'variable_use': 1})
        self.assertEqual(count("while ping -n 1 test.com; do echo $x; done"),
                         {'command(echo)': 1,
                         'command(ping)': 1,
                         'variable_use': 1,
                          'while_command': 1,})
        self.assertEqual(count("case $x in (hi) echo bye;; esac"),
                         {'command(echo)': 1,
                          'case_command': 1,
                          'variable_use': 1})

    def test_redir(self):
        self.assertEqual(count("echo hi > somefile.txt"),
                         {'command(echo)': 1,
                          'file_redirection': 1})


if __name__ == '__main__':
    unittest.main()
