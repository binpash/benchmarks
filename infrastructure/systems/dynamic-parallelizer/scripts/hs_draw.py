import sys
from dataclasses import dataclass
from datetime import datetime

import matplotlib.pyplot as plt

filename = sys.argv[1]

@dataclass
class Entry:
    state: str
    cnid: str
    t: float

@dataclass
class Bar:
    start_s: str
    end_s: str
    start_t: float
    end_t: float

class BarList:
    def __init__(self, barlist):
        self.barlist = barlist

    def __lt__(self, other):
        if self.barlist[-1].end_s != 'COMMIT' or other.barlist[-1].end_s != 'COMMIT':
            return False
        return self.barlist[-1].end_t < other.barlist[-1].end_t

    def __le__(self, other):
        if self.barlist[-1].end_s != 'COMMIT' or other.barlist[-1].end_s != 'COMMIT':
            return True
        return self.barlist[-1].end_t < other.barlist[-1].end_t

def get_timestamp(l:str):
    _, s, _ = l.split('|', maxsplit=2)
    datetime_obj = datetime.strptime(s, "%Y-%m-%d %H:%M:%S,%f")
    timestamp = datetime_obj.timestamp()
    return float(timestamp)

def parse_body(l: str):
    index = l.find('[STATE_LOG] ') + len('[STATE_LOG] ')
    body = l.strip()[index:]
    cnid, state = body.split(':')
    return cnid, state.strip()

event_logs = []
with open(filename) as f:
    for l in f:
        if 'STATE_LOG' in l:
            l = l.strip()
            t = get_timestamp(l)
            cnid, state = parse_body(l)
            event_logs.append(Entry(state, cnid, t))
x0 = None
bars = {}
total_state = {}
for event in event_logs:
    state, cnid, t = event.state, event.cnid, event.t
    if x0 is None:
        x0 = t
    if not cnid in bars:
        if not state in ['READY', 'EXE', 'SPEC_E']:
            breakpoint()
        assert state in ['READY', 'EXE', 'SPEC_E']
        bars[cnid] = []
    if state in ['EXE', 'SPEC_E']:
        total_state[cnid] = event
    elif state == 'SPEC_F':
        prev_event = total_state[cnid]
        assert prev_event.state == 'SPEC_E'
        bars[cnid].append(Bar(prev_event.state, event.state, prev_event.t, event.t))
        total_state[cnid] = event
    elif state == 'COMMIT':
        prev_event = total_state[cnid]
        assert prev_event.state in ['SPEC_F', 'EXE']
        bars[cnid].append(Bar(prev_event.state, event.state, prev_event.t, event.t))
        total_state[cnid] = event
    elif state == 'READY':
        if cnid in total_state:
            prev_event = total_state[cnid]
            assert prev_event.state in ['SPEC_E', 'SPEC_F']
            bars[cnid].append(Bar(prev_event.state, event.state, prev_event.t, event.t))
        del total_state[cnid]
    else:
        assert False

keys = list(bars.keys())
keys = sorted(keys, key=lambda key: BarList(bars[key]))

fig, ax = plt.subplots(figsize=(24, 0.3*len(keys)))

d = {'EXE': (0, 0, 0), 'SPEC_E': (0.8, 0.8, 0.8), 'SPEC_F': (0.4, 0.4, 0.4)}
for i, key in enumerate(keys):
    y_pos = i
    print(f'key: {key}')
    for bar in bars[key]:
        color = d[bar.start_s]
        if bar.end_s == 'READY':
            color = (0.8, color[1]*0.5, color[2]*0.5)
        print((bar.start_t-x0, bar.end_t-x0))
        ax.broken_barh([(bar.start_t-x0, bar.end_t-bar.start_t)], (y_pos-0.4, 0.8), color=color)

print(len(keys))
ax.set_yticks(list(range(len(keys))))
ax.set_yticklabels(keys)
plt.tight_layout()
plt.savefig('hs_plot.pdf')
