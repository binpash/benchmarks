#!/usr/bin/env python3

import pandas as pd
import fnmatch
import viz.syntax as stx
import viz.dynamic as dyn
import sys
import ast

from project_root import get_project_root

root = get_project_root()
data_path = root / 'infrastructure/target/dynamic_analysis.jsonl'
input_size_path = root / 'infrastructure/data/size_inputs.jsonl'
loc_data_path = root / 'infrastructure/target/lines_of_code.csv'
syscall_data_path = root / 'infrastructure/data/no_of_syscalls.csv'

def read_sys_results():
    df = pd.read_csv(syscall_data_path)
    df.rename(columns={'Benchmark-(small)': 'benchmark', 'System Calls': 'sys_calls'}, inplace=True)
    return df

benchmark_category_style = {
    'analytics': ('System admin.', 'Data analysis', '\\cite{dgsh:ieee:2017,posh:atc:2020,drake2014command}'),
    'bio': ('Data analysis', 'Biology', '\\cite{Cappellini2019,puritz2019bio594,ibrahim2021tera}'),
    'ci-cd': ('Continuous Integration', 'Build scripts', '\\cite{riker2022,makself}'),
    'covid': ('Data analysis', 'Data extraction', '\\cite{covid-mts-source}'),
    'file-mod': ('Automation Now', 'Misc.', '\\cite{cito2020empirical,dgsh:ieee:2017,posh:atc:2020}'),
    'inference': ('Machine learning', 'Data analysis', '\\cite{pagedout2025issue6,tunney2023bash}'),
    'ml': ('Machine learning', 'Data analysis', '\\cite{scikit-learn}'),
    'nlp': ('Machine learning', 'Text processing', '\\cite{unix-for-poets-church}'),
    'oneliners': ('Automation Now', 'Text processing', ''),
    'pkg': ('Continuous Integration', 'Automation Now', '\\cite{pacaur,vasilakis2021preventing}'),
    'repl': ('System admin.', 'Misc.', '\\cite{posh:atc:2020,vpsaudit}'),
    'unixfun': ('Misc.', 'Text processing', '\\cite{bhandari2020solutions}'),
    'weather': ('Data analysis', 'Data extraction', '\\cite{hadoop-guide-2009}'),
    'web-search': ('Misc.', 'Text processing', '\\cite{pash2021}'),
}

script_to_citation = {
        'oneliners/scripts/spell.sh': '\\cite{bentley-pearl-cacm-1985}',
        'oneliners/scripts/uniq-ips.sh': '\\cite{majkowski2020bloom}',
        'oneliners/scripts/top-n.sh': '\\cite{bentley-pearl-cacm-1986}',
        'oneliners...': '\\cite{unix-cacm-1974, wicked-cool-shell-scripts}'
}

def short_category(benchmark):
    dom, style, _ = benchmark_category_style[benchmark]
    def shorten(str):
        return ''.join([x[0].upper() for x in str.split(' ')])
    return shorten(dom) + '/' + shorten(style)

def citation(benchmark):
    citation = benchmark_category_style[benchmark][2]
    if 'cite' in citation: 
        return citation
    return ''

benchmark_input_description = {
    'analytics': 'log files',
    'bio': '\\xxx',
    'ci-cd': 'application sources',
    'covid': 'transit data',
    'file-mod': 'media files',
    'inference': '\\xxx',
    'ml': '\\xxx',
    'nlp': 'books',
    'oneliners': '\\xxx',
    'pkg': 'packages',
    'repl': None,
    'unixfun': 'challenge inputs',
    'weather': 'temperature data',
    'web-search': 'root webpages',
}

scripts_to_include = [
    'covid-mts/scripts/1.sh',
    'file-enc/scripts/encrypt_files.sh',
    'log-analysis/scripts/nginx.sh',
    'media-conv/scripts/img_convert.sh',
    'nlp/scripts/bigrams.sh',
    'oneliners/scripts/spell.sh',
    'oneliners/scripts/uniq-ips.sh',
    'oneliners/scripts/top-n.sh',
    'unix50/scripts/1.sh',
    'riker/scripts/redis/build.sh',
    # max-temp is just 1
    # sklearn is just 1
    # aurpkg is just 1
    # bio is just 1
    'web-index/scripts/ngrams.sh',
    'vps-audit/scripts/vps-audit.sh',
    'makeself/makeself/test/lsmtest/lsmtest.sh'
]

def script_name(script):
    script_name_map = {
            "encrypt_files.sh": "encrypt.sh",
            "img_convert.sh": "img-conv.sh",
    }
    if script in script_name_map:
        return script_name_map[script]
    else:
        return script

def script_citation(script):
    if script in script_to_citation:
        return script_to_citation[script]
    return ''


def count_unique_cmds(series):
    return len({node for node in series if 'command(' in node})

def format_number(n):
    if n == 0:
        return '\\textasciitilde 0'
    return f"{n:.1f}"

def count_constructs(series):
    return len(set(series))

def read_loc_data():
    loc_data = pd.read_csv(loc_data_path, header=None)
    loc_data.columns = ['script', 'loc']
    map_df = stx.get_map_df()
    loc_data = loc_data.merge(map_df, on='script')
    loc_data_bench = loc_data.groupby('benchmark').agg({
        'loc': 'sum',
        'script': 'count'
    }).reset_index()
    loc_data_bench.rename(columns={'script': 'number_of_scripts'}, inplace=True)
    return loc_data, loc_data_bench

def prettify_bytes_number(n):
    if n < 1024:
        value, unit = n, "B"
    elif n < 1024 * 1024:
        value, unit = n / 1024, "KB"
    elif n < 1024 * 1024 * 1024:
        value, unit = n / (1024 * 1024), "MB"
    else:
        value, unit = n / (1024 * 1024 * 1024), "GB"
    
    if value < 10:
        decimals = 2
    elif value < 100:
        decimals = 1
    else:
        decimals = 0
    
    color = 'black' if unit == 'GB' else 'gray'
    return f"{value:.{decimals}f}\\textcolor{{{color}}}{{{unit}}}"

def prettify_big_count(n):
    if type(n) is not int:
        return n
    if n < 1000:
        return str(n)
    elif n >= 1000 and n < 10**6:
        return f"{(n/1000):.1f}k"
    elif n >= 10**6 and n < 10**9:
        return f"{n/(10**6):.1f}m"
    else:
        return f"{n/(10**9):.1f}g"

def make_input_description(row):
    if row['input_description']:
        desc = prettify_bytes_number(row['input_size']) # + ' of ' + row['input_description']
        # return f"\\multirow{{2}}{{*}}{{\\parbox{{\\idw}}{{{desc}}}}}"
        return desc
    else:
        # Center the N/A
        return "\\multicolumn{1}{c}{N/A}"

def main():
    syntax_script, syntax_bench = stx.read_data(True)

    syntax_script_all_cmds, syntax_bench_all_cmds = stx.read_data(False)
    dyn_script,    dyn_bench = dyn.read_data()
    loc_data_script, loc_data_bench = read_loc_data()
    sys_results = read_sys_results()

    # fill any benchmarks missing from sys_results
    for benchmark in syntax_bench['benchmark']:
        if benchmark not in sys_results['benchmark'].values:
            new_row = pd.DataFrame([{'benchmark': benchmark, 'sys_calls': '\\xxx'}])
            sys_results = pd.concat([sys_results, new_row], ignore_index=True)
    sys_results.reset_index(drop=True, inplace=True)
    # replace sys_results file_descriptors numbers with those from children_num_fds in dyn_bench
    sys_results = sys_results.merge(dyn_bench[['benchmark', 'children_num_fds']], on='benchmark')
    sys_results['file_descriptors'] = sys_results['children_num_fds']

    syntax_script_all_cmds['unique_cmds'] = syntax_script_all_cmds['nodes'].apply(count_unique_cmds)
    syntax_bench_all_cmds['unique_cmds'] = syntax_bench_all_cmds['nodes'].apply(count_unique_cmds)
    syntax_script['constructs'] = syntax_script['nodes'].apply(count_constructs)
    syntax_bench['constructs'] = syntax_bench['nodes'].apply(count_constructs)
    
    # all_scripts = set(syntax_script['script'].unique())
    
    # missing_in_dyn = all_scripts - set(dyn_script['script'].unique())
    # missing_in_loc_data = all_scripts - set(loc_data_script['script'].unique())
    # missing_in_cmds = all_scripts - set(syntax_script_all_cmds['script'].unique())
    
    # print("Missing in dyn_script:", missing_in_dyn)
    # print("Missing in loc_data_script:", missing_in_loc_data)
    # print("Missing in syntax_script_all_cmds:", missing_in_cmds)

    dyn_bench['input_description'] = dyn_bench['benchmark'].apply(lambda x: benchmark_input_description[x])

    big_bench = syntax_bench.merge(dyn_bench, on='benchmark')\
        .merge(loc_data_bench, on='benchmark')\
        .merge(syntax_bench_all_cmds[['benchmark', 'unique_cmds']], on='benchmark')\
        .merge(sys_results, on='benchmark')
    
    big_script = syntax_script.merge(dyn_script, on='script')\
        .merge(loc_data_script, on='script')\
        .merge(syntax_script_all_cmds[['script', 'unique_cmds']], on='script')

    embedding_df = pd.read_csv(root / 'infrastructure/data/embeddings.csv')
    embedding_df['embedding'] = embedding_df['embedding'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
    # Embedding is a list of numbers, turn them into columns
    embedding_df = pd.concat([embedding_df['benchmark'], embedding_df['embedding'].apply(pd.Series)], axis=1)

    # Calculate summary statistics
    agg_order = ['min', 'max', 'mean']
    summary_names = [s.capitalize() for s in agg_order]
    summary_stats = big_bench[['loc', 'constructs', 'unique_cmds', 'number_of_scripts']].agg(agg_order).reset_index()
    summary_stats.rename(columns={
        'index': 'benchmark',
        'loc': 'loc',
        'constructs': 'constructs',
        'unique_cmds': 'unique_cmds',
        'number_of_scripts': 'number_of_scripts'
    }, inplace=True)

    # Add placeholder values for non-numeric columns
    summary_stats['benchmark'] = summary_names
    summary_stats['sys_calls'] = big_bench['sys_calls'].agg(agg_order).values if big_bench['sys_calls'].dtype == 'int64' else '\\xxx'
    summary_stats['file_descriptors'] = big_bench['file_descriptors'].agg(agg_order).values if big_bench['file_descriptors'].dtype == 'int64' else '\\xxx'
    summary_stats['input_description'] = None
    summary_stats['time_in_shell'] = big_bench['time_in_shell'].agg(agg_order).values
    summary_stats['time_in_commands'] = big_bench['time_in_commands'].agg(agg_order).values
    summary_stats['max_unique_set_size'] = big_bench['max_unique_set_size'].agg(agg_order).values
    summary_stats['io_chars'] = big_bench['io_chars'].agg(agg_order).values
    summary_stats['number_of_scripts'] = big_bench['number_of_scripts'].agg(agg_order).values
    # Ignore rows that have no input_description
    summary_stats['input_size'] = big_bench[big_bench['input_description'].notnull()]['input_size'].agg(agg_order).values
    summary_stats['input_description'] = '\\xxx' # Have something so that N/A doesn't show up

    # Append summary rows to big_bench
    # big_bench = pd.concat([big_bench, summary_stats], ignore_index=True)
        

    print("""
    \\begin{tabular}{@{}llrrllrrrrrrlrl@{}}
    \\toprule
    \\multirow{2}{*}{Benchmark/Script} & \\multicolumn{3}{c}{Surface} & \\multicolumn{2}{c}{Inputs} & \\multicolumn{2}{c}{Syntax} & \\multicolumn{4}{c}{Dynamic} & \\multicolumn{2}{c}{System} & \\multirow{2}{*}{Source} \\\\
        \\cline{2-4} \\cline{5-6} \\cline{7-8} \\cline{9-12} \\cline{13-14}
                                      & Dom     & \\#.sh     & LoC     & Small & Large & \\#Cons & \\#Cmd & $t_{S}$  & $t_{C}$  & Mem   & I/O & \\#SC & \\#FD &   \\\\
        \\midrule
    """)
    # generate a big latex table with the following columns:
    # benchmark, short_category, number of scripts, LOC, number of constructs, number of unique commands, input description, time in shell, time in commands, max memory, IO
    for _, row in big_bench.iterrows():
        numscripts_shown = 0
        numscripts = row['number_of_scripts']
        print(f"\\bs{{{row['benchmark']}}} & {short_category(row['benchmark'])} & {row['number_of_scripts']} & {row['loc']} & \\xxx & {make_input_description(row)} & {row['constructs']} & {row['unique_cmds']} & {format_number(row['time_in_shell'])} & {format_number(row['time_in_commands'])} & {prettify_bytes_number(row['max_unique_set_size'])} & {prettify_bytes_number(row['io_chars'])} & {prettify_big_count(row['sys_calls'])} & {row['file_descriptors']} & {citation(row['benchmark'])} \\\\")
        # now print the details of all scripts in the benchmark
        for _, row_script in big_script.iterrows():
            if row_script['benchmark'] == row['benchmark'] and any([fnmatch.fnmatch(row_script['script'], pattern) for pattern in scripts_to_include]):
                # all columns except leave blank benchmark, category, number of scripts, input description
                print(f"\\hspace{{0.5em}} \\ttt{{{script_name(row_script['script'].split('/')[-1])}}} & & & {row_script['loc']} & & & {row_script['constructs']} & {row_script['unique_cmds']} & {row_script['time_in_shell']:.1f} & {row_script['time_in_commands']:.1f} & {prettify_bytes_number(row_script['max_unique_set_size'])} & {prettify_bytes_number(row_script['io_chars'])} & & & {script_citation(row_script['script'])} \\\\")
                numscripts_shown += 1
        if numscripts_shown < numscripts and numscripts > 1:
            print(f"\\hspace{{0.5em}} \\ldots & & & & & & & & & & & & & & {script_citation(row['benchmark'] + '...')} \\\\")

    print("\\midrule")
    for aggregate in ['Min', 'Mean', 'Max']:
        row = summary_stats[summary_stats['benchmark'] == aggregate].iloc[0]
       
        def format_value(value):
            def round_whole(numstr):
                if numstr.endswith(".0"):
                    return numstr[:-2]
                return numstr
            if isinstance(value, (int, float)):
                return round_whole(f"{value:.1f}") if isinstance(value, float) else f"{int(value)}"
            return value  # For non-numeric values

        print(f"{{\\textbf{{\\centering {row['benchmark']}}}}} & & {format_value(row['number_of_scripts'])} & {format_value(row['loc'])} & \\xxx & \\xxx & {format_value(row['constructs'])} & {format_value(row['unique_cmds'])} & {format_value(row['time_in_shell'])} & {format_value(row['time_in_commands'])} & {prettify_bytes_number(row['max_unique_set_size'])} & {prettify_bytes_number(row['io_chars'])} & {prettify_big_count(row['sys_calls'])} & {format_value(row['file_descriptors'])} & \\\\")

    print("""
    \\bottomrule
  \\end{tabular}
""")
    
if __name__ == '__main__':
    main()
