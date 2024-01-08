#!/usr/bin/env python3
"""
Copyright 2024 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Circular-Contig-Extractor

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not,
see <https://www.gnu.org/licenses/>.
"""

import argparse
import collections
import gzip
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

__version__ = '0.1.0'


def get_arguments(args):
    parser = MyParser(description='Circular Contig Extractor', add_help=False,
                      formatter_class=MyHelpFormatter)

    positional_args = parser.add_argument_group('Positional arguments')
    positional_args.add_argument('in_gfa', type=pathlib.Path,
                                 help='Input assembly GFA file')

    clustering_args = parser.add_argument_group('Contig size settings')
    clustering_args.add_argument('--min', type=int, default=None,
                                 help='Minimum acceptable contig size in bp (default: no minimum size)')
    clustering_args.add_argument('--max', type=int, default=None,
                                 help='Maximum acceptable contig size in bp (default: no maximum size)')

    clustering_args = parser.add_argument_group('Query settings')
    clustering_args.add_argument('--query', type=pathlib.Path, default=None,
                                 help='Query reference sequence(s) in FASTA format (default: none)')
    clustering_args.add_argument('--mash', type=float, default=0.1,
                                 help='Maximum acceptable Mash distance to query sequence')

    other_args = parser.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version',
                            version='Circular Contig Extractor v' + __version__,
                            help="Show program's version number and exit")

    args = parser.parse_args(args)
    return args


def main(args=None):
    args = get_arguments(args)
    check_args(args)
    contigs, links = load_gfa(args.in_gfa)
    contigs = find_circular_contigs(contigs, links)
    if not contigs:
        sys.exit()
    contigs = trim_overlaps(contigs)
    if args.min is not None or args.max is not None:
        contigs = filter_by_size(contigs, args.min, args.max)
        if not contigs:
            sys.exit()
    if args.query is not None:
        contigs = filter_by_query(contigs, args.query, args.mash)
        if not contigs:
            sys.exit()
    output_contigs(contigs)


def check_args(args):
    if args.min is not None and args.min <= 0:
        sys.exit('Error: --min must be greater than 0')
    if args.max is not None and args.max <= 0:
        sys.exit('Error: --max must be greater than 0')
    if args.min is not None and args.max is not None and args.max < args.min:
        sys.exit('Error: --max must be greater than or equal to --min')
    if args.mash < 0:
        sys.exit('Error: --mash must be greater than or equal to 0.0')
    if args.mash >= 1:
        sys.exit('Error: --mash must be less than 1.0')
    if args.query is not None:
        check_file_exists(args.query)


def load_gfa(filename):
    print(f'\nLoading {filename}:', file=sys.stderr)
    contigs, links = [], []
    with get_open_func(filename)(filename, 'rt') as gfa_file:
        for line in gfa_file:
            parts = line.rstrip('\n').split('\t')
            if parts[0] == 'S':
                contigs.append(parts[1:3])
            if parts[0] == 'L':
                links.append(parts[1:6])
    print(f'  {len(contigs)} contig{"" if len(contigs) == 1 else "s"}', file=sys.stderr)
    print(f'  {len(links)} link{"" if len(links) == 1 else "s"}', file=sys.stderr)
    return contigs, links


def find_circular_contigs(contigs, links):
    print(f'\nFinding circular contigs:', file=sys.stderr)
    circular_links = {}
    for seg_a, strand_a, seg_b, strand_b, cigar in links:
        if seg_a == seg_b and strand_a == strand_b:
            circular_links[seg_a] = cigar
    for seg_a, strand_a, seg_b, strand_b, _ in links:
        if seg_a != seg_b or strand_a != strand_b:
            circular_links.pop(seg_a, None)
            circular_links.pop(seg_b, None)
    circular_contigs = []
    for name, seq in contigs:
        if name in circular_links:
            circular_contigs.append((name, seq, circular_links[name]))
            print(f'  {name}: {len(seq):,} bp', file=sys.stderr)
    if not circular_contigs:
        print('  no circular contigs found\n', file=sys.stderr)
    return circular_contigs


def trim_overlaps(contigs):
    print(f'\nTrimming overlaps:', file=sys.stderr)
    trimmed_contigs = []
    for name, seq, cigar in contigs:
        print(f'  {name}: ', file=sys.stderr, end='')
        overlap = get_overlap_from_cigar(cigar)
        if overlap is None:
            print(f'cannot determine overlap', file=sys.stderr)
        elif overlap == 0:
            print(f'no overlap', file=sys.stderr)
        else:
            print(f'trimming {overlap} bp', file=sys.stderr)
        trimmed_contigs.append((name, trim_seq(seq, overlap)))
    return trimmed_contigs


def get_overlap_from_cigar(cigar):
    match = re.match(r'^(\d+)M$', cigar)
    if match:
        return int(match.group(1))
    else:
        return None


def trim_seq(seq, trim_amount):
    if trim_amount is None or trim_amount == 0:
        return seq
    else:
        return seq[:-trim_amount]


def filter_by_size(contigs, min_size, max_size):
    print(f'\nFiltering by size:', file=sys.stderr)
    if min_size is not None:
        contigs = [s for s in contigs if len(s[1]) >= min_size]
    if max_size is not None:
        contigs = [s for s in contigs if len(s[1]) <= max_size]
    for name, seq in contigs:
        print(f'  {name}: {len(seq):,} bp', file=sys.stderr)
    if not contigs:
        print('  no contigs satisfy size parameters\n', file=sys.stderr)
    return contigs


def filter_by_query(contigs, query_filename, mash_threshold):
    print(f'\nFiltering by query sequence(s):', file=sys.stderr)
    mash_distances = get_all_mash_distances(contigs, query_filename)
    closest_distances = {contig_name: sorted(distances)[0]
                         for contig_name, distances in mash_distances.items()}
    matching_contigs = []
    for name, seq in contigs:
        closest_distance, closest_query = closest_distances[name]
        if closest_distance <= mash_threshold:
            matching_contigs.append((name, seq))
            print(f'  {name}: {closest_distance:.5f} Mash distance to {closest_query}', file=sys.stderr)
    if not matching_contigs:
        print('  no contigs match query sequence(s)\n', file=sys.stderr)
    return matching_contigs


def get_all_mash_distances(contigs, query_filename):
    mash_distances = collections.defaultdict(list)
    with tempfile.TemporaryDirectory() as temp_dir:
        query_fasta = pathlib.Path(temp_dir) / 'query.fasta'
        contig_fasta = pathlib.Path(temp_dir) / 'contig.fasta'
        for query_name, query_seq in iterate_fasta(query_filename):
            write_fasta(query_name, query_seq, query_fasta)
            for contig_name, contig_seq in contigs:
                write_fasta(contig_name, contig_seq, contig_fasta)
                dist = get_mash_distance(query_fasta, contig_fasta)
                mash_distances[contig_name].append((dist, query_name))
    return mash_distances


def write_fasta(name, seq, filename):
    with open(filename, 'wt') as f:
        f.write(f'>{name}\n')
        f.write(f'{seq}\n')


def get_mash_distance(fasta_a, fasta_b):
    mash_command = ['mash', 'dist', fasta_a, fasta_b]
    distances = []
    p = subprocess.Popen(mash_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stderr = p.stderr.read().strip()
    for line in p.stdout:
        parts = line.split('\t')
        distance = float(parts[2])
        distances.append(distance)
    if len(distances) != 1:
        sys.exit(f'Error: mash dist did not run successfully:\n{stderr}')
    assert len(distances) == 1
    return distances[0]


def output_contigs(contigs):
    print(f'\nOutputting contigs:', file=sys.stderr)
    for name, seq in contigs:
        print(f'  {name}: {len(seq):,} bp', file=sys.stderr)
        print(f'>{name}')
        print(seq)
    print('', file=sys.stderr)


def iterate_fasta(filename):
    """
    Takes a FASTA file as input and yields the contents as (name, seq) tuples.
    """
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    name_parts = name.split(maxsplit=1)
                    contig_name = name_parts[0]
                    yield contig_name, ''.join(sequence)
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line.upper())
        if name:
            name_parts = name.split(maxsplit=1)
            contig_name = name_parts[0]
            yield contig_name, ''.join(sequence)


def check_file_exists(filename):
    if filename.is_dir():
        sys.exit(f'Error: {filename} is a directory, not a file')
    if not filename.is_file():
        sys.exit(f'Error: {filename} does not exist')


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(str(filename), 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
DIM = '\033[2m'


class MyParser(argparse.ArgumentParser):
    """
    This subclass of ArgumentParser changes the error messages, such that if the script is run with
    no other arguments, it will display the help text. If there is a different error, it will give
    the normal response (usage and error).
    """
    def error(self, message):
        if len(sys.argv) == 1:  # if no arguments were given.
            self.print_help(file=sys.stderr)
            sys.exit(1)
        else:
            super().error(message)


class MyHelpFormatter(argparse.HelpFormatter):

    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ['COLUMNS'] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 3), 40)
        self.colours = get_colours_from_tput()
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        """
        Override this function to add default values, but only when 'default' is not already in the
        help text.
        """
        help_text = action.help
        if action.default != argparse.SUPPRESS and action.default is not None:
            if 'default: DEFAULT' in help_text:
                help_text = help_text.replace('default: DEFAULT', f'default: {action.default}')
        return help_text

    def start_section(self, heading):
        """
        Override this method to add bold underlining to section headers.
        """
        if self.colours > 1:
            heading = BOLD + heading + END_FORMATTING
        super().start_section(heading)

    def _format_action(self, action):
        """
        Override this method to make help descriptions dim.
        """
        help_position = min(self._action_max_length + 2, self._max_help_position)
        help_width = self._width - help_position
        action_width = help_position - self._current_indent - 2
        action_header = self._format_action_invocation(action)
        if not action.help:
            tup = self._current_indent, '', action_header
            action_header = '%*s%s\n' % tup
            indent_first = 0
        elif len(action_header) <= action_width:
            tup = self._current_indent, '', action_width, action_header
            action_header = '%*s%-*s  ' % tup
            indent_first = 0
        else:
            tup = self._current_indent, '', action_header
            action_header = '%*s%s\n' % tup
            indent_first = help_position
        parts = [action_header]
        if action.help:
            help_text = self._expand_help(action)
            help_lines = self._split_lines(help_text, help_width)
            first_line = help_lines[0]
            if self.colours > 8:
                first_line = DIM + first_line + END_FORMATTING
            parts.append('%*s%s\n' % (indent_first, '', first_line))
            for line in help_lines[1:]:
                if self.colours > 8:
                    line = DIM + line + END_FORMATTING
                parts.append('%*s%s\n' % (help_position, '', line))
        elif not action_header.endswith('\n'):
            parts.append('\n')
        for subaction in self._iter_indented_subactions(action):
            parts.append(self._format_action(subaction))
        return self._join_parts(parts)


def get_colours_from_tput():
    try:
        return int(subprocess.check_output(['tput', 'colors']).decode().strip())
    except (ValueError, subprocess.CalledProcessError, FileNotFoundError, AttributeError):
        return 1


if __name__ == '__main__':
    main()
