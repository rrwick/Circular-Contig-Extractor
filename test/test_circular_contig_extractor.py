"""
This module contains some tests for Circular Contig Extractor. To run them, execute
`python3 -m pytest` from the root Circular Contig Extractor directory.

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

import circular_contig_extractor
import collections
import glob
import gzip
import os
import pathlib
import pytest
import tempfile


def file_dir():
    return pathlib.Path(__file__).resolve().parent / 'files'


def test_compression_type_1():
    filename = file_dir() / 'uncompressed'
    compression_type = circular_contig_extractor.get_compression_type(filename)
    assert compression_type == 'plain'
    open_func = circular_contig_extractor.get_open_func(filename)
    assert open_func == open


def test_compression_type_2():
    filename = file_dir() / 'gzipped'
    compression_type = circular_contig_extractor.get_compression_type(filename)
    assert compression_type == 'gz'
    open_func = circular_contig_extractor.get_open_func(filename)
    assert open_func == gzip.open


def test_compression_type_3():
    filename = file_dir() / 'bzip2ed'
    with pytest.raises(SystemExit) as exit_message:
        circular_contig_extractor.get_compression_type(filename)
    assert 'cannot use bzip2' in str(exit_message.value)


def test_compression_type_4():
    filename = file_dir() / 'zipped'
    with pytest.raises(SystemExit) as exit_message:
        circular_contig_extractor.get_compression_type(filename)
    assert 'cannot use zip' in str(exit_message.value)


def test_help_1():
    with pytest.raises(SystemExit) as sysexit:
        circular_contig_extractor.main(['--help'])


def test_help_2():
    with pytest.raises(SystemExit) as sysexit:
        circular_contig_extractor.main([])


def test_check_args():
    Args = collections.namedtuple('Args', ['in_gfa', 'min', 'max', 'query', 'mash'])
    check_args = circular_contig_extractor.check_args
    check_args(Args(in_gfa='in', min=None, max=None, query=None, mash=0.1))
    check_args(Args(in_gfa='in', min=100, max=None, query=None, mash=0.1))
    check_args(Args(in_gfa='in', min=None, max=100, query=None, mash=0.1))
    check_args(Args(in_gfa='in', min=None, max=None, query=None, mash=0.9))
    with pytest.raises(SystemExit):
        check_args(Args(in_gfa='in', min=-100, max=None, query=None, mash=0.1))
    with pytest.raises(SystemExit):
        check_args(Args(in_gfa='in', min=None, max=-100, query=None, mash=0.1))
    with pytest.raises(SystemExit):
        check_args(Args(in_gfa='in', min=1000, max=100, query=None, mash=0.1))
    with pytest.raises(SystemExit):
        check_args(Args(in_gfa='in', min=None, max=None, query=None, mash=2.0))
    with pytest.raises(SystemExit):
        check_args(Args(in_gfa='in', min=None, max=None, query=None, mash=-0.1))
    with pytest.raises(SystemExit):
        check_args(Args(in_gfa='in', min=None, max=None, query=pathlib.Path('bad'), mash=0.5))


def test_check_file_exists():
    good = file_dir() / 'uncompressed'
    bad = file_dir() / 'not_a_file'
    directory = pathlib.Path(__file__).resolve().parent
    circular_contig_extractor.check_file_exists(good)
    with pytest.raises(SystemExit):
        circular_contig_extractor.check_file_exists(bad)
    with pytest.raises(SystemExit):
        circular_contig_extractor.check_file_exists(directory)


def test_load_gfa():
    filename = file_dir() / 'graph_1.gfa'
    contigs, links = circular_contig_extractor.load_gfa(filename)
    assert len(contigs) == 3
    assert len(links) == 3


def test_find_circular_contigs_1():
    filename = file_dir() / 'graph_1.gfa'
    contigs, links = circular_contig_extractor.load_gfa(filename)
    circular_contigs = circular_contig_extractor.find_circular_contigs(contigs, links)
    assert len(circular_contigs) == 1
    assert circular_contigs[0][0] == '2'


def test_find_circular_contigs_2():
    filename = file_dir() / 'graph_2.gfa'
    contigs, links = circular_contig_extractor.load_gfa(filename)
    circular_contigs = circular_contig_extractor.find_circular_contigs(contigs, links)
    assert len(circular_contigs) == 0


def test_trim_overlaps():
    contigs = [('1', 'ACGATCAGCACT', '0M'),
               ('2', 'ACGATCAGCACT', '5M'),
               ('3', 'ACGATCAGCACT', '*')]
    trimmed_contigs = circular_contig_extractor.trim_overlaps(contigs)
    assert trimmed_contigs == [('1', 'ACGATCAGCACT'), ('2', 'ACGATCA'), ('3', 'ACGATCAGCACT')]

def test_get_overlap_from_cigar():
    assert circular_contig_extractor.get_overlap_from_cigar('0M') == 0
    assert circular_contig_extractor.get_overlap_from_cigar('123M') == 123
    assert circular_contig_extractor.get_overlap_from_cigar('abc') is None
    assert circular_contig_extractor.get_overlap_from_cigar('2M1D4M') is None
    assert circular_contig_extractor.get_overlap_from_cigar('') is None


def test_trim_seq():
    assert circular_contig_extractor.trim_seq('ACACGACTACG', None) == 'ACACGACTACG'
    assert circular_contig_extractor.trim_seq('ACACGACTACG', 0) == 'ACACGACTACG'
    assert circular_contig_extractor.trim_seq('ACACGACTACG', 1) == 'ACACGACTAC'
    assert circular_contig_extractor.trim_seq('ACACGACTACG', 2) == 'ACACGACTA'
    assert circular_contig_extractor.trim_seq('ACACGACTACG', 3) == 'ACACGACT'
    assert circular_contig_extractor.trim_seq('ACACGACTACG', 4) == 'ACACGAC'
    assert circular_contig_extractor.trim_seq('ACACGACTACG', 5) == 'ACACGA'


def test_filter_by_size():
    contigs = [('1', 'ACGATC'), ('2', 'ACGATCAGC'), ('3', 'ACGATCAGCACT')]
    filtered_contigs = circular_contig_extractor.filter_by_size(contigs, None, None)
    assert filtered_contigs == [('1', 'ACGATC'), ('2', 'ACGATCAGC'), ('3', 'ACGATCAGCACT')]
    filtered_contigs = circular_contig_extractor.filter_by_size(contigs, 0, 100)
    assert filtered_contigs == [('1', 'ACGATC'), ('2', 'ACGATCAGC'), ('3', 'ACGATCAGCACT')]
    filtered_contigs = circular_contig_extractor.filter_by_size(contigs, 8, 10)
    assert filtered_contigs == [('2', 'ACGATCAGC')]
    filtered_contigs = circular_contig_extractor.filter_by_size(contigs, 9, 9)
    assert filtered_contigs == [('2', 'ACGATCAGC')]
    filtered_contigs = circular_contig_extractor.filter_by_size(contigs, None, 9)
    assert filtered_contigs == [('1', 'ACGATC'), ('2', 'ACGATCAGC')]
    filtered_contigs = circular_contig_extractor.filter_by_size(contigs, 9, None)
    assert filtered_contigs == [('2', 'ACGATCAGC'), ('3', 'ACGATCAGCACT')]
    filtered_contigs = circular_contig_extractor.filter_by_size(contigs, 100, 200)
    assert filtered_contigs == []


def test_filter_by_query_1():
    graph_filename = file_dir() / 'graph_3.gfa'
    query_filename = file_dir() / 'query.fasta'
    contigs, _ = circular_contig_extractor.load_gfa(graph_filename)
    matching_contigs = circular_contig_extractor.filter_by_query(contigs, query_filename, 0.1)
    assert len(matching_contigs) == 1
    assert matching_contigs[0][0] == '2'


def test_filter_by_query_2():
    graph_filename = file_dir() / 'graph_3.gfa'
    query_filename = file_dir() / 'query.fasta'
    contigs, _ = circular_contig_extractor.load_gfa(graph_filename)
    matching_contigs = circular_contig_extractor.filter_by_query(contigs, query_filename, 0.00001)
    assert len(matching_contigs) == 0


def test_get_mash_distance():
    fasta = file_dir() / 'query.fasta'
    empty = file_dir() / 'empty_file'
    assert circular_contig_extractor.get_mash_distance(fasta, fasta) == 0.0
    with pytest.raises(SystemExit):
        assert circular_contig_extractor.get_mash_distance(empty, empty)


def test_write_fasta():
    with tempfile.TemporaryDirectory() as temp_dir:
        fasta = pathlib.Path(temp_dir) / 'temp.fasta'
        circular_contig_extractor.write_fasta('a', 'ACGACTACGATC', fasta)
        result = list(circular_contig_extractor.iterate_fasta(fasta))
        assert result == [('a', 'ACGACTACGATC')]
