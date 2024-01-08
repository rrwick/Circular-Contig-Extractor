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


def test_compression_type_1():
    filename = pathlib.Path(__file__).resolve().parent / 'compression' / 'uncompressed'
    compression_type = circular_contig_extractor.get_compression_type(filename)
    assert compression_type == 'plain'
    open_func = circular_contig_extractor.get_open_func(filename)
    assert open_func == open


def test_compression_type_2():
    filename = pathlib.Path(__file__).resolve().parent / 'compression' / 'gzipped'
    compression_type = circular_contig_extractor.get_compression_type(filename)
    assert compression_type == 'gz'
    open_func = circular_contig_extractor.get_open_func(filename)
    assert open_func == gzip.open


def test_compression_type_3():
    filename = pathlib.Path(__file__).resolve().parent / 'compression' / 'bzip2ed'
    with pytest.raises(SystemExit) as exit_message:
        circular_contig_extractor.get_compression_type(filename)
    assert 'cannot use bzip2' in str(exit_message.value)


def test_compression_type_4():
    filename = pathlib.Path(__file__).resolve().parent / 'compression' / 'zipped'
    with pytest.raises(SystemExit) as exit_message:
        circular_contig_extractor.get_compression_type(filename)
    assert 'cannot use zip' in str(exit_message.value)


def test_help_1():
    with pytest.raises(SystemExit) as sysexit:
        circular_contig_extractor.main(['--help'])
        assert sysexit.code == 0


def test_help_2():
    with pytest.raises(SystemExit) as sysexit:
        circular_contig_extractor.main([])
        assert sysexit.code == 1


def test_check_args():
    Args = collections.namedtuple('Args', ['in_gfa', 'min', 'max', 'query', 'mash'])
    circular_contig_extractor.check_args(Args(in_gfa='in', min=None, max=None, query=None, mash=0.1))
    circular_contig_extractor.check_args(Args(in_gfa='in', min=100, max=None, query=None, mash=0.1))
    circular_contig_extractor.check_args(Args(in_gfa='in', min=None, max=100, query=None, mash=0.1))
    circular_contig_extractor.check_args(Args(in_gfa='in', min=None, max=None, query=None, mash=0.9))
    with pytest.raises(SystemExit):
        circular_contig_extractor.check_args(Args(in_gfa='in', min=-100, max=None, query=None, mash=0.1))
    with pytest.raises(SystemExit):
        circular_contig_extractor.check_args(Args(in_gfa='in', min=None, max=-100, query=None, mash=0.1))
    with pytest.raises(SystemExit):
        circular_contig_extractor.check_args(Args(in_gfa='in', min=1000, max=100, query=None, mash=0.1))
    with pytest.raises(SystemExit):
        circular_contig_extractor.check_args(Args(in_gfa='in', min=None, max=None, query=None, mash=2.0))


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
