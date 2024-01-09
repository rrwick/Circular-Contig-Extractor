# Circular Contig Extractor tests

`circular_contig_extractor.py` comes with a few automated tests to help with development and spotting bugs. You'll need [pytest](https://docs.pytest.org/en/latest/) installed to run them.

To run the tests, execute this command from the repo's root directory:
```bash
python3 -m pytest
```

Or if you have [coverage.py](https://coverage.readthedocs.io) installed, you can use it to get some code coverage stats:
```bash
coverage run --source . -m pytest && coverage report -m
```
