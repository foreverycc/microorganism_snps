#!/usr/bin/env python

"""Tests for `microorganism_snps` package."""

import sys
import os
import pytest
from argparse import Namespace
from unittest.mock import patch

# Add the src directory to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from microorganism_snps import microorganism_snps

def test_main():
    # Create a Namespace object with the test arguments
    test_args = Namespace(
        wkdir="./",
        outputbase="test_output",
        refGenome="NC_006273.2",
        inputFasta="./tests/cmv_50_strains.fasta",
        date="2024-10-16",
        minFragSize=100,
        minDepth=10,
        minFreq=0.1,
        segmentSize=1000
    )

    with patch('argparse.ArgumentParser.parse_args', return_value=test_args):
        try:
            # Call the main function
            result = microorganism_snps.main()
            # Add assertions here to check the result
            assert result is None  # or whatever you expect the result to be
        except Exception as e:
            print("An error occurred:")
            import traceback
            print(traceback.format_exc())
            raise

# You can add more test functions for other parts of your code
