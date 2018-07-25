#!/usr/bin/env python
import os
import numpy
import pyferret
import argparse


def pyferret_runner(pyferret_file):
    """
    This is a wrapper around pyferret to run a jnl file one line at a time, should help with debugging. Will do very little with the actual file other than trying to run it line by line.

    failures are raised, normally if this happens we should go back and debug the jnl file.
    """
    if os.path.exists(pyferret_file):
        pyferret_file = list(open(pyferret_file))
    else:
        raise Exception("fatal, pyferret command file {} does not exist!".format(pyferret_file))
    
    pyferret.init(enterferret=False)
    for line in pyferret_file:
        err_int, err_msg = pyferret.run(command="line")
        if not err_int == pyferret.FERR_OK:
            raise Exception(err_msg)
    pyferret.stop()
    return True

def fill_chl_on_time(files):
    


def single_bloom_approach():
    """
    """
    return True

def two_bloom_approach():
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser
    parser.add_argument("blooms", options=[1,2])
    parser.add_argument("filelist")
    parser.add_argument("file_directory")
    args = parser.parse_args()