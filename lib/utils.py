# standard library
from subprocess import run, PIPE
from typing import Dict, List, TypedDict
import os
import shutil



RUN_CMD_ONFAIL_EXITCODE = 22

class CMD_RETURN:
    ec: int
    stdout: str
    stderr: str

def run_cmd(cmd: List[str]):
    """A wrapper around subprocess.run that nicely fails on a non-zero exit code"""
    if len(cmd) == 0:
        raise ValueError('cmd has to be a non-empty list')

    res = run(cmd, stdout=PIPE, stderr=PIPE)

    if res.returncode != 0:
        error_message = f"command \"{cmd[0]}\" finished with exit code: {res.returncode}"
        stderr = res.stderr.decode('utf-8')
        if stderr:
            error_message += "\nand produced the following error message:\n"
            error_message += stderr
        raise ChildProcessError(error_message)

    return res.stdout.decode('utf-8').rstrip()

def run_cmd_rich(cmd: List[str]) -> CMD_RETURN:
    """A wrapper around subprocess.run, returns stdout, stderr, and exit code"""
    if len(cmd) == 0:
        raise ValueError('cmd has to be a non-empty list')

    res = run(cmd, stdout=PIPE, stderr=PIPE)

    r = CMD_RETURN()
    r.ec     = res.returncode
    r.stdout = res.stdout.decode('utf-8').rstrip()
    r.stderr = res.stderr.decode('utf-8').rstrip()

    return r

def run_bash(bash_code: str):
    """Safely runs a given bash code, and nicely fails on a non-zero exit code"""
    return run_cmd(['bash', '-c', bash_code])

def run_bash_rich(bash_code: str):
    """Runs a given bash code, returns stdout, stderr, and exit code"""
    return run_cmd_rich(['bash', '-c', bash_code])


def rm(file: str):
    try:
        return os.remove(file)
    except FileNotFoundError:
        return None

def rm_r(dir: str):
    try:
        return shutil.rmtree(dir)
    except FileNotFoundError:
        return None

def rm_rf(file_or_dir: str):
    if os.path.isdir(file_or_dir):  
        return rm_r(file_or_dir)
    elif os.path.isfile(file_or_dir):  
        return rm(file_or_dir)
    else:
        return None

def mv(file_or_dir: str, new_name: str):
    shutil.move(file_or_dir, new_name)
