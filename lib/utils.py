# standard library
from subprocess import run, PIPE
from typing import List
import os
import shutil



RUN_CMD_ONFAIL_EXITCODE = 22

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

def run_bash(bash_code: str):
    """Safely runs a given bash code, and nicely fails on a non-zero exit code"""
    return run_cmd(['bash', '-c', bash_code])

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
