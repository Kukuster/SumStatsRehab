# standard library
from subprocess import run, PIPE
from typing import List



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
