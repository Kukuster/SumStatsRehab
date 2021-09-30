# standard library
from typing import Union, List
import os
import shutil
import gzip
from zipfile import ZipFile

# third-party libraries
import requests
import magic


def file_exists(path: str):
    return os.path.isfile(path)

def get_file_size_bytes(path: str) -> int:
    return os.path.getsize(path)


def file_size_human(size_in_B: int, verbose: bool = True):
    """
    Only displays either in MB or B. Gotta update this (probably)
    """
    size_in_MiB: float = (size_in_B/1024)/1024
    if size_in_MiB < 1:
        size_in_MiB_str = "<1 MB"
    else:
        size_in_MiB_str = f"{int(size_in_MiB+1)} MB"
    if verbose:
        return f"{size_in_MiB_str} ({size_in_B} B)"
    else :
        return size_in_MiB_str


def download_file(url: str, path: str):
    with requests.get(url, stream=True) as r:
        with open(path, 'wb') as f:
            shutil.copyfileobj(r.raw, f)
    return path

def ungzip(gz: str, out: str):
    with gzip.open(gz, 'rb') as f_in:
        with open(out, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def largest_file_from(files: List[str]) -> Union[str, None]:
    real_files = list(filter(file_exists, files))
    if len(real_files) == 0:
        return None
    return max(real_files, key=lambda f: get_file_size_bytes(f))




def open_maybe_gz(path: str):
    if not file_exists(path):
        raise FileNotFoundError(f"Failed attempt to open file at path: {path}")

    mime: str = magic.from_file(path, mime=True)

    if mime == 'text/plain':
        return open(path, 'r')
    elif mime == 'application/gzip' or mime == 'application/x-gzip':
        return gzip.open(path, 'rb')
    elif mime == 'inode/x-empty':
        raise FileNotFoundError('The input file is empty!')
    else:
        raise ValueError(
            f"Got an unexpected type of the input file: {mime}")



def resolve_bare_text_file(maybe_archive_path: str, unpacked_text_file_path: str):
    """
    Returns the path for the bare (unpacked) dataset file.

    Takes a path to the input file which can be an archive,
    and a desired unpacked file path if unpacking takes place.
    """

    if not file_exists(maybe_archive_path):
        raise FileNotFoundError(f"Failed attempt to open file at path: {maybe_archive_path}")

    mime: str = magic.from_file(maybe_archive_path, mime=True)

    if mime == 'application/gzip' or mime == 'application/x-gzip':
        print("the SS file is a gzip. Unpacking")
        with gzip.open(maybe_archive_path, 'rb') as f_in:
            try:
                with open(unpacked_text_file_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            except EOFError as e:
                print("This error occured while trying to unpack files from gzip:")
                print(e)
                print("it is possible that the dataset file is broken!")
                print("it is also possible that it's not :/")
                print("Please check manually, maybe everything's fine.")
                raise EOFError
        return unpacked_text_file_path
    elif mime == 'application/zip':
        print("the SS file is a zip. Unpacking")
        with ZipFile(maybe_archive_path, 'r') as zipObj:
            # Extract all the contents of zip file in current directory
            extract_dir = os.path.dirname(os.path.abspath(maybe_archive_path))
            print(f"extract_dir = {extract_dir}")

            zipObj.extractall(path=extract_dir)
            extracted = [os.path.join(extract_dir, f) for f in zipObj.namelist()]
            print("extracted files:")
            print(extracted)

            extracted_file_path = largest_file_from(extracted)

            if extracted_file_path is not None:
                unpacked_text_file_path = extracted_file_path
            else:
                raise ValueError(
                    "Error extracting dataset file from a downloaded zip. Check manually.")

            print("treating the following file as the dataset")
            print(f"\"{extracted_file_path}\"")
        return unpacked_text_file_path
    elif mime == 'text/plain':
        # the file is the uncompressed dataset file
        return maybe_archive_path
    elif mime == 'inode/x-empty':
        raise FileNotFoundError('The file is empty!')
    else:
        raise ValueError(
            f"Got unexpected type of the input file: {mime}")

