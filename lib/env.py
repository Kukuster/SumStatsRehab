# standard library
import os
from typing import Literal



class WrongRuntimeEnvironmentVariable(Exception):
    pass


GWASSS_BUILD_NUMBER_ENV = 'build_num'


def get_build() -> Literal['hg38', 'hg19', 'hg18']:
    build = os.getenv(GWASSS_BUILD_NUMBER_ENV)
    if build == None:
        return 'hg38'
    elif build == 'hg38' or build.lower() == 'grch38':
        return 'hg38'
    elif build == 'hg19' or build.lower() == 'grch37':
        return 'hg19'
    elif build == 'hg18' or build.lower() == 'ncbi36':
        return 'hg18'
    else:
        raise WrongRuntimeEnvironmentVariable(f'got unknown GWAS SS build: \"{build}\"')


def set_build(build) -> Literal['hg38', 'hg19', 'hg18']:
    if build == None:
        return 'hg38'
    elif build == 'hg38' or build.lower() == 'grch38' or str(build) == '38':
        os.environ[GWASSS_BUILD_NUMBER_ENV] = 'hg38'
    elif build == 'hg19' or build.lower() == 'grch37' or str(build) == '37':
        os.environ[GWASSS_BUILD_NUMBER_ENV] = 'hg19'
    elif build == 'hg18' or build.lower() == 'ncbi36' or str(build) == '36':
        os.environ[GWASSS_BUILD_NUMBER_ENV] = 'hg18'
    else:
        raise WrongRuntimeEnvironmentVariable(f'got unknown GWAS SS build: \"{build}\"')

    return get_build()

