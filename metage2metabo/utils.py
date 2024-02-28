# Copyright (C) 2019-2024 Clémence Frioux & Arnaud Belcour - Inria Dyliss - Pleiade - Microcosme
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

import sys
import tarfile
import os
import logging

logger = logging.getLogger(__name__)

def get_basename(filepath):
    """Return the basename of given filepath.
    
    Args:
        filepath (str): path to a file
    
    Returns:
        str: basename

    >>> basename('~/an/interesting/file.txt')
    'file
    """
    return os.path.splitext(os.path.basename(filepath))[0]


def get_extension(filepath):
    """Get the extension of a filepath
    
    Args:
        filepath (str): path to a file
    
    Returns:
        str: extention of the file

    >>> extension('~/an/interesting/file.lp')
    'lp'
    >>> extension('nothing')
    ''
    >>> extension('nothing.important')
    'important'
    """
    return os.path.splitext(os.path.basename(filepath))[1][1:]


def is_valid_path(filepath):
    """Return True if filepath is valid
    
    Args:
        filepath (str): path to file
    
    Returns:
        bool: True if path exists, False otherwise
    """
    if filepath and not os.access(filepath, os.W_OK):
        try:
            open(filepath, 'w').close()
            os.unlink(filepath)
            return True
        except OSError:
            return False
    else:  # path is accessible
        return True


def is_valid_file(filepath):
    """Return True if filepath exists

    Args:
        filepath (str): path to file

    Returns:
        bool: True if path exists, False otherwise
    """
    try:
        open(filepath, 'r').close()
        return True
    except OSError:
        return False


def is_valid_dir(dirpath):
    """Return True if directory exists or can be created (then create it)
    
    Args:
        dirpath (str): path of directory

    Returns:
        bool: True if dir exists, False otherwise
    """
    if not os.path.isdir(dirpath):
        try:
            os.makedirs(dirpath)
            return True
        except OSError:
            return False
    else:
        return True

def check_program(program):
    """Check whether Pathway Tools is in the PATH

    Returns:
        bool: True if Pathway Tools is in the PATH, False otherwise
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    if is_exe(program):
        return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True

    return False


def file_or_folder(variable_folder_file):
    """Check if the variable is file or a folder

    Args:
        variable_folder_file (str): path to a file or a folder

    Returns:
        dict: {name of input file: path to input file}
    """
    file_folder_paths = {}

    check_file = False
    if os.path.isfile(variable_folder_file):
        filename = os.path.splitext(os.path.basename(variable_folder_file))[0]
        file_folder_paths[filename] = variable_folder_file
        check_file = True

    check_folder = False
    # For folder, iterate through all files inside the folder.
    if os.path.isdir(variable_folder_file):
        for file_from_folder in os.listdir(variable_folder_file):
            filename = os.path.splitext(os.path.basename(file_from_folder))[0]
            file_folder_paths[filename] = os.path.join(variable_folder_file, file_from_folder)
            check_folder = True

    if check_file is False and check_folder is False:
        logger.critical('ERROR: Wrong input, {0} does not exit'.format(variable_folder_file))
        sys.exit(1)

    return file_folder_paths

def check_absolute_path(directory, target):
    """ Check if the extracted element is inside the output directory.
    If not, it is a potential path traversal attempt.

    Args:
        directory (str): path to output directory for extraction.
        target (str): path of file contained in tar file.
    """
    abs_directory = os.path.abspath(directory)
    abs_target = os.path.abspath(target)

    prefix = os.path.commonprefix([abs_directory, abs_target])

    return prefix == abs_directory


def safe_tar_extract_all(tar_file, outdir):
    """ Perform a sanitized check to ensure no file outside the output folder will be modified.

    Args:
        tar_file (str): path to tar file.
        outdir (str): path to output directory for extraction.
    """
    tar = tarfile.open(tar_file, "r:gz")
    if sys.version_info >= (3, 12):
        tar.extractall(outdir, filter='data')
    else:
        tar.extractall(outdir)

    for member in tar.getmembers():
        member_path = os.path.join(outdir, member.name)
        if not check_absolute_path(outdir, member_path):
            raise Exception("Attempted Path Traversal in Tar File")

    tar.extractall(outdir)

    tar.close()