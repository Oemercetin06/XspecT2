"""This module contains functions that return file paths to the currently needed metagenome assembly."""

from xspect.definitions import get_xspect_tmp_path


def get_current_dir_file_path(dir_name):
    """Returns str of file path to the directory with the currently needed metagenome assembly.

    :param dir_name: Name of the current genus_metadata directory.
    :type dir_name: str
    :return: File path to the metagenome assembly.
    """
    return get_xspect_tmp_path() / dir_name
