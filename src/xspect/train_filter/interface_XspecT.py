import os
from pathlib import Path
from shutil import rmtree
from xspect.definitions import get_xspect_model_path, get_xspect_tmp_path


def make_paths(dir_name, genus):
    """Create paths to the concatenated sequences and to where the new bloomfilters will be saved.

    :param dir_name: Name of the parent directory.
    :type dir_name: str
    :param genus: Name of the genus.
    :type genus: str
    :return: The path to the sequence files and the bloomfilter directory.
    """
    # Path to concatenated sequences
    files_path = get_xspect_tmp_path() / dir_name / "concatenate"

    # Path for results.
    result_path = get_xspect_model_path() / genus

    return str(files_path), str(result_path)


def save_time_stats(time_stats, dir_name):
    """Saves the collected time measurements as a txt file.

    :param time_stats: The collected time measurements as a formatted string.
    :type time_stats: str
    :param dir_name: Name of the parent directory.
    :type dir_name: str
    """
    time_file = get_xspect_tmp_path() / dir_name / "time.txt"
    with open(time_file, "w+", encoding="utf-8") as f:
        f.write(time_stats)
