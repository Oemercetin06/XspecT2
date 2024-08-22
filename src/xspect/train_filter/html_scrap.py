""" HTML Scraping for the taxonomy check results from NCBI."""

__author__ = "Berger, Phillip"

import datetime
import pickle
import sys
import time
import requests
from loguru import logger
from xspect.definitions import get_xspect_root_path


class TaxonomyCheck:
    """Class to get the GCFs that passed the taxonomy check from NCBI."""

    _main_url = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/ANI_report_prokaryotes.txt"
    _test_url = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/"
    _main_path = get_xspect_root_path() / "taxonomy_check.txt"
    _test_path = get_xspect_root_path() / "tax_check_date.txt"
    _new_time: list
    _tax_check_ok = []

    def __init__(self):
        old_time = self._get_old_tax_date()
        self._new_time = self._get_new_tax_date()
        # Both Dates could be found.
        # Check if the html file was updated since the last time it was downloaded.
        # If yes than update the file.
        if self._new_time and old_time:
            if self._new_time == old_time:
                logger.info("File was not updated")
                self._get_old_file()
            else:
                logger.info("Updating file")
                self._update_tax_check()

        # The old date does not exist.
        # Get the html file for the taxonomy check results.
        elif self._new_time and not old_time:
            logger.info("No file was found. Creating new file")
            self._update_tax_check()

        elif not self._new_time and old_time:
            self._get_old_file()

        else:
            logger.error("Nothing was found")
            logger.error("Aborting")
            sys.exit()

    def _get_old_tax_date(self):
        try:
            with open(self._test_path, "rb") as f:
                old_time = pickle.load(f)
                return old_time
        except FileNotFoundError:
            return None

    def _get_new_tax_date(self):
        raw_response = requests.get(self._test_url, timeout=5)
        data = raw_response.text.split("\n")
        for line in data:
            if "ANI_report_prokaryotes.txt" in line:
                line_parts = line.split()
                date_parts = line_parts[-3].split("-")
                date = datetime.date(
                    int(date_parts[0]), int(date_parts[1]), int(date_parts[2])
                )
                time_parts = line_parts[-2].split(":")
                combined_time = datetime.time(int(time_parts[0]), int(time_parts[1]))
                new_time = [date, combined_time]

                return new_time

        return None

    def _update_tax_check(self):
        raw_response = requests.get(self._main_url, timeout=5)
        all_tax_checks = raw_response.text.split("\n")[1:-1]
        self._get_gcf_ok(all_tax_checks)
        self._save_time()
        self._save_file()

    def _get_gcf_ok(self, all_tax_checks: list):
        tax_check_ok = []
        for line in all_tax_checks:
            line_parts = line.split("\t")
            gcf = line_parts[1]
            tax_check_status = line_parts[-1]
            if tax_check_status == "OK":
                tax_check_ok.append(gcf)

        self._tax_check_ok = tax_check_ok

    def _save_time(self):
        with open(self._test_path, "wb") as f:
            pickle.dump(self._new_time, f)

    def _save_file(self):
        with open(self._main_path, "wb") as f:
            pickle.dump(self._tax_check_ok, f)

    def _get_old_file(self):
        with open(self._main_path, "rb") as f:
            self._tax_check_ok = pickle.load(f)

    @staticmethod
    def _get_current_time():
        return time.asctime(time.localtime()).split()[3]

    def ani_gcf(self):
        """Returns ANI GCFs that passed the taxonomy check."""
        return self._tax_check_ok
