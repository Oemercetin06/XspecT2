"""Module for connecting with the PubMLST database via API requests and downloading allele files."""

__author__ = "Cetin, Oemer"

import requests
import json
from xspect.mlst_feature.mlst_helper import (
    create_fasta_files,
    pick_species_number_from_db,
    pick_scheme_number_from_db,
    pick_scheme,
    scheme_list_to_dict,
)
from xspect.definitions import get_xspect_mlst_path, get_xspect_upload_path


class PubMLSTHandler:
    """Class for communicating with PubMLST and downloading alleles (FASTA-Format) from all loci."""

    base_url = "https://rest.pubmlst.org/db"

    def __init__(self):
        """Initialise a PubMLSTHandler object."""
        # Default values: Oxford (1) and Pasteur (2) schemes of A.baumannii species
        self.scheme_list = [
            self.base_url + "/pubmlst_abaumannii_seqdef/schemes/1",
            self.base_url + "/pubmlst_abaumannii_seqdef/schemes/2",
        ]
        self.scheme_paths = []
        self.scheme_mapping = {}

    def get_scheme_paths(self) -> dict:
        """
        Get the scheme paths in a dictionary.

        Returns:
            dict: A dictionary containing the scheme paths.
        """
        return scheme_list_to_dict(self.scheme_paths)

    def choose_schemes(self) -> None:
        """
        Changes the scheme list attribute to feature other schemes from another species.

        This function lets the user pick schemes to download all alleles that belong to it.
        The scheme has to be available in the database.
        """
        available_species = {}
        available_schemes = {}
        chosen_schemes = []
        counter = 1
        # retrieve all available species
        species_url = PubMLSTHandler.base_url
        for species_databases in requests.get(species_url).json():
            for database in species_databases["databases"]:
                if database["name"].endswith("seqdef"):
                    available_species[counter] = database["name"]
                    counter += 1
        # pick a species out of the available ones
        chosen_species = pick_species_number_from_db(available_species)

        counter = 1
        scheme_url = f"{species_url}/{chosen_species}/schemes"
        for scheme in requests.get(scheme_url).json()["schemes"]:
            # scheme["description"] stores the name of a scheme.
            # scheme["scheme"] stores the URL that is needed for downloading all loci.
            available_schemes[counter] = [scheme["description"], scheme["scheme"]]
            counter += 1

        # Selection process of available scheme from a species for download (doubles are caught!)
        while True:
            chosen_scheme = pick_scheme_number_from_db(available_schemes)
            (
                chosen_schemes.append(chosen_scheme)
                if chosen_scheme not in chosen_schemes
                else None
            )
            choice = input(
                "Do you want to pick another scheme to download? (y/n):"
            ).lower()
            if choice != "y":
                break
        self.scheme_list = chosen_schemes

    def download_alleles(self, choice: False) -> None:
        """
        Downloads every allele FASTA-file from all loci of the scheme list attribute.

        This function sends API-GET requests to PubMLST.
        It downloads all alleles based on the scheme_list attribute.
        The default schemes are the Oxford and Pasteur schemes of A.baumannii

        Args:
            choice (bool): The decision to download different schemes, defaults to False.
        """
        if choice:  # pick an own scheme if not Oxford or Pasteur
            self.choose_schemes()  # changes the scheme_list attribute

        for scheme in self.scheme_list:
            scheme_json = requests.get(scheme).json()
            # We only want the name and the respective featured loci of a scheme
            scheme_name = scheme_json["description"]
            locus_list = scheme_json["loci"]

            species_name = scheme.split("_")[1]  # name = pubmlst_abaumannii_seqdef
            scheme_path = get_xspect_mlst_path() / species_name / scheme_name
            self.scheme_mapping[str(scheme_path)] = scheme
            self.scheme_paths.append(scheme_path)

            for locus_url in locus_list:
                # After using split the last part ([-1]) of the url is the locus name
                locus_name = locus_url.split("/")[-1]
                locus_path = (
                    get_xspect_mlst_path() / species_name / scheme_name / locus_name
                )

                if not locus_path.exists():
                    locus_path.mkdir(exist_ok=True, parents=True)

                alleles = requests.get(f"{locus_url}/alleles_fasta").text
                create_fasta_files(locus_path, alleles)

    def assign_strain_type_by_db(self) -> None:
        """
        Sends an API-POST-Request to the database for MLST without bloom filters.

        This function sends API-POST requests to PubMLST.
        It is a different way to determine strain types based on a BLAST-Search.
        This function is only used for testing and comparing results.
        """
        scheme_url = (
            str(pick_scheme(scheme_list_to_dict(self.scheme_list))) + "/sequence"
        )
        fasta_file = get_xspect_upload_path() / "Test.fna"
        with open(fasta_file, "r") as file:
            data = file.read()
            payload = {  # Essential API-POST-Body
                "sequence": data,
                "filetype": "fasta",
            }
        response = requests.post(scheme_url, data=json.dumps(payload)).json()

        for locus, meta_data in response["exact_matches"].items():
            # meta_data is a list containing a dictionary, therefore [0] and then key value.
            # Example: 'Pas_fusA': [{'href': some URL, 'allele_id': '2'}]
            print(locus + ":" + meta_data[0]["allele_id"], end="; ")
        print("\nStrain Type:", response["fields"])

    def get_strain_type_name(self, highest_results: dict, post_url: str) -> str:
        """
        Send an API-POST request to PubMLST with the highest result of each locus as payload.

        This function formats the highest_result dict into an accepted input for the request.
        It gets a response from the site which is the strain type name.
        The name is based on the allele id with the highest score for each locus.
        Example of post_url for the oxford scheme of A.baumannii:
        https://rest.pubmlst.org/db/pubmlst_abaumannii_seqdef/schemes/1/designations

        Args:
            highest_results (dict): The allele ids with the highest kmer matches.
            post_url (str): The specific url for the scheme of a species

        Returns:
            str: The response (ST name or No ST found) of the POST request.
        """
        payload = {
            "designations": {
                locus: [{"allele": str(allele)}]
                for locus, allele in highest_results.items()
            }
        }

        response = requests.post(post_url + "/designations", json=payload)

        if response.status_code == 200:
            data = response.json()
            if "fields" in data:
                post_response = data["fields"]
                return post_response
            else:
                post_response = "No matching Strain Type found in the database. "
                post_response += "Possibly a novel Strain Type."
                return post_response
        else:
            post_response = "Error:" + str(response.status_code)
            post_response += response.text
            return post_response
