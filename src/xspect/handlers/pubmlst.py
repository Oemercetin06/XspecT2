"""Module for connecting with the PubMLST database via API requests and downloading allele files."""

__author__ = "Cetin, Oemer"

from pathlib import Path
import requests
from xspect.file_io import create_fasta_files


class PubMLSTHandler:
    """Class for communicating with PubMLST and downloading alleles (FASTA-Format) from all loci."""

    def __init__(self, base_url: str = "https://rest.pubmlst.org/db"):
        """Initialize a PubMLSTHandler object."""

        self.base_url = base_url

    def get_available_organisms(self) -> list:
        """
        Get a list of available species from PubMLST.

        Returns:
            list: A list of available species names.
        """
        available_species = []
        species_url = self.base_url
        for species_databases in requests.get(species_url, timeout=10).json():
            for database in species_databases["databases"]:
                if database["name"].endswith("seqdef"):
                    species_name = database["name"].split("_")[1]
                    available_species.append(species_name)
        return available_species

    def get_available_schemes(self, species: str) -> list:
        """
        Get a list of available schemes for a given species from PubMLST.

        Args:
            species (str): The species name.
        Returns:
            list: A list of available scheme names.
        """
        available_schemes = []
        scheme_url = f"{self.base_url}/pubmlst_{species}_seqdef/schemes"
        for scheme in requests.get(scheme_url, timeout=10).json()["schemes"]:
            available_schemes.append(scheme["description"])
        return available_schemes

    def get_scheme_url(self, species: str, scheme: str) -> str:
        """
        Get the scheme URL for a given species and scheme name.

        Args:
            species (str): The species name.
            scheme (str): The scheme name.
        Returns:
            str: The scheme URL.
        Raises:
            ValueError: If the scheme is not found for the given species.
        """
        scheme_url = f"{self.base_url}/pubmlst_{species}_seqdef/schemes"
        for scheme_entry in requests.get(scheme_url, timeout=10).json()["schemes"]:
            if scheme_entry["description"] == scheme:
                return f"{scheme_entry['scheme']}"
        raise ValueError(f"Scheme '{scheme}' not found for species '{species}'.")

    def download_alleles(self, species: str, scheme: str, scheme_path: Path) -> None:
        """
        Downloads every allele FASTA-file from all loci of the chosen scheme for a species.

        This function sends API-GET requests to PubMLST.
        It downloads all alleles based on the chosen scheme and species.

        Args:
            species (str): The species name.
            scheme (str): The scheme name.
            scheme_path (Path): The path where the scheme alleles will be stored.
        """
        scheme_url = self.get_scheme_url(species, scheme)

        scheme_json = requests.get(scheme_url, timeout=10).json()
        locus_list = scheme_json["loci"]

        for locus_url in locus_list:
            locus_name = locus_url.split("/")[-1]
            locus_path = scheme_path / locus_name

            if not locus_path.exists():
                locus_path.mkdir(exist_ok=True, parents=True)

            alleles = requests.get(f"{locus_url}/alleles_fasta", timeout=10).text
            create_fasta_files(locus_path, alleles)

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

        response = requests.post(post_url + "/designations", json=payload, timeout=10)

        if response.status_code == 200:
            data = response.json()
            if "fields" in data:
                post_response = data["fields"]
                return post_response
            post_response = "No matching Strain Type found in the database. "
            post_response += "Possibly a novel Strain Type."
            return post_response
        post_response = "Error:" + str(response.status_code)
        post_response += response.text
        return post_response
