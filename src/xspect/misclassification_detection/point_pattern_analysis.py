"""
Point pattern density analysis tool for the alignment-based misclassification detection.

Notes:
Developed by Oemer Cetin as part of Bsc thesis (2025), Goethe University Frankfurt am Main.
(An Integration of Alignment-Free and Alignment-Based Approaches for Bacterial Taxon Assignment)
"""

import numpy

__author__ = "Cetin, Oemer"


class PointPatternAnalysis:
    """Class for all point pattern density analysis procedures."""

    def __init__(self, points: list[int], length: int):
        """
        Initialise the class for point pattern analysis.

        This method sets up the required list with data points (sorted) and the length of the reference genome.
        All required intensity for the statistics is also calculated.

        Args:
            points (list): The start coordinates of mapped regions on the genome.
            length (int): The length of the reference genome.
        """
        if len(points) < 2:
            raise ValueError("Need at least 2 points.")
        self.sorted_points = numpy.sort(numpy.asarray(points, dtype=float))
        self.n = len(points)
        self.length = float(length)

    def ripleys_k(self) -> tuple[bool, float, float]:
        """
        Calculates the K-function for the given point distribution.

        This method calculates the K-function to describe the point distribution.
        The result is than compared with what would be expected under a completely random distribution.
        (Under complete randomness the K-function result is 2*r)

        Returns:
            tuple: A tuple containing the information whether points are clustered or not.
        """
        r = 0.01 * self.length
        left = 0
        right = 0
        total_neighbors = 0

        for i in range(self.n):
            while self.sorted_points[i] - self.sorted_points[left] > r:
                left += 1
            if right < i:
                right = i
            while (
                right + 1 < self.n
                and self.sorted_points[right + 1] - self.sorted_points[i] <= r
            ):
                right += 1
            total_neighbors += right - left
        k = (self.length / (self.n * (self.n - 1))) * total_neighbors
        return (k > 2 * r), k, 2 * r

    def ripleys_k_edge_corrected(self) -> tuple[bool, float, float]:
        """
        Calculates the K-function for the given point distribution with an edge correction factor.

        This method calculates the K-function to describe the point distribution.
        This time an additional factor is multiplied for each data point to account for edge effects.
        The result is than compared with what would be expected under a completely random distribution.
        (Under complete randomness the K-function result is 2*r)

        Returns:
            tuple: A tuple containing the information whether the points are clustered or not.
        """
        r = 0.01 * self.length
        left = 0
        right = 0
        total_weighted = 0

        for i in range(self.n):
            while self.sorted_points[i] - self.sorted_points[left] > r:
                left += 1
            if right < i:
                right = i
            while (
                right + 1 < self.n
                and self.sorted_points[right + 1] - self.sorted_points[i] <= r
            ):
                right += 1

            neighbors = right - left
            if neighbors > 0:
                a = max(0, self.sorted_points[i] - r)
                b = min(self.length, self.sorted_points[i] + r)
                overlap = b - a
                weight = (2 * r) / overlap if overlap > 0 else 0

                total_weighted += weight * neighbors

        k = (self.length / (self.n * (self.n - 1))) * total_weighted
        return (bool(k > 2 * r)), float(k), 2 * r
