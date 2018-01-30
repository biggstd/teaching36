# -*- coding: utf-8 -*-
"""A Python implementation of Conway's Game of Life.

This module demonstrantes an implementation of Conway's Game of Life in
Python. It is meant to be run, and its output viewed, from the terminal.

.. moduleauthor:: Tyler Biggs <biggstd@gmail.com>

"""


class Cell:
    """Model for a single cell within the Game of Life.
    """

    def __init__(self):
        """Initialization function for a cell.

        :param x:
            The X-coordinate.

        :param y:
            The Y-coordinate.

        :param state:
            The state of the cell.
        """

    def get_state(self):
        """Return the current state of the cell."""
        pass

    def set_state(self):
        """Set the state of this cell."""
        pass


def test_cell_class():
    """Test function for the Cell class."""

    # Initialize a Cell object.
    test_cell = Cell()

    # Test the get_state() function.
    test_cell.get_state()

    # Test the set_state() function.
    test_cell.set_state()


class Grid:
    """Holds a collection of Cell objects in a grid.

    """

    def __init__(self):
        """Grid initialization function.


        """
        pass

    def print_grid(self):
        """Prints the grid to the terminal."""
        pass

    def next_move(self):
        """Blank for now.
        """
        pass

    def set_cell(self, x, y, state):
        """Set the cell at index (x, y) to a new state.
        """
        pass

    def play(self):
        """
        """
        pass
