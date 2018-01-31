# -*- coding: utf-8 -*-
"""A Python implementation of Conway's Game of Life.

This module demonstrantes an implementation of Conway's Game of Life in
Python. It is meant to be run, and its output viewed, from the terminal.

.. moduleauthor:: Tyler Biggs <biggstd@gmail.com>

I am going to implement some of these classes with properties. Se the
links below for some reading on what properties are in Python.

https://docs.python.org/3.6/howto/descriptor.html#properties
https://stackoverflow.com/a/17330273
https://www.python-course.eu/python3_properties.php

"""

# Import the standard number / math package for data science
# and math in general. Numpy.
import numpy as np


class Cell:
    """Model for a single cell within the Game of Life.
    """

    def __init__(self, x_coord=0, y_coord=0, state=0):
        """Initialization function for a cell.

        :param x:
            The X-coordinate. Defaults to 0.

        :param y:
            The Y-coordinate. Defaults to 0.

        :param state:
            The state of the cell. Valid values are 1, or 0.
            Defaults to 0.
        """

        # Declare class attributes.
        self.x_coord = x_coord
        self.y_coord = y_coord

        # Declare private attributes for use with the state property.
        self._state = state

    @property
    def state(self):
        """Getter for the state property.

        The default deleter for this property will not be overridden.

        :returns:
            Returns the cells current state.
        """
        return self._state

    @state.setter
    def state(self, state):
        """Setter for the state property. The name of this function must be
        the same as the function after our `@propert` statement.

        :param state:
            The state attribute of the cell. This value can only be 1 or 0.
        """

        # Set the private class attribute.
        self._state = state

    def __str__(self):
        """Create / Override the default output to be generated for display
        using Python's `print()` function.
        """
        return f'x: {self.x_coord}, y: {self.y_coord}, state: {self.state}.'


def test_cell_class():
    """Test function for the Cell class."""

    # Initialize a Cell object.
    test_cell_111 = Cell(1, 1, 1)

    # Print it out.
    print(test_cell_111)

    # Test the state setter function.
    test_cell_111.state = 1
    test_cell_111.state = 0



def get_neighbors(x, y):
    """Returns the eight neighbors of a point upon a grid."""
    return [
        [x - 1, y - 1],  [x, y - 1],  [x + 1, y - 1],
        [x - 1, y    ],               [x + 1, y    ],
        [x - 1, y + 1],  [x, y + 1],  [x + 1, y + 1],
    ]


class Grid:
    """Holds a collection of Cell objects in a grid.

    Using a numpy array. Read more about it here:

    https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
    """



    def __init__(self, rows=15, cols=20):
        """Grid initialization function.

        :param rows:
            Number of rows in the grid.

        :param cols:
            Number of columns in the grid.
        """

        # Declare class attributes.
        self.rows = rows
        self.cols = cols
        self.grid = list()

        # Create the cell objects.
        for r in range(self.rows):
            # Create a temporary list for the current row.
            col_list = list()

            for c in range(self.cols):
                # Append a new Cell instance with the r, c index values.
                col_list.append(Cell(r, c, 0))

            # Append the column list to the grid list.
            self.grid.append(col_list)

        # Create a private attribute for a grid_state property.
        # This value will be filled out by the getter function
        # when it is called.
        self._grid_state = None

    def build_grid_state(self):
        """Getter for the grid_state property."""

        # Define a list for the current state.
        current_state = list()

        # Iterate through the cells as before.
        for r in range(self.rows):
            current_row = list()

            for c in range(self.cols):

                # This time append the state to the list.
                current_row.append(self.grid[r][c].state)

            current_state.append(current_row)

        return current_state

    def print_grid(self):
        """Prints the grid to the terminal."""
        # Declare some variables for printing.
        PAD = ' '
        CELL_DISPLAY = {
            1: '%',
            0: '.',
        }
        # Create a list to hold rows to print.
        row_list = list()

        # Create each row to be printed.
        for x in range(self.rows):

            row_string = str()

            for y in range(self.cols):
                # Get the Cell at the current index, and call
                # its state property getter like it is an attribute.
                row_string += CELL_DISPLAY[self.grid[x][y].state] + PAD

            # Add the row to the row list.
            row_list.append(row_string)

        # Join the items in row_list by a new line character, and return it.
        print('\n'.join(row for row in row_list))


    def next_move(self):
        """Blank for now.
        """
        # Build the indicies to be simulated. Trim 1 off.
        xx = range(1, self.rows - 1)
        yy = range(1, self.cols - 1)

        # Build the current state.
        current_state = self.build_grid_state()

        # Iterate through the buffered coordinates.
        for x in xx:
            for y in yy:
                # Get the current cells state.
                state = current_state[x][y]

                # Get the neighbor indicies.
                neighbor_i = get_neighbors(x, y)

                # Get the number of neighbors that are alive.
                alive_n = sum([current_state[i][j] for i, j in neighbor_i])

                # Look at the rules to deterine this Cells fate. You could
                # let it live forever. But you wont. You monster.
                if alive_n < 2 and state == 1:
                    self.grid[x][y].state = 0

                # Staying alive with a little help from my friends.
                elif (alive_n == 2 or 3) and state == 1:
                    self.grid[x][y].state = 1

                elif alive_n > 3 and state == 1:
                    self.grid[x][y].state == 1

                elif alive_n == 3 and state == 0:
                    self.grid[x][y].state = 1

    def set_cell_state(self, x, y, state):
        """Set the cell at index (x, y) to a new state."""
        # Easy, we access our setter function like it is an attribute.
        self.grid[x][y].state = state

    def play(self, ticks):
        """Force time upon this peacefull domain.

        :param ticks:
            The number of eons to simulate.
        """

        for age in range(ticks):
            self.next_move()
            self.print_grid()
