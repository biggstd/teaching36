# This script is a solution for the in-class assignment for
# weeks 4 and 5.
#
class cell(object):
    # The cell class is used to maintain the state of
    # any given cell in an nxm grid for Conway's game of life.
    x = 0
    y = 0
    state = 0

    def __init__(self, x, y, state):
        # Initialize the class variables
        self.state = state
        self.x = x
        self.y = y
        self.state = state

    def get_state(self):
        # A "getter" function to return the current state of the cell.
        return self.state

    def set_state(self, state):
        # A "setter" function to set the current state of the cell.
        self.state = state


def test_cell():
    # A unit test to make sure the cell class is working properly.
    temp_cell = cell(4,5, 1)
    print(temp_cell.get_state())
    temp_cell.set_state(0)
    print(temp_cell.get_state())

class grid(object):
    # The grid class serves as the grid for the conway game of life. It
    # consists of an nxm matrix of cells.
    grid = []
    rows = 10
    cols = 80

    def __init__(self, rows, cols):
        # Intialize the grid to the size requested by the end-user.
        # Fill each cell of the grid with a a cell object and
        # set its default tate to off.
        self.rows = rows
        self.cols = cols
        self.grid = []
        for i in range(0, rows):
            col_list = []
            for j in range(0, cols):
                col_list.append(cell(i, j, 0))
            self.grid.append(col_list)

    def print_grid(self):
        # Print the grid to the screen.
        for i in range(0, self.rows):
            for j in range(0, self.cols):
                state = self.grid[i][j].get_state()
                if state == 1:
                    print("X", end='')
                else:
                    print(" ", end='')
            print("\n", end='')


    def set_cell(self, x, y, state):
        # Sets any cell within the grid, using an x,y coordinate, to the
        # state requested by the caller.
        self.grid[x][y].set_state(state)


    def next_move(self):
        # Create a new grid of the current states.
        curr_state = []
        for i in range(0, self.rows):
            col_list = []
            for j in range(0, self.cols):
                col_list.append(self.grid[i][j].get_state())
            curr_state.append(col_list)

        # Get the neighborhood around all cells in the grid except along the
        # borders: range of 1:rows-1 and 1:cols-1.
        for i in range(1, self.rows - 1):
            for j in range(1, self.cols - 1):
                # Get the state of this cell and it's 8 neighbors.
                nw = curr_state[i-1][j-1]
                n = curr_state[i-1][j]
                ne = curr_state[i-1][j+1]
                w = curr_state[i][j-1]
                me = curr_state[i][j]
                e = curr_state[i][j+1]
                sw = curr_state[i+1][j-1]
                s = curr_state[i+1][j]
                se = curr_state[i+1][j+1]

                # Calculate the number of neighbors that are live.
                num_neighbors = ne + n + nw + e + w + se + s + sw

                # Any live cell with fewer than two live neighbours dies, as if
                # caused by underpopulation.
                if num_neighbors < 2 and me == 1:
                    self.grid[i][j].set_state(0)
                # Any live cell with two or three live neighbours lives on to
                # the next generation.
                elif (num_neighbors == 2 or num_neighbors == 3) and me == 1:
                    self.grid[i][j].set_state(1)
                # Any live cell with more than three live neighbours dies, as
                # if by overpopulation.
                elif num_neighbors > 3 and me == 1:
                    self.grid[i][j].set_state(0)
                # Any dead cell with exactly three live neighbours becomes a
                # live cell, as if by reproduction.
                elif num_neighbors == 3 and me == 0:
                    self.grid[i][j].set_state(1)

    def play(self, ticks = 100):
        for i in range(0, ticks):
            print(f"Step {i+1}")
            self.next_move()
            # self.print_grid()

def test_grid():
    # A unit test to make sure our grid is working.
    temp_grid = grid(30, 80)
    temp_grid.set_cell(29, 79, 1)
    temp_grid.print_grid()
#     # TODO: add a test for the grid.


#test_cell()
#test_grid()

# Intialize our grid, add some cells that are turned on and then
# play
gofl = grid(30, 80)
gofl.set_cell(14, 40, 1)
gofl.set_cell(15, 42, 1)
gofl.set_cell(16, 39, 1)
gofl.set_cell(16, 40, 1)
gofl.set_cell(16, 43, 1)
gofl.set_cell(16, 44, 1)
gofl.set_cell(16, 45, 1)
gofl.play(100)

gofl.print_grid()
