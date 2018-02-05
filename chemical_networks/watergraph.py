"""Contains the WaterGraph class. A container class for a networkx graph of
a water coordinate file.

"""

# Import basic Python tools.
import itertools
import collections
import random
import functools

# Import data science packages.
import numpy as np
import networkx as nx

# Import plotting tools.
import matplotlib.pyplot as plt
from plotly.offline import init_notebook_mode, plot, iplot
import plotly.graph_objs as go


class WaterGraph:
    """Models a single `snapshot.xyz` file.

    """

    ELEMENT_COLORS = {
        'H': 'red',
        'O': 'blue',
    }

    ELEMENT_SIZES = {
        'H': 25,
        'O': 60,
    }

    def __init__(self, snapshot_filepath):
        """
        """

        # Create the graph from the given snapshot file.
        self.graph = self.create_graph_from_xyz(snapshot_filepath)

    def create_graph_from_xyz(self, xyz_filepath):
        """Creates a new networkx graph object from the given file.

        :param xyz_filepath:
            The filepath, as a string, from which the network will
            be constructed.

        :returns:
            A networkx graph object.
        """

        # Open the file and read the lines therein.
        with open(xyz_filepath) as in_file:

            # Drop the first two lines of the .xyz file.
            data = in_file.readlines()[2:]

        # Create the empty network to be returned.
        graph = nx.Graph()

        # Iterate through the data, and generate a node for each line.
        for index, line in enumerate(data):

            # Split the string by whitespace.
            line = line.split()

            # Get the element type.
            node_element = line[0]

            # Get the x, y, z coordinates.
            x_val = np.float(line[1])
            y_val = np.float(line[2])
            z_val = np.float(line[3])

            # Create the node within the new graph.
            graph.add_node(
                # Save nodes by index.
                node_for_adding=index,

                # Remaining arguments passed here will be passed to
                # the node as key: value pair attributes.
                element=node_element,
                x=x_val,
                y=y_val,
                z=z_val,
                size=self.ELEMENT_SIZES[node_element],
                color=self.ELEMENT_COLORS[node_element]
            )

        # Return the graph generated.
        return graph

    def calculate_distance_matrix(self):
        """Calculate and return a distance matrix based a given graph.

        :param in_graph:
            A networkx graph object. The nodes in this graph must have
            x, y, and z attributes for this function to work.

        :returns:
            A distance matrix as a numpy array.
        """

        # Find the length / number of atoms in the system.
        atoms_total = len(self.graph)

        # Create a matrix of zeros the size of our edges.
        distance_matrix = np.zeros((atoms_total, atoms_total))

        # Get the x, y, z coordinates from the graph nodes.
        xyz_vals = [self.get_node_xyz(n) for n in self.graph]

        # We only need concern ourselves with the upper triagnle indicies.
        it, jt = np.triu_indices(atoms_total)

        # Iterate through the upper triangle indexes and calculate
        # the distance between each set of points.
        for i, j in [x for x in zip(it, jt)]:

            coord_a = xyz_vals[i]
            coord_b = xyz_vals[j]

            distance_matrix[i, j] = np.linalg.norm(coord_a - coord_b)

        # Add the transpose of the matrix to the original matrix.
        out = distance_matrix.T + distance_matrix
        np.fill_diagonal(out, np.diag(distance_matrix))
        return out

    def get_node_xyz(self, in_node):
        """Returns a tuple of x, y, z values from a given
        graph and node pair. These attributes must already be
        defined in the node.

        :param in_node:
            A networkx node object contained within the graph
            input above.

        :param graph:
            A networkx graph object.

        :returns:
            A tuple of the x, y, z coordinates of the node.
        """

        # Find the coordinates.
        xx = self.graph.node('x')[in_node]
        yy = self.graph.node('y')[in_node]
        zz = self.graph.node('z')[in_node]

        # Return an array of the values.
        return np.array((xx, yy, zz))

    def add_distance_array_to_nodes(self):
        """Stores a 1D array of distances to other atoms on the
        corresponding node within graph. The distance matrix
        creation function will be called by this function.

        The graph and distance_matrix objects must be indexed
        the same.

        :param graph:
            A networkx graph object. Node labels should be indexes.

        :returns:
            Modified the graph in place. Appends distance arrays
            to each node.
        """

        # Calculate the distance matrix for the given graph.
        dist_mat = self.calculate_distance_matrix()

        # Iterate through the graph nodes and each array within
        # the newly calculated distance matrix.
        for nd, da in zip(list(self.graph), dist_mat):

            # Add the distance array attribute to the node.
            self.graph.add_node(nd, distance_array=da)


    def add_conditional_edge(self, predicates):
        """Iterates through each node within a graph, checks a
        boolean predicate (which can examine another, or all of
        the other nodes) and adds an edge if the predicates all
        return true for a given node pair.

        :param predicates:
            A list of boolean functions that examine a pair of nodes.

        """

        # Get the list of all nodes contained within the given graph.
        node_list = list(self.graph)

        # Generate all possible combinations of nodes of length 2.
        node_combinations = itertools.combinations(node_list, 2)

        # Iterate over all of the node combinations.
        for n1, n2 in node_combinations:

            # Check all of the predicates, pass the two nodes being
            # examined, also pass the graph so other attributes can
            # be viewed by the predicate functions.
            if all(p(n1, n2, self.graph) for p in predicates):

                # All predicates have returned True if the function
                # gets to this point. So add the new edge.
                self.graph.add_edge(n1, n2)

    def calculate_dipole_vectors(self):
        """Calculate the dipole vectors for waters found within graph.

        This functinon assumes that the paths / edges between the have
        already be created.

        :param graph:
            A graph object. Must have the paths from oxygen to hydrogen
            already in place.

        :returns:
            Modifies the given graph in place.
        """

        # Iterate through those atoms within the graph.
        for nd in self.graph.nodes():

            # We only care about nodes that are oxygen.
            if self.graph.node('element')[nd] == 'O':

                # Get the index of this oxygen.
                oxy_idx = self.get_node_xyz(nd)

                # Ensure that we only consider oxygens with only two
                # neighbors.
                if len(self.graph[nd]) == 2:

                    # Get the neighbors of the oxygen node. These are the
                    # graph indixes of the hydrogens.
                    h_nbrs = [self.get_node_xyz(n) for n in self.graph[nd]]

                    # Calculate the mean of these two vectors. Append this
                    # value to the current node.
                    dpv = np.dot(
                        [np.array([hc for hc in h_nbrs[0] - oxy_idx])],
                        [np.array([hc for hc in h_nbrs[1] - oxy_idx])]
                    )
                    self.graph.add_node(node_for_adding=nd, dipole_vector=dpv)
