

def element_check(node1, node2, graph, element_set):
    """Examines the elements of node1 and node2 to see
    if they are equal to the set provided in element_set.

    :param node1:
        Index value of the first node.

    :param node2:
        Index value of the second ndoe.

    :param graph:
        A networkx graph object that contains node1 and node2.

    :param element_set:
        The set of elements we wish to compare to.

    :returns:
        True or False, based on the element type of the nodes.

    """

    # Get the element type from each node.
    node_1_element = graph.node('element')[node1]
    node_2_element = graph.node('element')[node2]

    # Build a set from the element types.
    node_set = set([node_1_element, node_2_element])

    # Compare a set of those elements to the given element_set.
    if node_set == element_set:
        return True

def distance_check(node1, node2, graph, low, high):
    """Examines node1 and ndoe2 within a given graph.

    :param node1:
        Index value of the first node.

    :param node2:
        Index value of the second ndoe.

    :param graph:
        A networkx graph object that contains node1 and node2.

    :param low:
        The lowest acceptable distance, not inclusive.

    :param high:
        The highest acceptable distance, not inclusive.

    :returns:
        True or False, based on the distance between nodes.
    """

    # Get the distance array from the first node.
    node1_distance_array = graph.node('distance_array')[node1]

    # The distance between the two nodes can be looked up from
    # the distance_array by using the index of the second ndoe.
    n1_n2_dist = node1_distance_array[node2]

    # Now compare the distance to the given ranges.
    if low < n1_n2_dist < high:
        return True
