import plotly.graph_objs as go
import networkx as nx


def draw_plotly_3d_scatter(in_graph):
    """COMMENT ME!

    """

    # Get the coordinate attributes back out of the graph.
    x_vals = [x for x in nx.get_node_attributes(in_graph, 'x').values()]
    y_vals = [x for x in nx.get_node_attributes(in_graph, 'y').values()]
    z_vals = [x for x in nx.get_node_attributes(in_graph, 'z').values()]

    # Create the color and size lists.
    c_vals = [x for x in nx.get_node_attributes(in_graph, 'color').values()]
    s_vals = [x for x in nx.get_node_attributes(in_graph, 'size').values()]

    # Add the points to a scatter plot.
    trace1 = go.Scatter3d(
        x=x_vals,
        y=y_vals,
        z=z_vals,
        mode='markers',
        hovertext = [x for x in in_graph.nodes()],
        marker=dict(
            sizemode='area',
            color=c_vals,
            size=s_vals,
            opacity=0.8
        ),
    )

    # Create empty list for line segments
    Xe=[]
    Ye=[]
    Ze=[]

    # Iterate through the edges, and draw line segments for each set.
    for edge in in_graph.edges():

        # Assign the node edge indecies.
        node_a = edge[0]
        node_b = edge[1]

        # Assign the x, y, z variables.
        xa = in_graph.node('x')[node_a]
        ya = in_graph.node('y')[node_a]
        za = in_graph.node('z')[node_a]

        xb = in_graph.node('x')[node_b]
        yb = in_graph.node('y')[node_b]
        zb = in_graph.node('z')[node_b]

        # Add these values to a corresponding list.
        # The third value is none, otherwise a plane would be drawn.
        Xe += [xa, xb, None]
        Ye += [ya, yb, None]
        Ze += [za, zb, None]

    # Add the node edges / line segments to the plot.
    trace2 = go.Scatter3d(
        x=Xe,
        y=Ye,
        z=Ze,
        mode='lines',
        line=go.Line(color='rgb(125,125,125)', width=1),
        hoverinfo='none'
    )

    # Show the plot.
    data = [trace1, trace2]
    fig = go.Figure(data=data)
    # iplot(fig)
    return fig
