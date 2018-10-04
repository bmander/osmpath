"""osmpath - Shortest paths using OpenStreetMap."""

from osmread import parse_file, Way, Node
from collections import Counter
import pyproj

__version__ = '0.1.0'
__author__ = 'Brandon Martin-Anderson <badhill@gmail.com>'
__all__ = []

geod = pyproj.Geod(ellps='clrk66')

def chop(ary, endpoints):
    """Returns a generator of segments of `ary` bounded by values `endpoints`.

    For example: `chop([1,2,3,4,5,6,7], [2,5,7])` -> `[[2,3,4,5],[5,6,7]]`
    """

    i_beg = None
    
    for i, el in enumerate(ary):
        if el in endpoints:
            if i_beg is not None:
                yield ary[i_beg:i+1]
                
            i_beg = i

def cons(ary):
    """Returns all pairs of consecutive items in `ary`"""

    return zip(ary[:-1], ary[1:])

def parse_osm_file(filename, way_filter=None):
    nodes = {}
    ways = {}

    referenced_nodes = set()
    for _, entity in enumerate( parse_file(filename) ):
        
        if isinstance( entity, Node ):
            nodes[ entity.id ] = (entity.lon, entity.lat)
        elif isinstance( entity, Way ):
            if way_filter and not way_filter(entity):
                continue
                
            ways[ entity.id ] = entity
            referenced_nodes.update( entity.nodes )

    # filter nodes to those referenced by a way
    nodes = {ndid:node for ndid,node in nodes.items() if ndid in referenced_nodes}

    return nodes, ways

def is_oneway(way):
    return way.tags.get("oneway") in {"yes","true","1"}

def find_vertex_nodes(highways):
    endpoint_nodes = set()

    vcount = Counter()

    for highway in highways:
        endpoint_nodes.add( highway.nodes[0] )
        endpoint_nodes.add( highway.nodes[-1] )

        # keep a running total of how many times nodes appear
        vcount.update( highway.nodes )

    # nodes that appear more than once are intersection nodes
    intersection_nodes = {vid for (vid, ct) in vcount.items() if ct>1}

    vertex_nodes = intersection_nodes | endpoint_nodes

    return vertex_nodes

def geo_len(pts):
    pair_len = lambda pt0, pt1: geod.inv( pt0[0], pt0[1], pt1[0], pt1[1] )[2]
    return sum([pair_len(x,y) for (x,y) in cons(pts)])

def get_edges(nodes, highways):
    edges = []

    vertex_nodes = find_vertex_nodes( highways.values() )

    for highway in highways.values():
        
        for i, seg in enumerate( chop(highway.nodes, vertex_nodes) ):
            pts = [nodes[nd] for nd in seg]
            
            fromv = seg[0]
            tov = seg[-1]
            seglen = geo_len(pts)
            edge_id_forward = (highway.id,i+1)
            edge_id_backward = (highway.id,-(i+1))
            
            edges.append( (fromv, tov, (edge_id_forward, seglen) ) )
            if not is_oneway(highway):
                edges.append( (tov, fromv, (edge_id_backward, seglen) ) )

    return edges