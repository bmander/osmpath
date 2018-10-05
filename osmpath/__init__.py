"""osmpath - Shortest paths using OpenStreetMap."""

from osmread import parse_file, Way, Node
from collections import Counter
import pyproj
from .util import cons, chop

__version__ = '0.1.0'
__author__ = 'Brandon Martin-Anderson <badhill@gmail.com>'
__all__ = []

geod = pyproj.Geod(ellps='clrk66')

def is_oneway(way):
    return way.tags.get("oneway") in {"yes","true","1"}

def geo_len(pts):
    pair_len = lambda pt0, pt1: geod.inv( pt0[0], pt0[1], pt1[0], pt1[1] )[2]
    return sum([pair_len(x,y) for (x,y) in cons(pts)])

class OSM:
    def __init__(self):
        self.nodes = {}
        self.ways = {}

    @classmethod
    def _find_vertex_nodes(cls, highways):
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

    @classmethod
    def parse(cls, filename, way_filter=None):
        ret = cls()

        referenced_nodes = set()
        for _, entity in enumerate( parse_file(filename) ):
            
            if isinstance( entity, Node ):
                ret.nodes[ entity.id ] = (entity.lon, entity.lat)
            elif isinstance( entity, Way ):
                if way_filter and not way_filter(entity):
                    continue
                    
                ret.ways[ entity.id ] = entity
                referenced_nodes.update( entity.nodes )

        # filter nodes to those referenced by a way
        ret.nodes = {ndid:node for ndid,node in ret.nodes.items() if ndid in referenced_nodes}

        return ret

    def get_edges(self):
        edges = []

        vertex_nodes = self._find_vertex_nodes( self.ways.values() )

        for highway in self.ways.values():
            
            for i, seg in enumerate( chop(highway.nodes, vertex_nodes) ):
                pts = [self.nodes[nd] for nd in seg]
                
                fromv = seg[0]
                tov = seg[-1]
                seglen = geo_len(pts)
                edge_id_forward = (highway.id,i+1)
                edge_id_backward = (highway.id,-(i+1))
                
                edges.append( (fromv, tov, (edge_id_forward, seglen) ) )
                if not is_oneway(highway):
                    edges.append( (tov, fromv, (edge_id_backward, seglen) ) )

        return edges



