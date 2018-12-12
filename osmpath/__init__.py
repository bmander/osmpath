"""osmpath - Shortest paths using OpenStreetMap."""

from osmread import parse_file, Way, Node
from collections import Counter
import pyproj
from .util import cons, chop
import numpy as np
from scipy.sparse import csr_matrix, dok_matrix
from scipy.sparse.csgraph import dijkstra
import json

__version__ = '0.1.0'
__author__ = 'Brandon Martin-Anderson <badhill@gmail.com>'
__all__ = []

geod = pyproj.Geod(ellps='clrk66')

def is_oneway(way):
    return way.tags.get("oneway") in {"yes","true","1"}

def geo_len(pts):
    pair_len = lambda pt0, pt1: geod.inv( pt0[0], pt0[1], pt1[0], pt1[1] )[2]
    return sum([pair_len(x,y) for (x,y) in cons(pts)])

class ShortestPaths:
    def __init__(self, graph, origins, preds):
        self.graph = graph
        self.preds = preds

        self.orig_ix = dict(zip(origins, range(len(origins))))

    def get_path(self, orig, dest):
        i = self.orig_ix[orig]
        spt = self.preds[i]

        path = [self.graph.vid_ix[dest]]

        # if the previous-vertex entry for the destination vertex
        # is negative, it means it wasn't reached
        if spt[ path[-1] ]<0:
            return None

        while True:
            last = spt[ path[-1] ]
            if last<0:
                break
            path.append( last )
        
        return [self.graph.ix_vid[x] for x in reversed(path)]

class Graph:
    def __init__(self, edges):
        """
        Args:
            edges (list): List of tuples, with format `(from_vertex, to_vertex, 
                weight)`. `from_vertex` and `to_vertex` can be any hashable
                object. `weight` is a number.
        """

        # transpose edge list into columns
        orig = [edge[0] for edge in edges]
        dest = [edge[1] for edge in edges]
        weights = [edge[2] for edge in edges]

        # convert vertex ids into matrix indices
        # sorted for the sake of consistency
        # it imposes some cost, but the graph should be small enough such that
        # O(n*log n) is easy, so this should be fine.
        vertex_ids = sorted(list( set(orig) | set(dest) ))

        self.vid_ix = dict(zip(vertex_ids, range(len(vertex_ids))))
        self.ix_vid = {v:k for k,v in self.vid_ix.items()}

        ii = [self.vid_ix[vid] for vid in orig]
        jj = [self.vid_ix[vid] for vid in dest]

        # build matrix
        n = len(vertex_ids)
        mat = dok_matrix((n,n))
        mat[ii,jj] = weights

        # compress before returning
        self.mat = mat.tocsr()

    def get_shortest_paths(self, origins):
        """
        Args:
            origins (list): List of vertex ids.
        Returns:
            shortest path tree object
        """

        orig_ix = [self.vid_ix[x] for x in origins]
        _, preds = dijkstra( self.mat, indices=orig_ix, return_predecessors=True)

        return ShortestPaths(self, origins, preds)

class Fiz:
    def __init__(self, edges):
        # set of all node ids
        nds = set([x[0] for x in edges]) | set([x[1] for x in edges])

        # node id -> index lookup
        i_nd = dict( zip(nds, range(len(nds))) )
        # index -> node id lookup
        nd_i = {v:k for k,v in i_nd.items()}

        # distance matrix
        dist = csr_matrix((len(nds),len(nds)))
        i = np.vectorize(i_nd.get)( [x[0] for x in edges] )
        j = np.vectorize(i_nd.get)( [x[1] for x in edges] )
        weight = [edge[1] for _, _, edge in edges]
        dist[i,j] = weight

        # edge details directory
        edge_details = {}
        for fromv, tov, edge in edges:
            edge_details[(fromv,tov)] = edge

        self.i_nd = i_nd
        self.nd_i = nd_i
        self.dist = dist
        self.edges = edge_details

class OSMPathPlanner:
    def __init__(self):
        self.nodes = {}
        self.ways = {}
        self.vertex_nodes = set()

    @classmethod
    def from_osm(cls, filename, way_filter=None, verbose=False):
        ret = cls()

        referenced_nodes = Counter()

        # while we're parsing the file, it's easy to keep track of the OSM
        # nodes that correspond to graph vertices
        ret.vertex_nodes = set()

        for i, entity in enumerate( parse_file(filename) ):
            if verbose:
                if i%100000==0:
                    print( "{} entities read".format(i) )
            
            if isinstance( entity, Node ):
                ret.nodes[ entity.id ] = (entity.lon, entity.lat)
            elif isinstance( entity, Way ):
                if way_filter and not way_filter(entity):
                    continue
                    
                ret.ways[ entity.id ] = entity
                referenced_nodes.update( entity.nodes )

                ret.vertex_nodes.add( entity.nodes[0] )
                ret.vertex_nodes.add( entity.nodes[-1] )

        # filter nodes to those referenced by a way
        ret.nodes = {ndid:node for ndid,node in ret.nodes.items() if ndid in referenced_nodes}

        ret.vertex_nodes.update( [nd for nd,ct in referenced_nodes.items() if ct>1] )

        return ret

    def _get_edges(self):
        """Returns: edge tuples with format `(vertex1, vertex2, edge)`. Each vertex
        is an OSM node id. The `edge` object is a tuple with format `(edge_spec, dist)`,
        where `dist` is the edge distance in meters. `edge_spec` is unique tuple,
        with the format `(way_id, (index0, index1))` where way_id is a valid
        OSM way ID, index0 is the index of the first node in the edge, and
        index1 is the **inclusive** index of the last node in the edge."""

        edges = []

        for highway in self.ways.values():
            
            for j0, j1 in chop(highway.nodes, self.vertex_nodes):
                nds = highway.nodes[j0:j1+1]
                pts = [self.nodes[nd] for nd in nds]
                
                fromv = highway.nodes[j0]
                tov = highway.nodes[j1]
                seglen = geo_len(pts)
                edge_id_forward = (highway.id,(j0,j1))
                edge_id_backward = (highway.id,(j1,j0))
                
                edges.append( (fromv, tov, (edge_id_forward, seglen) ) )
                if not is_oneway(highway):
                    edges.append( (tov, fromv, (edge_id_backward, seglen) ) )

        return edges

    def serialize(self, fn):
        edges = self._get_edges()

        geoms = {}
        for way in self.ways.values():
            pts = np.array([self.nodes[x] for x in way.nodes])
            pts = (pts*1e7).astype(int)
            
            firstrow = pts[0:1,:]
            diffs = np.diff( pts, axis=0 )
            geoms[way.id] = np.vstack( [firstrow, diffs] ).tolist()

        vertex_coords = [(k,v) for k,v in self.nodes.items() if k in self.vertex_nodes]

        with open(fn,"w") as fp:
            json.dump({"edges":edges, "geoms":geoms, 'nodes':vertex_coords}, fp)

    @classmethod
    def deserialize(self, fn):
        with open(fn) as fp:
            data = json.load(fp)
            graph = Graph(data["edges"])
            geoms = data["geoms"]
            geoms = {k:np.array(v).cumsum(axis=0)/1e7 for k,v in geoms.items()}
            nodes = data["nodes"]

            return geoms, nodes, graph

    def get_graph(self):
        edges = self._get_edges()
        return Graph(edges)


if __name__=='__main__':
    gg = Graph([("a", "b", 1), ("b", "c", 2), ("a", "c", 2)])
    spts = gg.get_shortest_paths(['a','b'])
    spts.get_path("b","a")
