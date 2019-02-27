"""osmpath - Shortest paths using OpenStreetMap."""

from osmread import parse_file, Way, Node
from collections import Counter
import pyproj
from .util import cons, chop
import numpy as np
from scipy.sparse import csr_matrix, dok_matrix
from scipy.sparse.csgraph import dijkstra, connected_components
import json
from scipy.spatial import KDTree

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

def get_largest_connected_component(g):
    # Find all connected components. We're interested in a strongly connected
    # component, because in order to consider a place "in the same street grid"
    # you need to be able to both get there and get away from there.
    _, comps = connected_components( g.mat, connection="strong" )

    # find the connected component with the most vertices
    big_comp = np.bincount(comps).argmax()
    
    # get index of mostly common connected component, and map to new indices
    ix = np.nonzero( comps==big_comp )[0]
    oldix_newix = dict(zip(ix, range(len(ix))))

    # cut down the graph to just the interesting indices, and update the
    # vertex_ix -> index dictionaries.
    ret = Graph([])
    ret.mat = g.mat[:,ix][ix,:]
    ret.ix_vid = {oldix_newix[k]:v for k,v in g.ix_vid.items() if k in oldix_newix}
    ret.vid_ix = {v:k for k,v in ret.ix_vid.items()}

    return ret

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

    def get_shortest_paths(self, origins, verbose=False, chunk=30):
        """
        Args:
            origins (list): List of vertex ids.
        Returns:
            shortest path tree object
        """

        orig_ix = [self.vid_ix[x] for x in origins]

        # find {chunk} shortest path trees at a time
        predchunks = []
        for i in range(0,len(orig_ix),chunk):
            if verbose==True:
                print( "finding SPTs %d-%d"%(i,i+chunk) )

            _, preds = dijkstra( self.mat, indices=orig_ix[i:i+chunk], 
                                 return_predecessors=True)
            predchunks.append( preds )
        preds = np.vstack( predchunks )

        return ShortestPaths(self, origins, preds)

class SpatialIndex:
    def __init__(self, ids, points):
        self.ix_id = dict(zip( range(len(ids)), ids ))
        self.tree = KDTree(points)

    def query(self, points):
        dists, ix = self.tree.query( points )
        return dists, [self.ix_id[x] for x in ix]

class OSMPathPlanner:

    @classmethod
    def _parse_osm(cls, filename, way_filter=None, verbose=False):
        referenced_nodes = Counter()

        # while we're parsing the file, it's easy to keep track of the OSM
        # nodes that correspond to graph vertices
        vertex_nodes = set()
        nodes = {}
        ways = {}

        for i, entity in enumerate( parse_file(filename) ):
            if verbose:
                if i%100000==0:
                    print( "{} entities read".format(i) )
            
            if isinstance( entity, Node ):
                nodes[ entity.id ] = (entity.lon, entity.lat)
            elif isinstance( entity, Way ):
                if way_filter and not way_filter(entity):
                    continue
                    
                ways[ entity.id ] = entity
                referenced_nodes.update( entity.nodes )

                vertex_nodes.add( entity.nodes[0] )
                vertex_nodes.add( entity.nodes[-1] )

        # filter nodes to those referenced by a way
        nodes = {ndid:node for ndid,node in nodes.items() if ndid in referenced_nodes}

        vertex_nodes.update( [nd for nd,ct in referenced_nodes.items() if ct>1] )

        return ways, nodes, vertex_nodes

    @classmethod
    def _get_edges(cls, ways, nodes, vertex_nodes):
        """Returns: edge tuples with format `(vertex1, vertex2, edge)`. Each vertex
        is an OSM node id. The `edge` object is a tuple with format `(edge_spec, dist)`,
        where `dist` is the edge distance in meters. `edge_spec` is unique tuple,
        with the format `(way_id, (index0, index1))` where way_id is a valid
        OSM way ID, index0 is the index of the first node in the edge, and
        index1 is the **inclusive** index of the last node in the edge."""

        edges = []

        for highway in ways.values():
            
            for j0, j1 in chop(highway.nodes, vertex_nodes):
                nds = highway.nodes[j0:j1+1]
                pts = [nodes[nd] for nd in nds]
                
                fromv = highway.nodes[j0]
                tov = highway.nodes[j1]
                seglen = geo_len(pts)
                edge_id_forward = (highway.id,(j0,j1))
                edge_id_backward = (highway.id,(j1,j0))
                
                edges.append( (fromv, tov, (edge_id_forward, seglen) ) )
                if not is_oneway(highway):
                    edges.append( (tov, fromv, (edge_id_backward, seglen) ) )

        return edges

    @classmethod
    def from_osm(cls, filename, way_filter=None, verbose=False):
        ways, nodes, vertex_nodes = cls._parse_osm(filename, way_filter, verbose)

        edges = cls._get_edges(ways, nodes, vertex_nodes)

        return cls(nodes, ways, edges)

    def __init__(self, nodes, ways, edges):
        # node keys are integers, but json keys are strings; ungarble them
        self.nodes = {int(k):v for k,v in nodes.items()}
        self.ways = {int(k):Way(*v) for k,v in ways.items()}
        self.edges = edges

        self.edge_index = {(a,b):c for a,b,(c,_) in edges}

        edge_weights = [(a, b, c) for (a, b, (_, c)) in edges]
        self.graph = Graph(edge_weights)

        # only spatially index vertex nodes
        vids = set([x[0] for x in edge_weights]) | set([x[1] for x in edge_weights])
        vertex_nodes = {k:v for k,v in self.nodes.items() if k in vids}
        self.index = SpatialIndex(list(vertex_nodes.keys()), 
                                  list(vertex_nodes.values()))
        
    def serialize(self, fn):
        with open(fn,"w") as fp:
            json.dump({"edges":self.edges, 
                       "nodes":self.nodes, 
                       'ways':self.ways}, fp, indent=2)

    @classmethod
    def deserialize(cls, fn):
        with open(fn) as fp:
            data = json.load(fp)
            nodes = data["nodes"]
            ways = data["ways"]
            edges = data["edges"]

            ret = cls(nodes, ways, edges)
            return ret

    def _vertex_pairs_to_edges(self, vertex_path):
        ret = []
        for vpair in cons(vertex_path):
            ret.append( self.edge_index[vpair] )
        return ret

    def get_shortest_paths(self, orig, dest, verbose=False):
        # convert locations to vertex ids
        _, orig_vids = self.index.query( orig )
        _, dest_vids = self.index.query( dest )

        assert len(orig_vids) == len(dest_vids)

        unique_orig_vids = list(set(orig_vids))
        spts = self.graph.get_shortest_paths(unique_orig_vids, verbose)

        for o, d in zip(orig_vids, dest_vids):
            vertex_path = spts.get_path(o,d)
            if vertex_path is None:
                yield None
            else:
                yield self._vertex_pairs_to_edges( vertex_path )

    def get_edge_geom( self, edge ):
        way_id, (ix1, ix2) = edge

        assert ix1 != ix2

        # nd_ix is the inclusive indices of the segment to slice from the
        # way.
        if ix2 > ix1:
            segslice = slice(ix1, ix2+1)
        else:
            if ix2==0:
                segslice = slice(ix1, None, -1)
            else:
                segslice = slice(ix1, ix2-1, -1)
        
        nds = self.ways[way_id].nodes[segslice]
        return np.array( [self.nodes[nd] for nd in nds] )

    def get_path_geom( self, path ):
        return np.vstack( [self.get_edge_geom( x ) for x in path] )


if __name__=='__main__':
    gg = Graph([("a", "b", 1), ("b", "c", 2), ("a", "c", 2)])
    spts = gg.get_shortest_paths(['a','b'])
    spts.get_path("b","a")
