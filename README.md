osmpath
===========

Shortest paths using OpenStreetMap. In-process; efficient for batch queries. Ideal for city-scale data analysis.

```python
>>> from osmpath import OSMPathPlanner
>>> def is_bikeable(way):
...   nogo = {'motorway', 'motorway_link', 'footway', 'service'}
...   return ("highway" in way.tags) and (way.tags["highway"] not in nogo)
>>> osm = OSMPathPlanner.from_osm("path/to/sf.osm", way_filter=is_bikeable, verbose=True)
>>> # specify more than one origin/destination pair
>>> orig = [[-122.45893333,   37.80118333],
...         [-122.41249   ,   37.79199167],
...         [-122.41979667,   37.77634   ]]
>>> dest = [[-122.41249667,   37.791975  ],
...         [-122.41979667,   37.77634   ],
...         [-122.41767833,   37.77632   ]]
>>> paths = list( osm.get_shortest_paths( orig, dest ) )
>>> len(paths)
3
>>> paths[0] # each edge in a path is an OSM way_id and an inclusive range of node indices.
[[219171857, [0, 6]],
 [219171857, [6, 7]],
 [529868660, [0, 1]],
...
 
>>> osm.get_path_geom( paths[0] )
array([[-122.459098 ,   37.8018337],
       [-122.459182 ,   37.8018022],
       [-122.4592759,   37.801769 ],
       [-122.4593613,   37.8017409],
...


```

Authors
-------

`osmpath` was written by `Brandon Martin-Anderson <badhill@gmail.com>`_.
