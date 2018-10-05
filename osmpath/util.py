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