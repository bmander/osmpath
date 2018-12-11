def chop(ary, endpoints):
    """Returns inclusive slice indices of segments of `ary` bounded by values `endpoints`.

    For example: `chop([1,2,3,4,5,6,7], [2,5,7])` -> `[(1,4),(4,6)]`
    """

    # ary indices that are endpoints
    ix = [i for i,item in enumerate(ary) if item in endpoints]

    return [(i0, i1) for i0, i1 in cons(ix)]

def cons(ary):
    """Returns all pairs of consecutive items in `ary`"""

    return zip(ary[:-1], ary[1:])