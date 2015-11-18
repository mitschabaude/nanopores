'''
python class to handle n-dimensional rectangular boxes and unions of them

note: very naive implementation of union, suitable for a couple 100 boxes
'''

# TODO: collaps pathological boxes,
#       reduce number of boxes in disjoint union when there is no intersection

import itertools

class Box(object):
    
    def __init__(self, a, b):
        # a, b should be tuples
        assert len(a) == len(b)
        self.dim = len(a)
        a, b = _sort(a, b)
        self.a = Point(a)
        self.b = Point(b)
        self.intervals = zip(self.a, self.b)

    def __str__(self):
        return "Box(%s)" %("x".join(["[%s, %s]" %i for i in self.intervals]),)
        
    def __cmp__(self, other):
        return cmp((self.a, self.b), (other.a, other.b))
        
    def __add__(self, other):
        # return union of two Boxes as a BoxCollection
        assert isinstance(other, Box)
        return box_union(self, other)
                        
        
      
class BoxCollection(list):
    " list of disjoint boxes "

    def __str__(self):
        return "Union{ %s }" %(", ".join(str(box) for box in self),)
        
    def __add__(self, other):
        # box + boxcollection
        new = BoxCollection()
        if isinstance(other, Box):
            for box in self:
                new.extend([b for b in box + other if not b in new])
            return new
            
        if isinstance(other, BoxCollection):
            for box in other:
                new.extend([b for b in self + box if not b in new])
            return new
        
        
class Point(object):
    " n-dimensional point "
    tol = 1e-10
    
    def __init__(self, a):
        self.tup = tuple(Float(x) for x in a)
    
    def __cmp__(self, other):
        " test for equality up to tolerance "
        tol = self.tol
        for s, t in zip(self.tup, other.tup):
            c = cmp(s, t)
            if c != 0:
                return c
        return 0
        
    def __eq__(self, other):
        return cmp(self, other) == 0

    def __ne__(self, other):
        return cmp(self, other) != 0
        
    def __len__(self):
        return len(self.tup)
        
    def __getitem__(self, i):
        return self.tup[i]
        
    def __iter__(self):
        return iter(self.tup)
        

class Float(float):
    tol = 1e-10
    
    def __cmp__(self, other):
        return _cmp(self, other, self.tol)
        
    def __eq__(self, other):
        return cmp(self, other) == 0

    def __ne__(self, other):
        return cmp(self, other) != 0
    
        
def _sort(a, b):
    a_ = tuple(map(min, zip(a,b)))
    b_ = tuple(map(max, zip(a,b)))
    return a_, b_
    
    
def _cmp(s, t, tol):
    " test floats for approximate equality "
    if s > t+tol:
        return 1
    elif s < t-tol:
        return -1
    else:
        return 0
        
        
def _sorted(I, J):
    return (min(I,J), max(I,J))
    
            
def interval_union(I, J):
    # takes two intervals and returns list of disjoint intervals with the same union,
    # plus information to which of the initial intervals every new interval belongs
    # interval := tuple (a,b) with a <= b
    # code: 0 := I, 1 := J, 2 := both
    i, j = 0, 1
    if I > J:
        i, j = j, i
        I, J = J, I
        
    a,b = I
    c,d = J
    if b <= c:
        return [I, J], [i, j]
    elif a == c and b == d:
        return [I], [2]
    elif a == c:
        return [(a,b), (b,d)], [2, j]
    elif b == d:
        return [(a,c), (c,b)], [i, 2]
    elif b < d:
        return [(a,c), (c,b), (b,d)], [i, 2, j]
    elif d < b:
        return [(a,c), (c,d), (d,b)], [i, 2, i]
        
        
def _admissible(l):
    return not (any(x==0 for x in l) and any(x==1 for x in l))
    
def _box_from_intervals(ints):
    a = tuple(x for x,y in ints)
    b = tuple(y for x,y in ints)
    return Box(a, b)
        
def box_union(A, B):
    assert A.dim == B.dim
    intss = []
    infos = []
    
    # get disjoint intervals in every dimension
    for i in range(A.dim):
        I = A.intervals[i]
        J = B.intervals[i]
        ints, info = interval_union(I, J)
        intss.append(ints)
        infos.append(info)
        
    # form all possible products of intervals that correspond to a subbox
    # = all combinations with info (0,0,0), (1,1,1), (0,0,2), (1,2,1) etc.
    collection = BoxCollection()
    
    for ints, info in zip(itertools.product(*intss), itertools.product(*infos)):
        #print ints, info
        if _admissible(info):
            box = _box_from_intervals(ints)
            collection.append(box)
            
    return collection
    
    
