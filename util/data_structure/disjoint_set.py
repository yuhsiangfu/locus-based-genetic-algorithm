"""
Disjoint set data structure with rank & path compression

@auth: Yu-Hsiang Fu
@date: 2014/12/11
@update: 2018/03/26
"""


class disjoint_set:
    # initialization
    def __init__(self, N):
        self._parent = list(range(N))
        self._rank   = [0] * N


    def find_set(self, x):
        # find root
        if x != self._parent[x]:
            self._parent[x] = self.find_set(self._parent[x])

        return self._parent[x]


    def union(self, x, y):
        # find roots of x and y
        ri = self.find_set(x)
        rj = self.find_set(y)

        # if it's the same root, then return
        if ri == rj:
            return

        # rank
        if self._rank[ri] > self._rank[rj]:
            self._parent[rj] = ri
        else:
            self._parent[ri] = rj
            if self._rank[ri] == self._rank[rj]:
                self._rank[rj] = self._rank[rj] + 1
