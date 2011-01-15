from __future__ import division

import math
import random
import fractions


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

def lcm(a, b):
    return (a * b) / gcd(a, b)

def _trial_div(n): ## From http://blog.dreamshire.com/2009/03/26/94/
    if n == 1: return 1
    for p in (2, 3, 5):
        if n % p == 0: return p
    dif = [6, 4, 2, 4, 2, 4, 6, 2]
    m = 7; i = 1
    while m <= n and m*m <= n:
        if n % m == 0:
            return m
        m += dif[i % 8]
        i += 1
    return n
 
def factor(n): ## From http://blog.dreamshire.com/2009/03/26/94/
    if n in [-1, 0, 1]: return
    if n < 0: n = -n
    while n != 1:
        p = _trial_div(n)
        e = 1
        n /= p
        while n % p == 0:
            e += 1; n /= p
        yield p, e

def n_root(n, r=2, c=100, x_k=None):
    if x_k is None:
        x_k = int(n ** (1 / r))
    for _ in xrange(c):
        x_k1 = ((r - 1) * x_k + n / (x_k ** (r - 1))) / r
        x_k = x_k1
    return x_k

def random_root(base_range=xrange(1, 101), index_range=xrange(2, 4), coeff_range=xrange(-1, 10), not_null=False):
    kwargs = vars()
    b, i, c = map(random.choice, (base_range, index_range, coeff_range))
    if b == 0 and not_null:
        return random_root(**kwargs)
    return Root(b, i, c)

class Root(object):
    def __init__(self, b, index=2, coeff=1):
        self._b = b
        self._index = index
        self._coeff = coeff
        if isinstance(b, Root):
            if b._coeff != 1:
                b = b.move_in()
            self._b = b._b
            self._index *= b._index

    @ property
    def b(self):
        return self._b

    @ b.setter
    def b(self, v):
        self._b = v

    @ property
    def index(self):
        return self._index

    @ index.setter
    def index(self, v):
        self._index = v

    @ property
    def coeff(self):
        return self._coeff

    @ coeff.setter
    def coeff(self, v):
        self._coeff = v

    @ property
    def root(self):
        if isinstance(self._b, fractions.Fraction):
            return fractions.Fraction.from_float(n_root(self._b.numerator, self._index)) / \
                   fractions.Fraction.from_float(n_root(self._b.denominator, self._index))

        return n_root(self._b, self._index)

    def copy(self):
        return Root(self._b, self._index, self._coeff)

    def move_in(self):
        if self._coeff == 1:
            return self.copy()
        i = self._coeff ** self._index
        return Root(i * self._b, self._index)

    def move_out(self):
        if self._coeff != 1:
            return self.move_in().move_out()
        coeff, b = 1, 1
        for f, e in factor(self._b):
            if e < self._index:
                b *= f
                continue
            q, r = divmod(e, self._index)
            coeff *= f * q
            b *= f ** r
        if b == 1 and self._b != 1:
            return self.copy()
        return Root(b, self._index, coeff)

    def same_index(self, other):
        i = lcm(self._index, other._index)
        return Root(i / self._b, i), Root(i / other._b, i)

    def similar(self, other):
        return self._b == other._b and \
               self._index == other._index

    def __repr__(self):
        if self._coeff == 1:
            return 'Root(b={0}, index={1})'.format(self._b, self._index)
        return 'Root(b={0}, index={1}, coeff={2})'.format(self._b, self._index,
                                                          self._coeff)

    def __add__(self, other):
        a, b = self.move_out(), other.move_out()
        if a.similar(b):
            return Root(a._b, a._index, a._coeff + b._coeff)
        return NotImplemented

    def __sub__(self, other):
        a, b = self.move_out(), other.move_out()
        if a.similar(b):
            return Root(a._b, a._index, a._coeff - b._coeff)
        return NotImplemented

    def __mul__(self, other):
        r1, r2 = self.same_index(other)
        return Root(r1._b * r2._b, r1._index)

    def __div__(self, other):
        r1, r2 = self.same_index(other)
        return Root(r1._b / r2._b, r1._index)

    def __pow__(self, exp):
        return Root(self._b ** exp, self._index)

    def __eq__(self, other):
        return self._b == other._b and \
               self._index == other._index and \
               self._coeff == other._coeff

    def __ne__(self, other):
        return not self == other