from collections import Counter

class Fragment(Counter):
    """
    >>> fragment01 = Fragment('*AB')
    >>> fragment02 = Fragment('*ABC')
    >>> fragment01 < fragment02
    True
    >>> Fragment('*') <= Fragment('*A')
    True

    """
    def __init__(self, key):


        super(Fragment, self).__init__(key)
        self.key = key

    def __hash__(self):

        return hash(self.key)
        
    def __le__(self, other):
        
        for key, value in self.items():
            if key not in other or self[key]> other[key]:
                return False
        return True
        
    def __lt__(self, other):
        
        return self != other and self <= other
        
    def __gt__(self, other):
        
        return other < self and not self == other
    
    def __ge__(self, other):
        
        return other <= self
    
    def difference(self, other):
        
        assert other <= self, 'given fragment must be a supersets'
        
        s = {}
        for key, value in self.items():
            if other[key] != value:
                s[key] = value - other[key]

        return Fragment(s)
        
    def __len__(self):
        
        return sum(self.values())
    
    def __repr__(self):

        result = ''
        for key in sorted(self):
            result += key * self[key]
        return result
    
    def __str__(self):
        
        return repr(self)

if __name__ == '__main__':
    import doctest
    doctest.testmod()