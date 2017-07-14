# -*- coding: utf-8 -*-
"""
Created on Sat Jul 18 14:20:14 2015

@author: yaron
"""

class Collection_base(object):
    '''
    Manages a collection of items. An item is assigned a number when added to 
    the collection, and this number serves as a unique key throughout the 
    lifetime of the collection. Colelction can be given a name on
    initialization.
    
    Supports iteration via 
       for uniqueKey, y in collection_base_object:
           do_something(y)
   Supports access via
      collection_base_object[unique_id]
    '''    
    def __init__(self, name):
        self.name = name
        self.__items = {}
        self.__next_key = 1

    def __len__(self):
        '''
        Number of objects in this this class contains.
        '''
        return len(self.__items)

    def __iter__(self):
        '''
        Allows the following:
        >>> for x, y in mol_set:
        ...    do_something(y)

        '''
        for key, value in self.__items.items():
            yield key, value

    def __getitem__(self, id):
        '''
        allows Collection_base[x] to return the item with the given id
        '''
        return self.__items[id]

    def add(self, x):
        '''
        x is given unique id and added to collection
        
        Returns:
        integer that was assigned as the unique id
        '''
        my_id = self.__next_key
        self.__items[my_id] = x
        self.__next_key = my_id + 1
        return my_id

    def ids(self):
        '''
        Returns a list of all ids currently in the collection
        '''
        return self.__items.keys()
