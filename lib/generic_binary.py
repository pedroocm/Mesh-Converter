# coding: utf-8
from math import floor
def binarySearch(array, value, compFuncs):
    arrayLength = len(array)
    if arrayLength == 1:
        if value == array[0]:
            return 0
        else:
            return -1
    notFound = True
    length = arrayLength
    middle = 0
    initialIndex = 0
    finalIndex = arrayLength - 1
    eq = compFuncs[0]
    lt = compFuncs[1]
    gt = compFuncs[2]
    while notFound:
        #print(middle)
        middle = floor( (initialIndex + finalIndex) / 2)
        #print("Passou")
        if eq(value, array[middle]):
            return middle
        elif lt(value, array[middle]):
            if middle == initialIndex:
                return -1
            else:
                initialIndex = initialIndex
                finalIndex = middle-1
                length = (finalIndex - initialIndex) + 1
        elif gt(value, array[middle]):
            if  middle == finalIndex:
                return -1
            else:
                initialIndex = middle+1
                finalIndex = finalIndex
                length = (finalIndex - initialIndex) + 1
                
    
        #print(length)
