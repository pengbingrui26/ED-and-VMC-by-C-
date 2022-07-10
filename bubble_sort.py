import numpy as np

def bubble_sort(arr):
    brr = arr.copy()
    for i in range(len(brr)-1, -1, -1):
        for j in range(i):
            if brr[j] > brr[j+1]:
                brr[j], brr[j+1] = brr[j+1], brr[j]
    return brr 

a = [2,1,4,3,9,-1]
b = bubble_sort(a)
print(b)
