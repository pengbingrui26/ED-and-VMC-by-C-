import numpy as np

def comb(n, m, M): # enumerate all arrays with length m from the array arr whose length is n
     if m<=1:
         res = [ [x] for x in range(n-M+1) ]
         return res
     else:
         res = comb(n, m-1, M)
         res_new = []
         for xx in res:
             for a in range(xx[-1] + 1, n):
                 res_new.append(xx+[a])
         return res_new

#aa = comb(6, 3, 3)
#print(aa)

def comb_new(empty_arr, n, m, M): # enumerate all arrays with length m from the array arr whose length is n
     if m<=1:
         for x in range(n-M+1):
             empty_arr.append([x])
         return empty_arr
     else:
         res = comb_new(empty_arr, n, m-1, M)
         res_new = []
         for xx in res:
             for a in range(xx[-1] + 1, n):
                 empty_arr.append(xx+[a])
         empty_arr

print(comb_new([], 4, 3, 3))


