#----------------------------------------------------------------------------------------------------#

cdef inline int cy_healpy_lm2ind(int el, int em, int L):
        return em*(2*L-1-em)/2+el

cdef inline int cy_healpy_ind2lm(int ind, int L):

    cdef int n_elements = L, m=0, el

    while (True):
        if ind >= n_elements:
            ind -= n_elements
            n_elements -= <int>1
            m += 1
        else:
            el = ind + m
            return (el, m)

cdef inline int cy_mw_elm2ind( int el, int m):

  return el * el + el + m 