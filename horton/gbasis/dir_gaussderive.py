#!/usr/bin/python

import numpy as np

def dir_gaussian_derivative(r, atom, n, alpha, norm, order=0, derx=0, dery=0, derz=0):
    
    # Evaluate distance between point and nucleus
    poly = np.zeros([4])
    poly[0] = 1.0
    poly[1] = r[0] - atom[0]
    poly[2] = r[1] - atom[1]
    poly[3] = r[2] - atom[2]

    # Evaluate exponential part 
    pre0 = np.exp(-alpha*(poly[1]**2 + poly[2]**2 + poly[3]**2))*norm
 
    # Evaluate exponential part for +1 order in derivative
    pre0_h = -2.0*alpha
    
    # Size of array needed
    workdim = np.zeros([order + 1], dtype = int)
    #elem = np.zeros([order + 1], dtype = int)
    workdim[0] = 1  
    #elem[0] = 1  
    #print "wasabi order", order
    for iorder in range(1,order+1):
        workdim[iorder] += workdim[iorder-1] + 2**iorder 
        #elem[i] = 6**i 
        print "workdim[", iorder, "]=", workdim[iorder]
        #print "elem[", i, "]=", elem[i]
   
    # Evaluate gaussian function
    tmp = poly[1]**n[0]
    tmp *= poly[2]**n[1]
    tmp *= poly[3]**n[2]
    offset = np.empty([workdim[order]], dtype = int)
    offset[0] = (workdim[order] - 1)/2 # The middle of the tree
    work_cart = np.empty([workdim[order]])
    print "offset[0]=", offset[0]
    print "workdim[", order, "]=", workdim[order]
    work_cart[offset[0]] = pre0*tmp
    #tmp_n = np.copy(n)
 
    # Iterate over derivatives of gaussian function
    idx = 0
    for i in range(0,order):
        print "wasabi i", i+1, " of order", order
        nelem = workdim[order - i - 1]
        #print "wasabi nelem", nelem
        for j in range(0, 2**i):
            print "wasabi j", j+1, " of workdim[i-1]+1", workdim[i]+1
            # Derive toward x
            if derx>i:
                offset[idx+1+i+j] =offset[idx] + (nelem + 1)/2
                print "offset[", idx+1+i+j, "]=", offset[idx+1+i+j]
                offset[idx+2+i+j] =offset[idx] - (nelem + 1)/2
                print "offset[", idx+2+i+j, "]=", offset[idx+2+i+j]
                work_cart[offset[idx+1+i+j]] = pre0_h*poly[1]*work_cart[offset[idx]]
                if (n[0]-i) > 0:
                    work_cart[offset[idx+2+i+j]] = (n[0]-i)/poly[1]*work_cart[offset[idx]]
                else: 
                    work_cart[offset[idx + 2+i+j]] = 0.
            # Derive toward y
            elif dery>(i-derx):
                offset[idx+1+i+j] =offset[idx] + (nelem + 1)/2
                print "offset[", idx+1+i+j, "]=", offset[idx+1+i+j]
                offset[idx+2+i+j] =offset[idx] - (nelem + 1)/2
                print "offset[", idx+2+i+j, "]=", offset[idx+2+i+j]
                work_cart[offset[idx+1+i+j]] = pre0_h*poly[2]*work_cart[offset[idx]]
                if (n[1]-i) > 0:
                    work_cart[offset[idx+2+i+j]] = (n[1]-i)/poly[2]*work_cart[offset[idx]]
                else: 
                    work_cart[offset[idx + 2+i+j]] = 0.
            # Derive toward z
            else:
                offset[idx+1+i+j] =offset[idx] + (nelem + 1)/2
                print "offset[", idx+1+i+j, "]=", offset[idx+1+i+j]
                offset[idx+2+i+j] =offset[idx] - (nelem + 1)/2
                print "offset[", idx+2+i+j, "]=", offset[idx+2+i+j]
                work_cart[offset[idx+1+i+j]] = pre0_h*poly[3]*work_cart[offset[idx]]
                if (n[2]-i) > 0:
                    work_cart[offset[idx+2+i+j]] = (n[2]-i)/poly[3]*work_cart[offset[idx]]
                else: 
                    work_cart[offset[idx+2+i+j]] = 0.

            print "work_cart[",offset[idx],"]=",work_cart[offset[idx]]
            print "work_cart[",offset[idx+1+i+j],"]=",work_cart[offset[idx+1+i+j]]
            print "work_cart[",offset[idx+2+i+j],"]=",work_cart[offset[idx+2+i+j]]
            idx += 1
 
    # Add terms to obtain actual derivatives
    result = 0.
    print"derx=", derx
    print"dery=", dery
    print"derz=", derz
    idx_offset = np.zeros(order)

    if (order==0):
        result = work_cart[offset[0]]

        
    return result

