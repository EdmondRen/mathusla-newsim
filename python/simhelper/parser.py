import numpy as np

def unpack_cov(cov: list, dim: int = 3) -> np.ndarray:
    ncov = len(cov)//(dim*dim)
    if (len(cov)%(dim*dim)!=0) :
        warning.warn("Check the dimension!")

    return np.reshape(cov, (ncov, dim,dim))


def unpack_at(concatlist, divider=-1):
    lists = []
    n = 0
    for val in concatlist:
        if val == -1:
            n += 1
        else:
            while len(lists) <= n:
                lists.append([])
            lists[n].append(val)
    return lists