# some basic functions

import numpy as np

def match(ref_ids, new_ids, uniq_ref_only=True):
    """
    Mapping new_ids to ref_ids. ref_ids can have repeated values, but new_ids 
    can only have unique ids or values. Therefore, new_ids[RT_idx] will be 
    the same as ref_ids. Note, 
    
    Parameters
    ----------
    ref_ids : array_like or list
        ids for reference with type of int, float, or string
    new_ids : array_like or list
        ids waiting to map.
        
    Returns
    -------
    RV_idx : array_like, the same length of ref_ids
        The index for new_ids mapped to ref_ids. If an id in ref_ids does not 
        exist in new_ids, then return a None for that id. 
    Examples
    --------
    >>> x1 = [5, 9, 1]
    >>> x2 = [1, 2, 5, 7, 9]
    >>> hilearn.match(x1, x2)
    array([2, 4, 0])
    >>> hilearn.match(x2, x1)
    array([2, None, 0, None, 1], dtype=object)
    >>> RT_idx = hilearn.match(x2, x1)
    >>> idx1 = np.where(RT_idx != None)[0]
    >>> idx1
    array([0, 2, 4])
    >>> idx2 = RT_idx[idx1].astype(int)
    >>> idx2
    array([2, 0, 1])
    """
    idx1 = np.argsort(ref_ids)
    idx2 = np.argsort(new_ids)
    RT_idx1, RT_idx2 = [], []
    
    i, j = 0, 0
    while i < len(idx1):
        if j == len(idx2) or ref_ids[idx1[i]] < new_ids[idx2[j]]:
            RT_idx1.append(idx1[i])
            RT_idx2.append(None)
            i += 1
        elif ref_ids[idx1[i]] == new_ids[idx2[j]]:
            RT_idx1.append(idx1[i])
            RT_idx2.append(idx2[j])
            i += 1
            if uniq_ref_only: j += 1
        elif ref_ids[idx1[i]] > new_ids[idx2[j]]:
            j += 1
            
    origin_idx = np.argsort(RT_idx1)
    RT_idx = np.array(RT_idx2)[origin_idx]
    return RT_idx
