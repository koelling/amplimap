from cpython.mem cimport PyMem_Malloc, PyMem_Free

#check in jupyter with: %%cython --annotate

cdef unsigned int mismatch_count_with_max_cy(unsigned int my_len, char* a, char* b, unsigned int max_mismatches = 1):
    """Count number of mismatches between two strings, but stop once max_mismatches has been reached and just return max_mismatches+1."""
    cdef unsigned int found_mm = 0

    cdef unsigned int i
    for i in range(my_len):
        if a[i] != b[i] and a[i] != b'.' and b[i] != b'.':
            found_mm += 1
            if found_mm > max_mismatches:
                return found_mm

    return found_mm

#a is full string, b is prefix
def mismatch_count_with_max(a, b, max_mismatches = 1):
    """Python wrapper for mismatch_count_with_max_cy."""
    comparable_len = min(len(a), len(b))
    found_mm = mismatch_count_with_max_cy(comparable_len, a[:comparable_len], b[:comparable_len], max_mismatches)

    #if len(a) > len(b) [normal case, where b is prefix of a] we just ignore the rest
    #if len(b) > len(a) we add to mismatch
    if len(b) > len(a):
        found_mm += len(b) - len(a)

    return found_mm

cdef list find_prefixes_with_mismatches_remove(char* str, tuple p, unsigned int k):
    """Search for overlap between a single string and a list (tuple) of prefix strings at the same time.

    Args:
        str: The string to search for (eg. read).
        p: Tuple of prefix strings (eg. expected primer arms).
        k: Maximum number of mismatches to allow.

    Returns:
        list of tuples (prefix index, number of mismatches)
    """
    cdef char str_char
    cdef char prefix_char
    cdef unsigned int char_ix
    cdef unsigned int prefix_ix
    cdef unsigned int n = len(p)
    cdef unsigned int max_search_len = 0
    cdef unsigned int str_len = len(str)
    cdef bytes prefix_str
        
    result = []
    
    #TODO: these should really be allocated outside the main loop already! at the moment we realloc this for every read...
    cdef char **prefixes = <char **>PyMem_Malloc(n * sizeof(char *))
    cdef unsigned int *lengths = <unsigned int *>PyMem_Malloc(n * sizeof(unsigned int))
    cdef unsigned int *mismatches = <unsigned int *>PyMem_Malloc(n * sizeof(unsigned int))
    if not prefixes or not mismatches or not lengths:
        raise MemoryError()
    try:
        for prefix_ix in range(n):
            prefix_str = p[prefix_ix]
            prefixes[prefix_ix] = prefix_str
            lengths[prefix_ix] = len(prefix_str)
            mismatches[prefix_ix] = 0
            
            if lengths[prefix_ix] > max_search_len:
                max_search_len = lengths[prefix_ix]
                
        if max_search_len > str_len:
            max_search_len = str_len

        #print('max len =', max_search_len)
        for char_ix in range(max_search_len):
            str_char = str[char_ix]
            #print('#', char_ix, ':', str_char)

            for prefix_ix in range(n):
                if mismatches[prefix_ix] <= k:
                    #have we reached the end?
                    if lengths[prefix_ix] <= char_ix:
                        #then this as a match
                        #print('found: ', prefix_ix, 'with mm =', mismatches[prefix_ix], 'after', char_ix, 'chars')

                        result.append((prefix_ix, mismatches[prefix_ix]))
                        mismatches[prefix_ix] = k+1
                    else:
                        prefix_char = prefixes[prefix_ix][char_ix]
                        #print(prefix_char)
                        #do we have a mismatch
                        if str_char != prefix_char and str_char != b'.' and prefix_char != b'.':
                            mismatches[prefix_ix] += 1

        #do one last loop through prefixes to find full-length cases (where lengths[prefix_ix] >= max_search_len)
        for prefix_ix in range(n):
            #cases where search len < prefix len are never a match, since search string was too short
            #NOTE: may want to change this by instead adding length diff to mismatch count, but that would be inconsistent with how we handle R2
            if lengths[prefix_ix] <= max_search_len:
                if mismatches[prefix_ix] <= k:            
                    result.append((prefix_ix, mismatches[prefix_ix]))
                        
    finally:
        # return the previously allocated memory to the system
        PyMem_Free(prefixes)
        PyMem_Free(lengths)
        PyMem_Free(mismatches)
                
    return result

def find_probe_from_reads(
        tuple first_arms_probes,
        tuple first_arms_seq,
        dict second_arms,
        char* first_read_str, char* second_read_str,
        unsigned int umi_one, unsigned int umi_two,
        unsigned int mismatches,
        debug = False):
    """
    Find probes for which the primer arms match the starts of the first and second read.
    """
    cdef unsigned int match_ix, probe_ix, mismatches_first
    cdef char* second_prefix
    
    matches = []
    matches_first = find_prefixes_with_mismatches_remove(first_read_str[umi_one:], first_arms_seq, mismatches)    
    if debug:
        print('R1:', first_read_str[umi_one:])
        for fas in first_arms_seq:
            print(fas)
        print('Matches for first arm:', matches_first)

    for match_ix in range(len(matches_first)):
        probe_ix = matches_first[match_ix][0]
        mismatches_first = matches_first[match_ix][1]
        probe = first_arms_probes[probe_ix]
        second_prefix = second_arms[probe]

        if debug:
            print('R2:', second_read_str[umi_two:], 'len =', len(second_read_str) - umi_two, len(second_prefix))
            print('Testing second arm:', probe, second_prefix)

        if len(second_read_str) - umi_two >= len(second_prefix):
            mismatches_second = mismatch_count_with_max_cy(len(second_prefix), second_read_str[umi_two:], second_prefix, mismatches)  

            if debug:
                print('Mismatches in second arm:', mismatches_second)        
            if mismatches_second <= mismatches:
                matches.append([probe, mismatches_first, mismatches_second])
        elif debug:
            print('Too short!')

    return matches