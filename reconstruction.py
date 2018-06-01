import sys
import pprint
from fragment import Fragment
        
# graph from start to stop
# convert fragments from string to Counter

fragments=['*'] + ['*A', '*AB', '*ABC', '*ABCD', '*CCDDD', '*ABCDE', '*AABEE', '*ACCDDD', '*AABCDE', '*AABCEE', '*ACCCDDD', '*AABCCDE', '*AABBCEE', '*AAAADDE', '*AABBCCDE' ]
fragments=[Fragment(fragment) for index, fragment in enumerate(fragments)]
#fragments=[Fragment(fragment) for index, fragment in enumerate(fragments) if index % 3 in {0, 1}]
    
# Create graph

graph=dict()
for fragment in fragments:
    graph[fragment]=set()
    
for fragment1 in fragments:
    
    for fragment2 in fragments:
        
        # make edge if fragment1 is a subset of fragment2
        if fragment1 < fragment2:
            
            # determine difference between both fragments
            diff=fragment2.difference(fragment1)
            
            if Fragment('*') <= diff:
                diff = diff.difference(Fragment('*'))
                
#             if Fragment('o') <= diff:
#                 diff = diff.difference(Fragment('o'))
            
            # add new edge
            graph[fragment1].add((fragment2, 1/len(diff)))

# find the longest path
def longest_path(initial, graph):
    
    # step 1: assign to every node a tentative distance value: set it to zero for 
    # our initial node and to infinity for all other nodes
    longest_paths = {fragment:(None, None) for fragment in graph}
    current = initial
    longest_paths[current] = (str(current), 0.0) 
    
    # step 2: set the initial node as current. Mark all other nodes unvisited; 
    # create a set of all the unvisited nodes called the unvisited set
    processed = set()
    unprocessed = set(graph)
    
    while current:
        
        # step 3: for the current node, consider all of its unvisited neighbors and 
        # calculate their tentative distances; compare the newly calculated tentative 
        # distance to the current assigned value and assign the larger one; for example, 
        # if the current node A is marked with a distance of 6, and the edge connecting 
        # it with a neighbor B has length 2, then the distance to B (through A) will be 
        # 6 + 2 = 8. If B was previously marked with a distance smaller than 8 then 
        # change it to 8; otherwise, keep the current value
        for neighbour, edge_weight in graph[current]:
            if neighbour in unprocessed:
                new_distance = longest_paths[current][1] + edge_weight
                if (
                    longest_paths[neighbour][1] is None or 
                    new_distance > longest_paths[neighbour][1]
                ):
                    
                    diff = str(neighbour.difference(current))
                    if len(diff) > 1:
                        diff = '[' + diff + ']'
                    longest_paths[neighbour] = (
                        longest_paths[current][0] + diff,
                        new_distance
                    )
                    
        # step 4: when we are done considering all of the neighbors of the current node, 
        # mark the current node as visited and remove it from the unvisited set; a 
        # visited node will never be checked again
        processed.add(current)
        unprocessed.remove(current)
        
        # step 5: if the destination node has been marked visited (when planning a route 
        # between two specific nodes) or if the largest tentative distance among the 
        # nodes in the unvisited set is infinity (when planning a complete traversal; 
        # occurs when there is no connection between the initial node and remaining 
        # unvisited nodes), then stop; the algorithm has finished; otherwise, select the 
        # unvisited node that is marked with the largest tentative distance, set it as 
        # the new "current node", and go back to step 3.

        largest_fragment, largest_distance = None, None
        for fragment in unprocessed:
            if largest_distance is None or longest_paths[fragment][1] > largest_distance:
                largest_fragment, largest_distance = fragment, longest_paths[fragment][1]

        current = largest_fragment if largest_distance is not None else None

    # find the fragment with the largest path
    largest_path, largest_distance = None, None
    for fragment in processed:
        if largest_distance is None or longest_paths[fragment][1] > largest_distance:
            largest_path, largest_distance = longest_paths[fragment]
    return largest_path

print('forward path: {}'.format(longest_path(Fragment('*'), graph)[1:]))

    
#########################################################
#########################################################
#########################################################

# Graph from stop to start
# convert fragments from string to Counter

fragments = ['*'] + ['*B', '*BC', '*ABC', '*ABCE', '*ABCDE', '*AD', '*ADD', '*ADDE', '*AADDE', '*CCCCE', '*AACDDE', '*AAABEE', '*CCCCCE', '*ABCCDE', '*ABBCCDE', '*CCDDDD', '*AABBCCDE']
fragments = [Fragment(fragment) for index, fragment in enumerate(fragments)]
#fragments = [Fragment(fragment) for index, fragment in enumerate(fragments) if index % 3 in {0, 1}]

# Create graph
graph=dict()
for fragment in fragments:
    graph[fragment]=set()
    
for fragment1 in fragments:
    
    for fragment2 in fragments:
        
        if fragment1<fragment2:
            diff=fragment2.difference(fragment1)
            
            if Fragment('*') <= diff:
                diff = diff.difference(Fragment('*'))
                
            if Fragment('o') <= diff:
                diff = diff.difference(Fragment('o'))
            
            graph[fragment1].add((fragment2, 1/len(diff)))

# Find the longest path
print('backward path: {}'.format(longest_path(Fragment('*'), graph)[1:][::-1].replace(']', '?').replace('[', ']').replace('?', '[')))