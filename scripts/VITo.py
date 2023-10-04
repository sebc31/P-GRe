# -*- coding: utf-8 -*-
import collections
import itertools

'''Various Intervals TOols contains a set of functions used to work with list of lists, i.e. "intervals"'''

def addInterval(new,index,old_intervals,new_intervals):
    if new == 0:
        new_intervals.append([old_intervals[index][0], old_intervals[index][1]])
    else:
        new_intervals.append([new, old_intervals[index][1]])

def sortIntervals(parentDic):  # sort intervals by their left position, then right position, and
    # keep only unique intervals
    if str(type(parentDic))=="<class 'collections.defaultdict'>":
        interval_list = []
        for prot in parentDic:
            for intervals in parentDic[prot]['intervals'].values():
                interval_list.append(intervals)
        interval_list = sorted(interval_list, key=lambda x: (x[0], x[1]))
    else:
        interval_list = sorted(parentDic, key=lambda x: (x[0], x[1]))
    return list(k for k, _ in itertools.groupby(interval_list))

def mergeOverlapWOFrameShift(interval_list): # Merge overlapping hit into "enlarged hits". Enlarged hits with a length not divisible
    # by 3 are supposed to contains frame-shift
    interval_fs = []  # Contains the enlarged hits
    new_start = 0  # Act both as an integer to stock the start of the new enlarged hits and has a boolean, if >0 hits
    # are being tested to be merged
    frame_shift = False
    for i in range(0, len(interval_list)):
        if i != len(interval_list) - 1:  # This condition is necessary as iterations stop before the last intervals
            if interval_list[i][1] >= interval_list[i + 1][0]:  # Compare each intervals to the next one and look
                # for overlap
                if new_start == 0:  # If a new enlarged hit is not currently being created, create a new one
                    new_start = interval_list[i][0]
                # If length of the potential enlarged hit is not divisible by 3, classified as frame_shift
                if (max(interval_list[i][1], interval_list[i + 1][1]) - new_start +1) % 3 == 0:
                    frame_shift = False
                else:
                    frame_shift = True
            # If no overlap is found between a hit (or an on-going enlarged hit) and the next one, save the hit/on-
            #going enlarged hit as a new one. If a frame-shift has occured, the next hit is not integrated into the
            # on-going enlarged hit, and the overlapping part will remain.
            else:
                addInterval(new_start, i, interval_list, interval_fs)
                new_start = 0
            if frame_shift:
                interval_fs.append([new_start, interval_list[i][1]])
                new_start = 0
                frame_shift = False
        # For last iteration
        else:
            addInterval(new_start, i, interval_list, interval_fs)
            new_start = 0
    return interval_fs

def extendedIntervalToNext(intervals_list):
    intervals_extended=[]
    for i in range(0,len(intervals_list)-1):
        intervals_extended.append([intervals_list[i][0],intervals_list[i+1][0]])
        if i==len(intervals_list)-2:
            intervals_extended.append(intervals_list[i+1])
    return intervals_extended

def protToDNAInterv(interval_list,add=0):
    for intervals in interval_list:
        intervals[0]=intervals[0]*3-2+add
        intervals[1]=intervals[1]*3+add
    return interval_list

def substract(interval,to_remove):
    new_interv=[]
    for i in range(0,len(to_remove)):
        if i==0:
            new_start = interval[0]
            new_end = to_remove[i][0]-1
        else:
            new_start = to_remove[i-1][1]
            new_end = to_remove[i][0]-1
        new_interv.append([new_start, new_end])
        if i==len(to_remove)-1:
            new_interv.append([to_remove[-1][1],interval[1]])
    return new_interv
