#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 12:40:18 2021

@author: mpalermo
"""

import queue
from multiprocessing import Process, Queue, Value
import os

class JobProcessor:
    def __init__(self, jobs_list, maxcores=None, cores_per_job=1):
        
        # Set maximum total number of cores used. If user does not provide it,
        # just deduce from os settings
        self.maxcores= maxcores
        
        if maxcores == None:
            self.maxcores = len(os.sched_getaffinity(0))

        self.cores_per_job = cores_per_job # number of cores used per job 

        # Create job queue
        self.jobs = [] # list of jobs to create jobs queue and track workers
        self.jobs_queue = Queue()        
        [ self.jobs_queue.put_nowait(job) for job in jobs_list ]
        
        
        self.number_of_jobs = self.jobs_queue.qsize()
        
        print(f"Running {self.number_of_jobs} jobs on {self.maxcores} cores, {self.cores_per_job} cores per job")
        
        self.running_jobs = Queue() # queue for *counting* running jobs
        self.running_jobs_value = Value('i', 0)
        self._results_queue = Queue() # queue containing jobs results
        

    def _job_completion_tracker(self, user_function):
        '''
        Takes user function and decorate it by setting the number of cores
        available to the jobs, gather the results in a common queue and update
        the number of available cores once the functions terminates.

        Parameters
        ----------
        user_function : method
            User function to be run 

        Returns
        -------
        decorated_function : method
            Function decorated for the inner working of the program.
            .

        '''
        def decorated_function():
            '''
            Decorated function to be returned.
            '''
            
            os.environ["OMP_NUM_THREADS"] = str(self.cores_per_job) #set maximum cores per job
            result = user_function() # run user function and catch the output
            self._results_queue.put(result) # store the result in queue
            counter = self.running_jobs.get_nowait() # update number of running jobs
            print(f"Job {counter}/{self.number_of_jobs} completed")
            
        return decorated_function
        
    def run(self):
        
        job_counter = 0
        used_cores = self.cores_per_job*self.running_jobs.qsize() # cores currently in use
        
        while True:
            # Check if adding a new job would exceed the available num of cores
            if used_cores+self.cores_per_job <= self.maxcores:
                try:
                    job = self.jobs_queue.get_nowait() # get job from queue
                    decorated_job = self._job_completion_tracker(job)
                    worker = Process(target=decorated_job, args=())
                    self.jobs.append(worker)
                    worker.start()
                    job_counter += 1
                    print(f"Starting job {job_counter}")
                    self.running_jobs.put_nowait(job_counter)
                
                except queue.Empty:
                    if self._results_queue.qsize() == self.number_of_jobs:
                        print("No more jobs to execute")
                        break
                    else:
                        pass
                except Exception as error:
                    print("An exception has been caught.")
                    print(error)
           
        for worker in self.jobs:
            worker.join()
            
        dump = []
        
        for _ in range(self.number_of_jobs):
            dump.append(self._results_queue.get_nowait())
        
        return dump
