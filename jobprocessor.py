#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 12:40:18 2021

@author: mpalermo
"""

import queue
from multiprocessing import Process, Queue
import os
from typing import Callable, List


class JobProcessor:
    def __init__(
        self, jobs_list: List[Callable], maxcores: int = -1, cores_per_job: int = 1
    ):

        # Set maximum total number of cores used. If user does not provide it,
        # just deduce from os settings

        self.maxcores: int

        if maxcores == -1:
            self.maxcores = len(os.sched_getaffinity(0))
        else:
            self.maxcores = maxcores

        self.cores_per_job: int = cores_per_job  # number of cores used per job

        if not all((isinstance(self.maxcores, int), self.maxcores > 0)):
            raise TypeError("maxcores must be a positive integer")
        if not all((isinstance(cores_per_job, int), cores_per_job > 0)):
            raise TypeError("cores_per_job must be a positive integer")

        # Create job queue
        self.workers: List[
            Process
        ] = []  # list of functions (jobs) to create jobs queue
        self.jobs_queue: Queue = Queue()

        for job in jobs_list:
            self.jobs_queue.put_nowait(job)

        self.number_of_jobs: int = self.jobs_queue.qsize()

        print(
            f"Running {self.number_of_jobs} jobs on {self.maxcores} cores,"
            f"{self.cores_per_job} cores per job"
        )

        self.running_jobs: Queue = Queue()  # queue for *counting* running jobs
        self._results_queue: Queue = Queue()  # queue containing jobs results

    def _job_completion_tracker(self, user_function: Callable) -> Callable:
        """
        Take user function and decorate it by setting the number of cores
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

        """

        def decorated_function():
            """
            Decorated function to be returned.
            """

            os.environ["OMP_NUM_THREADS"] = str(
                self.cores_per_job
            )  # set maximum cores per job
            result: object = user_function()  # run function and catch output
            self._results_queue.put(result)  # store the result in queue

            # The logic of the following line is broken - the counter does not
            # work but it accomplishes updating the number of running jobs,
            # so for now it's fine
            counter: int = (
                self.running_jobs.get_nowait()
            )  # update number of running jobs
            print(f"Job {counter}/{self.number_of_jobs} completed")

        return decorated_function

    def run(self) -> List:

        job_counter: int = 0

        while True:
            used_cores: int = (
                self.cores_per_job * self.running_jobs.qsize()
            )  # cores currently in use

            # Check if adding a new job would exceed the available num of cores
            if used_cores + self.cores_per_job <= self.maxcores:
                try:
                    job: Callable = self.jobs_queue.get_nowait()  # get job from queue
                    decorated_job: Callable = self._job_completion_tracker(job)
                    worker: Process = Process(target=decorated_job, args=())
                    self.workers.append(worker)
                    worker.start()
                    job_counter += 1
                    print(f"Starting job {job_counter}")
                    self.running_jobs.put_nowait(job_counter)

                except queue.Empty:
                    if self._results_queue.qsize() == self.number_of_jobs:
                        print("No more jobs to execute")
                        break

                except Exception as error:
                    print("An exception has been caught.")
                    print(error)

        for worker in self.workers:
            worker.join()

        dump: List = []

        for _ in range(self.number_of_jobs):
            dump.append(self._results_queue.get_nowait())

        return dump
