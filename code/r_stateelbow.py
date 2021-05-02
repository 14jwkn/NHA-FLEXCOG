#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Usage: 
    r_stateelbow.py <roiname> <width> <sublist>
    
Arguments:
    <roiname> ROI set name
    <width> Width
    <sublist> File name of subject list
    
"""

#Imports necessary packages.
import os
from docopt import docopt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from sklearn_extra.cluster import KMedoids
import time
import warnings
import scipy.sparse as sp
from collections.abc import Iterable
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics.pairwise import pairwise_distances
from yellowbrick.utils import KneeLocator, get_param_names
from yellowbrick.style.palettes import LINE_COLOR
from yellowbrick.cluster.base import ClusteringScoreVisualizer
from yellowbrick.exceptions import YellowbrickValueError, YellowbrickWarning

#Yellowbrick KMeans Elbow Visualizer source code, with metric changed from 
#Euclidean to Manhattan. 
#https://www.scikit-yb.org/en/latest/api/cluster/elbow.html#:~:text=The%20elbow%20method%20runs%20k,point%20to%20its%20assigned%20center.
try:
    from sklearn.metrics import calinski_harabasz_score as chs
except ImportError:
    from sklearn.metrics import calinski_harabaz_score as chs

__all__ = ["KElbowVisualizer", "KElbow", "distortion_score", "kelbow_visualizer"]

def distortion_score(X, labels, metric="cityblock"):
    
    # Encode labels to get unique centers and groups
    le = LabelEncoder()
    labels = le.fit_transform(labels)
    unique_labels = le.classes_

    # Sum of the distortions
    distortion = 0

    # Loop through each label (center) to compute the centroid
    for current_label in unique_labels:
        # Mask the instances that belong to the current label
        mask = labels == current_label
        instances = X[mask]

        # Compute the center of these instances
        center = instances.mean(axis=0)

        # NOTE: csc_matrix and csr_matrix mean returns a 2D array, numpy.mean
        # returns an array of 1 dimension less than the input. We expect
        # instances to be a 2D array, therefore to do pairwise computation we
        # require center to be a 2D array with a single row (the center).
        # See #370 for more detail.
        if not sp.issparse(instances):
            center = np.array([center])

        # Compute the square distances from the instances to the center
        distances = pairwise_distances(instances, center, metric=metric)
        distances = distances ** 2

        # Add the sum of square distance to the distortion
        distortion += distances.sum()

    return distortion

KELBOW_SCOREMAP = {
    "distortion": distortion_score,
    "silhouette": silhouette_score,
    "calinski_harabasz": chs,
}

class KElbowVisualizer(ClusteringScoreVisualizer):
    def __init__(
        self,
        estimator,
        ax=None,
        k=10,
        metric="distortion",
        timings=True,
        locate_elbow=True,
        **kwargs
    ):
        super(KElbowVisualizer, self).__init__(estimator, ax=ax, **kwargs)

        # Get the scoring method
        if metric not in KELBOW_SCOREMAP:
            raise YellowbrickValueError(
                "'{}' is not a defined metric "
                "use one of distortion, silhouette, or calinski_harabasz"
            )

        # Store the arguments
        self.scoring_metric = KELBOW_SCOREMAP[metric]
        self.metric = metric
        self.timings = timings
        self.locate_elbow = locate_elbow

        # Convert K into a tuple argument if an integer
        if isinstance(k, int):
            self.k_values_ = list(range(2, k + 1))
        elif (
            isinstance(k, tuple)
            and len(k) == 2
            and all(isinstance(x, (int, np.integer)) for x in k)
        ):
            self.k_values_ = list(range(*k))
        elif isinstance(k, Iterable) and all(
            isinstance(x, (int, np.integer)) for x in k
        ):
            self.k_values_ = list(k)
        else:
            raise YellowbrickValueError(
                (
                    "Specify an iterable of integers, a range, or maximal K value,"
                    " the value '{}' is not a valid argument for K.".format(k)
                )
            )

        # Holds the values of the silhoutte scores
        self.k_scores_ = None
        # Set Default Elbow Value
        self.elbow_value_ = None

    def fit(self, X, y=None, **kwargs):
        """
        Fits n KMeans models where n is the length of ``self.k_values_``,
        storing the silhouette scores in the ``self.k_scores_`` attribute.
        The "elbow" and silhouette score corresponding to it are stored in
        ``self.elbow_value`` and ``self.elbow_score`` respectively.
        This method finishes up by calling draw to create the plot.
        """

        self.k_scores_ = []
        self.k_timers_ = []
        self.kneedle = None
        self.knee_value = None

        if self.locate_elbow:
            self.elbow_value_ = None
            self.elbow_score_ = None

        for k in self.k_values_:
            # Compute the start time for each  model
            start = time.time()

            # Set the k value and fit the model
            self.estimator.set_params(n_clusters=k)
            self.estimator.fit(X, **kwargs)

            # Append the time and score to our plottable metrics
            self.k_timers_.append(time.time() - start)
            self.k_scores_.append(self.scoring_metric(X, self.estimator.labels_))

        if self.locate_elbow:
            locator_kwargs = {
                "distortion": {
                    "curve_nature": "convex",
                    "curve_direction": "decreasing",
                },
                "silhouette": {
                    "curve_nature": "concave",
                    "curve_direction": "increasing",
                },
                "calinski_harabasz": {
                    "curve_nature": "concave",
                    "curve_direction": "increasing",
                },
            }.get(self.metric, {})
            elbow_locator = KneeLocator(
                self.k_values_, self.k_scores_, **locator_kwargs
            )
            if elbow_locator.knee is None:
                self.elbow_value_ = None
                self.elbow_score_ = 0
                warning_message = (
                    "No 'knee' or 'elbow' point detected, "
                    "pass `locate_elbow=False` to remove the warning"
                )
                warnings.warn(warning_message, YellowbrickWarning)
            else:
                self.elbow_value_ = elbow_locator.knee
                self.elbow_score_ = self.k_scores_[
                    self.k_values_.index(self.elbow_value_)
                ]

        self.draw()

        return self


    def draw(self):
        """
        Draw the elbow curve for the specified scores and values of K.
        """
        # Plot the silhouette score against k
        self.ax.plot(self.k_values_, self.k_scores_, marker="D")
        if self.locate_elbow is True and self.elbow_value_ is not None:
            elbow_label = "elbow at $k={}$, $score={:0.3f}$".format(
                self.elbow_value_, self.elbow_score_
            )
            self.ax.axvline(
                self.elbow_value_, c=LINE_COLOR, linestyle="--", label=elbow_label
            )

        # If we're going to plot the timings, create a twinx axis
        if self.timings:
            self.axes = [self.ax, self.ax.twinx()]
            self.axes[1].plot(
                self.k_values_,
                self.k_timers_,
                label="fit time",
                c="g",
                marker="o",
                linestyle="--",
                alpha=0.75,
            )

        return self.ax


    def finalize(self):
        """
        Prepare the figure for rendering by setting the title as well as the
        X and Y axis labels and adding the legend.

        """
        # Get the metric name
        metric = self.scoring_metric.__name__.replace("_", " ").title()

        # Set the title
        self.set_title("{} Elbow for {} Clustering".format(metric, self.name))

        # Set the x and y labels
        self.ax.set_xlabel("k")
        self.ax.set_ylabel(metric.lower())

        # set the legend if locate_elbow=True
        if self.locate_elbow is True and self.elbow_value_ is not None:
            self.ax.legend(loc="best", fontsize="medium", frameon=True)

        # Set the second y axis labels
        if self.timings:
            self.axes[1].grid(False)
            self.axes[1].set_ylabel("fit time (seconds)", color="g")
            self.axes[1].tick_params("y", colors="g")

KElbow = KElbowVisualizer

def kelbow_visualizer(
    model,
    X,
    y=None,
    ax=None,
    k=10,
    metric="distortion",
    timings=True,
    locate_elbow=True,
    show=True,
    **kwargs
):
    klass = type(model)

    # figure out which kwargs correspond to fit method
    fit_params = get_param_names(klass.fit)

    fit_kwargs = {key: kwargs.pop(key) for key in fit_params if key in kwargs}

    oz = KElbow(
        model,
        ax=ax,
        k=k,
        metric=metric,
        timings=timings,
        locate_elbow=locate_elbow,
        **kwargs
    )
    oz.fit(X, y, **fit_kwargs)

    if show:
        oz.show()
    else:
        oz.finalize()

    return oz

#Gets the window vectors for the subject.
def makedFCvec(subject):
    
    #Gets the vectorized dFC matrices. 
    subvec = pd.read_csv('../outputs/r_vectorized/'+roiname+'/'+width+
                             '/'+subject+'/vector.csv',header=None)
    
    #Converts to numpy array for future functions.
    subvec = subvec.to_numpy()
    
    #Return window vectors for the subject.
    return subvec

#Finds the exemplar windows which have the local maxima of variance over all
#connections.
def findexemplars(labsubwins):
    
    #Gets subject ID from first cell. 
    subject = labsubwins[0,0]
    
    #Separates IDs from windows for processing. Remove unused variables to 
    #save memory.
    subwins = labsubwins[:,1:]
    del labsubwins
    
    #For each window, append the variance of the connectivity values for 
    #that window to a list.
    subwinvar = []
    for winval in subwins:
        subwinvar.append(winval.var())
    
    #Find the window indices of the local maxima of variance from the list 
    #of variances.
    varwinindices = find_peaks(np.array(subwinvar))
    
    #Append the windows with local maxima of variance and their subject-labelled
    #indices to a list.
    exemwins = subwins[varwinindices[0]] #Windows.
    exemindices = np.insert(varwinindices[0],0,subject) #Label before indices.
    
    #Returns labelled exemplar windows.
    return [exemwins,exemindices]

#Produces the graph to find the elbow, using the Manhattan distance 
#metric and the build algorithm to select initial medoids which is constant
#for a set of data. Elbow is labelled on the graph using an algorithm.
def elbow():
    model = KMedoids(metric='cityblock') #Gets model.
    visualizer = KElbowVisualizer(model,k=(2,20)) #Produce visualizer to test 2-20 k.
    visualizer.fit(allwins) #Fits model to data for visualizing.
    visualizer.show(outpath=basepath+'elbow.jpg') #Saves.

#Required for parallel processing function.
if __name__ == '__main__':
    __spec__ = None
    
    #Catches command-line arguments.
    args = docopt(__doc__)
    
    #Imports ROI name.
    roiname = args['<roiname>']
    
    #Imports window width.
    width = args['<width>']
    
    #Reads in chosen list of subject IDs.
    sublist = args['<sublist>']
    with open(sublist) as f:
        subjects = [subject.rstrip() for subject in f]
    
    #Combine subject window vectors. Remove unused variables to save memory.
    #Convert list to a numpy array to save memory.
    subdFCvec = []
    for onesub in subjects: #For each subject.
        onesubvec = makedFCvec(onesub) #Produce subject window vectors.
        subdFCvec.append(onesubvec)
    del onesubvec
    subdFCvec = np.array(subdFCvec)
    
    #Sets output folder path. If output folder doesn't exist, creates it.
    basepath = '../outputs/r_stateelbow/'+roiname+'/'+width+'/'
    if not os.path.exists(basepath):   
        os.makedirs(basepath)
    
    #Find exemplars. Removes unused variables to save memory.
    exemlist = []
    for subwindows in subdFCvec: #For each subject's window vectors.
        subexem = findexemplars(subwindows) #Find exemplar windows of high variance.
        exemlist.append(subexem)
    del subwindows
    del subexem
    del subdFCvec
    
    #Unpacks windows and indices. Removes unused variables to save memory.
    exemvec = []
    exemindices = []
    for exempack in exemlist: #For each exemplar pack.
        exemvec.append(exempack[0]) #Appends exemplar windows.
        exemindices.append(exempack[1]) #Appends exemplar indices.
    del exempack
    del exemlist
    
    #Saves the window indices with local maxima of variance as a csv file. 
    #Removes unused variables to save memory.
    pd.DataFrame(exemindices).to_csv(basepath+'varwinindex.csv',header=None,
                                       index=0)
    del exemindices
    
    #Takes the exemplar windows from each subject and removes subject level.
    #Removes unused variables to save memory. Convert list to a numpy array 
    #to save memory. 
    allwins = []
    for subval in exemvec: 
        for winval in subval:
            allwins.append(winval)
    del subval
    del winval
    del exemvec
    allwins = np.array(allwins)
    
    #Produces elbow graph.
    elbow()
    