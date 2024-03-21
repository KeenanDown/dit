"""Tools for computing weakly-shared information between variables.
That is, information that is produced by common generators.
For more information see the preprint:
https://arxiv.org/abs/2305.07554
"""
# Specify all functions defined in this module.
__all__=[
    'weakly_shared_content',
    'weakly_shared',
    'strongly_shared_content',
    'strongly_shared'
]

# Do some imports.
import warnings
from ..npdist import Distribution
from .content import content, measure
from .combinatorics import get_n_atoms, upper_set

def weakly_shared_content(dist, list_of_rv_groups, order = 2):
    """
    Returns the weakly-shared information as a set derived from the logarithmic decomposition.

    Parameters
    ----------
    dist : dit.Distribution
        The distribution to be examined.
    list_of_rv_groups : list
        A list of lists, containing random variable groupings for computing shared information.
    order : int, string
        Default value is 2. If set to n, redundancy is computed using n-atom upper sets.
        If given as "even", then will be generated for all even n-atoms. Likewise for "odd".
    log_base : int, float
        The base of the logarithm used for computing information.

    Returns
    -------
    shared_information : set
        The information atoms weakly shared between the groups of variables.

    Examples
    --------
    To see what information atoms X, Y and Z weakly share, do [["X"], ["Y"], ["Z"]]
    For information atoms they jointly share with Z, do [["X", "Y"], ["Z"]].
    """
    # Check the inputs are the correct types.
    if not isinstance(dist, Distribution):
        raise TypeError("'dist' must be a dit distribution.")
    if not isinstance(list_of_rv_groups, list):
        raise TypeError("'list_of_rv_groups' must be a list of lists of random variables.")
    if not isinstance(order, (int, str)):
        raise TypeError("'order' must be either an integer, 'even' or 'odd'.")
    if isinstance(order, str):
        if order not in ["even", "odd"]:
            raise ValueError("'order' must be either an integer, 'even' or 'odd'.")
    # Compute the relevant content.
    relevant_contents = {frozenset(content(dist, group)) for group in list_of_rv_groups}
    # Initialise a large set.
    intersection_set = relevant_contents.pop()
    # Loop over an intersection.
    for content_set in relevant_contents:
        intersection_set = intersection_set.intersection(content_set)
    # Take the intersection_set and find the relevant atoms.
    filtered_atoms = get_n_atoms(set(intersection_set), order)
    # Compute the upper set of these atoms.
    set_to_measure = upper_set(dist, filtered_atoms)
    # Return the set_to_measure
    return set_to_measure

def weakly_shared(dist, list_of_rv_groups, order = 2, log_base = 2):
    """
    Returns the weakly-shared information derived from the logarithmic decomposition
    shared between the groups of random variables in 'list_of_rv_groups'.

    Parameters
    ----------
    dist : dit.Distribution
        The distribution to be examined.
    list_of_rv_groups : list
        A list of lists, containing random variable groupings for computing shared information.
    order : int, string
        Default value is 2. If set to n, redundancy is computed using n-atom upper sets.
        If given as "even", then will be generated for all even n-atoms. Likewise for "odd".
    log_base : int, float
        The base of the logarithm used for computing information.

    Returns
    -------
    shared_information : float
        The information shared between the groups of variables.

    Examples
    --------
    To see what information X, Y and Z weakly share, do [["X"], ["Y"], ["Z"]]
    For information they give together about Z, do [["X", "Y"], ["Z"]].
    """
    # Check the inputs are the correct types.
    if not isinstance(dist, Distribution):
        raise TypeError("'dist' must be a dit distribution.")
    if not isinstance(list_of_rv_groups, list):
        raise TypeError("'list_of_rv_groups' must be a list of lists of random variables.")
    if not isinstance(order, (int, str)):
        raise TypeError("'order' must be either an integer, 'even' or 'odd'.")
    if not isinstance(log_base, (int, float)):
        raise TypeError("'log_base' must be an int or a float.")
    if isinstance(order, str):
        if order not in ["even", "odd"]:
            raise ValueError("'order' must be either an integer, 'even' or 'odd'.")
    # Get the set.
    set_to_measure = weakly_shared_content(dist, list_of_rv_groups, order)
    # Measure the set.
    shared_information = measure(dist, set_to_measure, log_base)
    # Return this.
    return shared_information

def strongly_shared_content(dist, list_of_rv_groups, order = 2):
    """
    Returns the strongly-shared information derived from the logarithmic decomposition
    shared between the groups of random variables in 'list_of_rv_groups'.

    All of these atoms exhibit a pure "co-change" between input variables.
    Parameters
    ----------
    dist : dit.Distribution
        The distribution to be examined.
    list_of_rv_groups : list
        A list of lists, containing random variable groupings for computing shared information.
    order : int, string
        Default value is 2. If set to n, redundancy is computed using n-atom upper sets.
        If given as "even", then will be generated for all even n-atoms. Likewise for "odd".
    log_base : int, float
        The base of the logarithm used for computing information.

    Returns
    -------
    shared_information : float
        The information shared between the groups of variables.

    Examples
    --------
    To see what information X, Y and Z strongly share, do [["X"], ["Y"], ["Z"]]
    For information XY strongly shares with Z, do [["X", "Y"], ["Z"]].
    """

