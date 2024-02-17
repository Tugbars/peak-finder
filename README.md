# Peak Finding Algorithm Snippet


## Description:
This repository contains a stripped-down version of a peak finding algorithm initially developed as part of a larger embedded systems project. The code demonstrates my ability to write efficient algorithms and manage data processing tasks within constrained environments. Although the code is isolated from its original context and lacks the complete functionality and integrations of the full system, it serves as a showcase of algorithm development skills, particularly in identifying peaks within a dataset, calculating their prominences, and organizing these values in a vector for further analysis.

## Key Features:
Peak Detection: Efficiently scans a dataset to identify peaks based on predefined criteria.
Prominence Calculation: Determines the prominence of each peak, offering insights into the relative height or significance of peaks within the data.

## Note to Viewers:
This code snippet is extracted from a larger, complex embedded system project. As such, it operates independently of the project's other components, with most of its functions, helper functions, and linkage methods to the main application removed. This standalone version is intended to illustrate my proficiency in algorithm writing and problem-solving within embedded systems, rather than serve as a fully functional application.

## Usage:
The algorithm is designed for educational and demonstration purposes and can be adapted or expanded to fit specific data analysis needs. Users are encouraged to explore the code and consider potential applications or improvements within their projects.

## Peak Finding and Prominence Calculation:

The core algorithm first differentiates the input dataset to find points where the slope changes (potential peaks). It then evaluates these points to confirm they are peaks based on the change in direction (from increasing to decreasing).
For each identified peak, the algorithm calculates its prominence, which is the height of the peak relative to the highest contour line encircling it but not containing any higher peak.
The algorithm includes an option to consider or exclude the dataset's endpoints as potential peaks, controlled by the includeEndpoints parameter.
The extrema parameter allows the algorithm to focus on maxima (extrema > 0) or minima (extrema < 0), adjusting the dataset's sign accordingly.
Selection of the Most Significant Peak:

After identifying peaks and calculating their prominences, the comparator function is used to select the peak with the highest prominence, considering the extrema direction.
