#include "mes_peakfinder.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <iterator>
#include "mqs_def.h"

struct peak OldfindPeaks(std::vector<float> x0, bool includeEndpoints, float extrema);

///////////////////////
// STATIC VARIABLES ///
///////////////////////
static std::vector<float> dataVec;
static std::vector<peak> peakHolder;

/////////////////////
// related functions //
///////////////////////

peak comparator(std::vector<peak> minmag, int extrema)
{
	float max = 0.0f;
	float min = 0.0f;
	int index = 0;
	for (int i = 0; i < minmag.size(); i++)
	{
		//printf("mag %f prom %f loc %d\n", minmag[i].magnitude, minmag[i].prominence, minmag[i].location);
		if (minmag[i].magnitude > minmag[i].prominence + 3.0f)
		{
			if ((float)minmag[i].prominence > max && extrema == 1.0f)
			{
				max = (float)minmag[i].prominence;
				index = i;
			}
			if ((float)minmag[i].prominence < min && extrema == -1.0f)
			{
				min = (float)minmag[i].prominence;
				index = i;
			}
		}
	}
	printf("comparator: index is %d, mag is %f, prominence is %f\n", minmag[index].location, minmag[index].magnitude, minmag[index].prominence);
	return minmag[index];
}

/**
 * @brief Calculates the difference between consecutive elements in a vector.
 * 
 * @param in Input vector of floats.
 * @param out Output vector containing the difference between consecutive elements of `in`.
 */
void diff(std::vector<float> in, std::vector<float>& out)
{
	out = std::vector<float>(in.size() - 1);
	for (int i = 1; i < in.size(); ++i)
		out[i - 1] = in[i] - in[i - 1];
}

/**
 * @brief Calculates the element-wise product of two vectors.
 * 
 * @param a First input vector of floats.
 * @param b Second input vector of floats.
 * @param out Output vector containing the element-wise product of `a` and `b`.
 */
void vectorElementsProduct(std::vector<float> a, std::vector<float> b, std::vector<float>& out)
{
	out = std::vector<float>(a.size());

	for (int i = 0; i < a.size(); ++i)
		out[i] = a[i] * b[i];
}

/**
 * @brief Finds indices of vector elements greater than a specified threshold.
 * 
 * @param in Input vector of floats.
 * @param threshold Float threshold for comparison.
 * @param indices Output vector of indices where `in[i]` is greater than `threshold`.
 */
void findIndicesMoreThan(std::vector<float> in, float threshold, std::vector<int>& indices)
{
	for (int i = 0; i < in.size(); ++i)
		if (in[i] > threshold)
			indices.push_back(i + 1);
}

/**
 * @brief Finds indices of vector elements less than a specified threshold.
 * 
 * @param in Input vector of floats.
 * @param threshold Float threshold for comparison.
 * @param indices Output vector of indices where `in[i]` is less than `threshold`.
 */
void findIndicesLessThan(std::vector<float> in, float threshold, std::vector<int>& indices)
{
	for (int i = 0; i < in.size(); ++i)
		if (in[i] < threshold)
			indices.push_back(i + 1);
}

/**
 * @brief Selects elements from an input vector based on specified indices.
 * 
 * @param in Input vector from which to select elements.
 * @param indices Vector of indices of elements to select from `in`.
 * @param out Output vector containing selected elements.
 */
void selectElementsFromIndices(std::vector<float> in, std::vector<int> indices, std::vector<float>& out)
{
	for (int i = 0; i < indices.size(); ++i)
		out.push_back(in[indices[i]]);
}

/**
 * @brief Selects elements from an input vector of integers based on specified indices.
 * 
 * @param in Input vector of integers from which to select elements.
 * @param indices Vector of indices of elements to select from `in`.
 * @param out Output vector of integers containing selected elements.
 */
void selectElementsFromIndices(std::vector<int> in, std::vector<int> indices, std::vector<int>& out)
{
	for (int i = 0; i < indices.size(); ++i)
		out.push_back(in[indices[i]]);
}

/**
 * @brief Assigns a sign to each element of the input vector: 1 for positive, -1 for negative, and 0 for zero values.
 * 
 * @param in Input vector of floats.
 * @param out Output vector of integers with signs corresponding to elements of `in`.
 */
void signVector(std::vector<float> in, std::vector<int>& out)
{
	out = std::vector<int>(in.size());

	for (int i = 0; i < in.size(); ++i)
	{
		if (in[i] > 0)
			out[i] = 1;
		else if (in[i] < 0)
			out[i] = -1;
		else
			out[i] = 0;
	}
}

/**
 * @brief Multiplies each element of the input vector by a scalar value.
 * 
 * @param scalar Scalar value to multiply with each element of `in`.
 * @param in Input vector of floats.
 * @param out Output vector of floats, each element of `in` multiplied by `scalar`.
 */
void scalarProduct(float scalar, std::vector<float> in, std::vector<float>& out)
{
	out = std::vector<float>(in.size());

	for (int i = 0; i < in.size(); ++i)
		out[i] = scalar * in[i];
}

/**
 * @brief Finds the previous bin with a magnitude greater than or equal to the current bin's magnitude.
 * 
 * @param bins Input vector of bin magnitudes.
 * @param index Index of the current bin.
 * @return Index of the previous bin with magnitude >= current bin's magnitude, or 0 if not found.
 */
int find_prev_bin_with_magnitude(std::vector<float> bins, int index)
{
	float magnitude = bins[index];
	for (int bin_index = index - 1; bin_index > 1; bin_index--)
	{
		// looking for a bigger peak
		if (bins[bin_index] >= magnitude && bins[bin_index] < 900000.0f)
		{
			return bin_index;
		}
	}
	return 0; // couldn't find
}

/**
 * @brief Finds the next bin with a magnitude greater than or equal to the current bin's magnitude.
 * 
 * @param bins Input vector of bin magnitudes.
 * @param index Index of the current bin.
 * @param num_bins Total number of bins.
 * @return Index of the next bin with magnitude >= current bin's magnitude, or last bin index if not found.
 */
int find_next_bin_with_magnitude(std::vector<float> bins, int index, int num_bins)
{
	float magnitude = bins[index];
	for (int bin_index = index + 1; bin_index < num_bins; bin_index++)
	{
		if (bins[bin_index] >= magnitude && bins[bin_index] < 900000.0f)
		{
			return bin_index;
		}
	}
	return num_bins - 1; // couldn't find
}

/**
 * @brief Finds the minimum magnitude within a specified range of bins.
 * 
 * @param bins Input vector of bin magnitudes.
 * @param index_left Starting index of the range.
 * @param index_right Ending index of the range.
 * @return Minimum magnitude within the specified range.
 */
float get_min_magnitude_in_range(std::vector<float> bins, int index_left, int index_right)
{
	float min_magnitude = 1e6;
	for (int bin_index = index_left; bin_index < index_right; bin_index++)
	{
		if (bins[bin_index] < min_magnitude)
		{
			min_magnitude = bins[bin_index];
		}
	}
	return min_magnitude;
}

/**
 * @brief Calculates the prominence of a peak at a given location.
 * 
 * @param x0 Input vector of data.
 * @param index_holder Vector of peak indices.
 * @param tempLoc Index of the current peak whose prominence is to be calculated.
 * @return The prominence of the peak at `tempLoc`.
 */
float prominenceCalc(std::vector<float> x0, std::vector<int> index_holder, int tempLoc)
{
	// read MATLAB's prominence algorithm explanation to learn more about the calculations below
	int index_left = find_prev_bin_with_magnitude(x0, index_holder[tempLoc]);
	int index_right = find_next_bin_with_magnitude(x0, index_holder[tempLoc], x0.size());
	float min_contour_left = get_min_magnitude_in_range(x0, index_left, index_holder[tempLoc]);
	float min_contour_right = get_min_magnitude_in_range(x0, index_holder[tempLoc], index_right);
	// reference level for the prominence
	float max_min_contour = std::max(min_contour_left, min_contour_right);
	if (max_min_contour == 0)
		max_min_contour = 1e-9;
	float prominence = x0[index_holder[tempLoc]] - max_min_contour;
	return prominence;
}


/**
 * @brief Change of Slope Detection Method for Peak Identification
 *
 * This function implements a heuristic method for detecting peaks in a dataset by calculating the product of 
 * consecutive differences (first derivatives) and identifying locations where this product is negative. This 
 * approach is akin to finding points of inflection where the slope of the function changes sign, without directly 
 * computing the second derivative but leveraging a concept related to it.
 * 
 * In the context of digital signal processing (DSP) and numerical analysis, this method is seen as a simplified 
 * or heuristic approach for peak detection. It is particularly effective when the goal is to detect changes in the 
 * trend of the data rather than identifying exact inflection points through the zero points of the second derivative, 
 * which is a common approach for continuous functions.
 * 
 * The method identifies potential peaks by exploiting the fact that at a peak (either a maximum or a minimum), the 
 * direction of the slope changes. By calculating the product of consecutive slopes (differences), a negative result 
 * indicates a directional change, characteristic of a peak. This computational approach is simpler than a true second 
 * derivative test, making it advantageous in discrete datasets where direct derivative calculations are more 
 * challenging and susceptible to noise.
 * 
 * This technique, described as a change of slope detection method, balances computational efficiency with the need 
 * for identifying significant changes in the data's trend. It is particularly suitable for embedded systems, where 
 * limited resources necessitate efficient algorithms.
 *
 * @param dataVec A std::vector<float> containing the dataset to analyze.
 * @return std::vector<int> containing the indices of detected peaks, based on the change of slope detection method.
 */
struct peak OldfindPeaks(std::vector<float> x0, bool includeEndpoints, float extrema)
{
	float selectivity = 0.0f;
	int vector_length = x0.size();

	scalarProduct(extrema, x0, x0); // maxima or minima? extrema determines that. 

	std::vector<float> derivative;
	// compare scalar products of each item to the one before and put it inside dx vector
	diff(x0, derivative);
	// clears machine epsilons.
	std::replace(derivative.begin(), derivative.end(), 0.0f, EPS); // not usable for now.
	std::vector<float> derivative0(derivative.begin(), derivative.end() - 1); // from 0 to n-1 vector
	std::vector<float> derivative1(derivative.begin() + 1, derivative.end()); // from 1 to n vector
	std::vector<float> second_derivative;

	vectorElementsProduct(derivative0, derivative1, second_derivative); //multiply instead of following 
	std::vector<int> index_holder;
	// Find where the derivative changes sign
	findIndicesLessThan(second_derivative, 0, index_holder); // 0 is for threshold.

	//*************************************************************************************************************************************************
	// FIND THE MAXIMA
	std::vector<float> x;
	float leftMin;
	int minMagIdx;
	float minMag;
	// this gives the indexes of the values which are negative. it makes a list out of the indexes of the negative leaning derivatives

	if (includeEndpoints)
	{
		selectElementsFromIndices(x0, index_holder, x);			// select the peaks from the x0 and output as x.
		x.insert(x.begin(), x0[0]);								// add start point
		x.insert(x.end(), x0[x0.size() - 1]);					// add end point to the vector which holds all the peaks
		index_holder.insert(index_holder.begin(), 1);			// prepare index holder start point
		index_holder.insert(index_holder.end(), vector_length); // prepare index holder end point
		minMagIdx = std::distance(x.begin(), std::min_element(x.begin(), x.end()));
		minMag = x[minMagIdx];
		leftMin = minMag;
	}
	else
	{
		selectElementsFromIndices(x0, index_holder, x); // based on the indices it keeps pushing the values.
		if (x.size() > 2)
		{
			minMagIdx = std::distance(x.begin(), std::min_element(x.begin(), x.end())); // the distance of the minima inside the selected 2nd derivative array. the index of it.
			minMag = x[minMagIdx];														// magnitude of the minima
			leftMin = x[0] < x0[0] ? x[0] : x0[0];
		}
	}

	int len = x.size();
	// holds the location of the peak for next comparison temporarily.
	int tempLoc = 0;
	// an int to traverse through minimas and keep their indices in record.
	int index_indicator;

	// Skip the first point if it is smaller so we always start on maxima
	if (x[0] >= x[1])
		index_indicator = 0;
	else
		index_indicator = 1;
	float tempMag = minMag; // correct initialization.

	if (len == 2)
	{
		if (x[index_indicator] > x[index_indicator - 1])
		{
			tempLoc = index_indicator;
			tempMag = x[index_indicator];
			float max_prominence3 = prominenceCalc(x0, index_holder, tempLoc);		 // calculate the prominence of the minima/maxima
			struct peak maxtemp = { tempMag, index_holder[tempLoc], max_prominence3 }; // gather all the properties of the minima in a struct.
			peakHolder.push_back(maxtemp);											 // pysh
		}
	}
	if (len > 2)
	{
		bool foundPeak = false;
		if (includeEndpoints)
		{
			// Deal with first point a little differently since tacked it on
			// Calculate the sign of the derivative since we tacked the first
			//  point on it does not neccessarily alternate like the rest.
			std::vector<float> xSub0(x.begin(), x.begin() + 3); // subvector
			std::vector<float> xDiff;							// subvector
			diff(xSub0, xDiff);

			std::vector<int> signDx;

			signVector(xDiff, signDx); // finds the differences in mutlak form.

			if (signDx[0] <= 0) // The first point is larger or equal to the second
			{
				if (signDx[0] == signDx[1]) // Want alternating signs
				{
					x.erase(x.begin() + 1);
					index_holder.erase(index_holder.begin() + 1);
					len = len - 1;
				}
			}
			else // First point is smaller than the second
			{
				if (signDx[0] == signDx[1]) // Want alternating signs
				{
					x.erase(x.begin());
					index_holder.erase(index_holder.begin());
					len = len - 1;
				}
			}
		}
		// Skip the first point if it is smaller so we always start on maxima
		if (x[0] >= x[1])
			index_indicator = 0;
		else
			index_indicator = 1;

		// check the corresponding number of the len.
		while (index_indicator < len)
		{
			index_indicator = index_indicator + 1; // This is a peak
			// Reset peak finding if we had a peak and the next peak is bigger
			// than the last or the left min was small enough to reset.
			if (foundPeak)
			{
				tempMag = minMag;
				foundPeak = false;
			}

			// compare the new peak to the previous peak, and compare the new peak to the dip of the previous peak + selectivity.
			if (x[index_indicator - 1] > tempMag && x[index_indicator - 1] > leftMin + selectivity)
			{
				tempLoc = index_indicator - 1;
				tempMag = x[index_indicator - 1];
			}

			// Make sure we don't iterate past the length of our vector
			if (index_indicator == len)
				break; // We assign the last point differently out of the loop
			// move on to the valley in order to check if the tempMag(the peak which we are assessing) is higher than the peak in that valley.
			index_indicator = index_indicator + 1; // Move onto the valley

			// Selectivity has not been utilized
			if (!foundPeak && tempMag > selectivity + x[index_indicator - 1])
			{
				foundPeak = true; // We have found a peak
				leftMin = x[index_indicator - 1];
				float prominence2 = prominenceCalc(x0, index_holder, tempLoc);
				// put the values inside the struct.
				struct peak temp = { tempMag, index_holder[tempLoc] + 1, prominence2 };
				peakHolder.push_back(temp);
			}
			else if (x[index_indicator - 1] < leftMin) // New left minima
				// this is the peak one before the previous maxima.
				leftMin = x[index_indicator - 1];
		}

		// Check end point
		if (includeEndpoints)
		{
			// if peak was found but this is bigger
			if (x[x.size() - 1] > tempMag && x[x.size() - 1] > leftMin + selectivity)
			{
				tempLoc = len - 1;
				tempMag = x[x.size() - 1];
				float prominence_end_point = prominenceCalc(x0, index_holder, tempLoc);
				struct peak end_temp = { tempMag, index_holder[tempLoc], prominence_end_point };
				peakHolder.push_back(end_temp);
			}
			// if no peak ever found.
			else if (!foundPeak && tempMag > minMag) // Check if we still need to add the last point
			{
				/*	float prominence_end_point = prominenceCalc(x0, index_holder, tempLoc);
					struct peak end_temp = { tempMag, index_holder[tempLoc], prominence_end_point };
					peakHolder.push_back(end_temp);*/
			}
		}
		else if (!foundPeak)
		{
			////both conditions above are not met, yet we are checking if it is worth to add endpoints as peaks.
			float minAux = x0[x0.size() - 1] < x[x.size() - 1] ? x0[x0.size() - 1] : x[x.size() - 1];
			if (x[x.size() - 1] > tempMag && x[x.size() - 1] > leftMin + selectivity)
			{
				tempLoc = len - 1;
				tempMag = x[x.size() - 1];
				float prominence_end_point = prominenceCalc(x0, index_holder, tempLoc);
				struct peak end_temp = { tempMag, index_holder[tempLoc], prominence_end_point };
				peakHolder.push_back(end_temp);
			}
			else if (!tempMag > minAux + selectivity) // Check if we still need to add the last point
			{
				/*	float prominence_end_point = prominenceCalc(x0, index_holder, tempLoc);
					struct peak end_temp = { tempMag, index_holder[tempLoc], prominence_end_point };
					peakHolder.push_back(end_temp);*/
			}
		}
	}
	return comparator(peakHolder, extrema); // do I have to store peakHolder statically? probably not. This is the next thing to change.
}
