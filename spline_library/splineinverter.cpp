#include "splineinverter.h"

#include "spline.h"
#include "spline_library/utils/splinesample_adaptor.h"
#include "spline_library/utils/nanoflann.hpp"

#include <algorithm>

#define PI 3.14159265359

template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}
template<typename T> T abs(T val) {
    return val * sign(val);
}

//inner class used to provide an abstraction between the spline inverter and nanoflann
class SplineInverter::SampleTree
{
    typedef SplineSampleAdaptor<SplineSamples3D> AdaptorType;
    typedef nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<double, AdaptorType>,
            AdaptorType, 3>
        TreeType;

public:
    SampleTree(const SplineSamples3D &samples)
        :adaptor(samples), tree(3, adaptor)
    {
        tree.buildIndex();
    }

    double findClosestSample(const Vector3D &queryPoint) const
    {
        //build a query point
        double queryPointArray[3] = {queryPoint.x(), queryPoint.y(), queryPoint.z()};

        // do a knn search
        const size_t num_results = 1;
        size_t ret_index;
        double out_dist_sqr;
        nanoflann::KNNResultSet<double> resultSet(num_results);
        resultSet.init(&ret_index, &out_dist_sqr );
        tree.findNeighbors(resultSet, &queryPointArray[0], nanoflann::SearchParams());

        return adaptor.derived().pts.at(ret_index).t;
    }

private:
    AdaptorType adaptor;
    TreeType tree;
};

SplineInverter::SplineInverter(const std::shared_ptr<Spline> &spline, int samplesPerT)
    :spline(spline), sampleStep(1.0 / double(samplesPerT)), slopeTolerance(0.01)
{
    SplineSamples3D samples;

	//first step is to populate the splineSamples map
	//we're going to have sampled T values sorted by x coordinate
	double currentT = 0;
	double maxT = spline->getMaxT();
	while(currentT < maxT)
	{
        auto sampledPoint = spline->getPosition(currentT);

        SplineSamples3D::Point3D currentSample(sampledPoint.x(), sampledPoint.y(), sampledPoint.z(), currentT);
        samples.pts.push_back(currentSample);

		currentT += sampleStep;
	}

	//if the spline isn't a loop and the final t value isn't very very close to maxT, we have to add a sample for maxT
    double lastT = samples.pts.at(samples.pts.size() - 1).t;
    if(!spline->isLooping() && abs(lastT / maxT - 1) > .0001)
	{
		auto sampledPoint = spline->getPosition(maxT);

        SplineSamples3D::Point3D currentSample(sampledPoint.x(), sampledPoint.y(), sampledPoint.z(), maxT);
        samples.pts.push_back(currentSample);
    }

    //populate the sample kd-tree
    sampleTree = std::unique_ptr<SampleTree>(new SampleTree(samples));
}

SplineInverter::~SplineInverter()
{

}

double SplineInverter::findClosestT(const Vector3D &queryPoint) const
{
    double closestSampleT = sampleTree->findClosestSample(queryPoint);

    double sampleDistanceSlope = getDistanceSlope(queryPoint, closestSampleT);

    //if the slope is very close to 0, just return the sampled point
    if(abs(sampleDistanceSlope) < slopeTolerance)
        return closestSampleT;

    //if the spline is not a loop there are a few special cases to account for
    if(!spline->isLooping())
    {
        //if closest sample T is 0, we are on an end. so if the slope is positive, we have to just return the end
        if(abs(closestSampleT) < .0001 && sampleDistanceSlope > 0)
            return closestSampleT;

        //if the closest sample T is max T we are on an end. so if the slope is negative, just return the end
        if(abs(closestSampleT / spline->getMaxT() - 1) < .0001 && sampleDistanceSlope < 0)
            return closestSampleT;
    }


    //step forwards or backwards in the spline until we find a point where the distance slope has flipped sign
    //because "currentsample" is the closest point,  the "next" sample's slope MUST have a different sign
    //otherwise that sample would be closer
    //note: this assumption is only true if the samples are close together

    double a = closestSampleT;

    //if sample distance slope is positive we want to move backwards in t, otherwise forwards
    double b = closestSampleT - sampleStep * sign(sampleDistanceSlope);

    double aValue = sampleDistanceSlope;
    double bValue = getDistanceSlope(queryPoint, b);

    //we know that the actual closest point is now between littleT and bigT
    //use the circle projection method to find the actual closest point, using the Ts as bounds
    return brentsMethod(queryPoint, a, aValue, b, bValue);
}

double SplineInverter::brentsMethod(const Vector3D &queryPoint, double contrapoint, double contrapointValue, double currentGuess, double currentGuessValue) const
{
    // http://en.wikipedia.org/wiki/Brent%27s_method
    // this algorithm is sort of a hybrid between the secant method and bisection method
    // the bisection method is too slow but is guaranteed to converge, the secant method is fast but not guaranteed to converge
    // this combines the good of both and leaves out the bad of both, but it's much more complicated to read

    //if a is a better guess than b, swap them - we want b to be better than a
    if(abs(contrapointValue) < abs(currentGuessValue))
    {
        std::swap(contrapoint, currentGuess);
        std::swap(contrapointValue, currentGuessValue);
    }

    bool mflag = true;
    double prevGuess = contrapoint;

    double prevGuessValue = contrapointValue;

    //this gets set at the end of the loop, but is only referenced if mflag is false
    //and mflag starts out true, so i can't be referenced while uninitialized
    double oldGuess;

    double minDelta = 0.001;

    while(currentGuessValue > slopeTolerance && abs(currentGuess - contrapoint) > slopeTolerance)
    {
        //s will be the next guess for the actual t value
        double nextGuess;

        if(contrapointValue != prevGuessValue && currentGuessValue != prevGuessValue)
        {
            //if f(a),f(b) and f(c) are all distinct, use inverse quadractic interpolation
            nextGuess = contrapoint * currentGuessValue * prevGuessValue / ((contrapointValue - currentGuessValue) * (contrapointValue - prevGuessValue));
            nextGuess += currentGuess * contrapointValue * prevGuessValue / ((currentGuessValue - contrapointValue) * (currentGuessValue - prevGuessValue));
            nextGuess += prevGuess * contrapointValue * currentGuessValue / ((prevGuessValue - contrapointValue) * (prevGuessValue - currentGuessValue));
        }
        else
        {
            //use secant method
            nextGuess = currentGuess - currentGuessValue * (currentGuess - contrapoint) / (currentGuessValue - contrapointValue);
        }


        //determine if we can use the s that we jsut calculated. if not we have to use bisection method :(
        //condition numbers correspond to the pseudocode in the wiki article
        if(     (nextGuess < (3 * contrapoint + currentGuess) / 4 || nextGuess  > currentGuess) //condition 1
                || (mflag && (abs(nextGuess - currentGuess) >= abs(currentGuess - prevGuess)/2)) //condition 2
                || (!mflag && (abs(nextGuess - currentGuess) >= abs(prevGuess - oldGuess)/2)) //condition 3
                || (mflag && (abs(currentGuess - prevGuess) < minDelta)) //condition 4
                || (!mflag && (abs(nextGuess - currentGuess) < minDelta)) //condition 5
                )
        {
            //one of these was true so we have to use bisection
            nextGuess = (contrapoint + currentGuess) / 2;
            mflag = true;
        }
        else
        {
            mflag = false;
        }

        double sValue = getDistanceSlope(queryPoint, nextGuess);

        oldGuess = prevGuess;
        prevGuess = currentGuess;

        //choose the next a: if a and s have the same sign then s is the new a
        //if a and s have different signs then we keep the existing a
        if(sign(sValue) == sign(contrapointValue))
        {
            contrapoint = nextGuess;
            contrapointValue = sValue;
        }
        else
        {
            currentGuess = nextGuess;
            currentGuessValue = sValue;
        }

        //if a is a better guess than b, swap them - we want b to be better than a
        if(abs(contrapointValue) < abs(currentGuessValue))
        {
            std::swap(contrapoint, currentGuess);
            std::swap(contrapointValue, currentGuessValue);
        }
    }

    return currentGuess;
}

double SplineInverter::getDistanceSlope(const Vector3D &queryPoint, double t) const
{
    auto result = spline->getTangent(t);

	//get the displacement from the spline at T to the query point
	Vector3D displacement = result.position - queryPoint;

	//find projection of spline velocity onto displacement
    return Vector3D::dotProduct(displacement.normalized(), result.tangent);
}
