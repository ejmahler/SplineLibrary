#include "splineinverter.h"

#include "spline.h"
#include "spline_library/utils/splinesample_adaptor.h"
#include "spline_library/utils/nanoflann.hpp"
#include "spline_library/utils/optimization.h"

template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
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

    //define a lambda to compute the slope of the distance to the querypoint at T
    auto splineInstance = spline;
    auto distanceSlopeFunction = [splineInstance, queryPoint](double t) {
        auto result = splineInstance->getTangent(t);

        //get the displacement from the spline at T to the query point
        Vector3D displacement = result.position - queryPoint;

        //find projection of spline velocity onto displacement
        return Vector3D::dotProduct(displacement.normalized(), result.tangent);
    };

    double sampleDistanceSlope = distanceSlopeFunction(closestSampleT);

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
    double bValue = distanceSlopeFunction(b);

    //we know that the actual closest point is now between littleT and bigT
    //use the circle projection method to find the actual closest point, using the Ts as bounds
    return Optimization::brentsMethod(distanceSlopeFunction, a, aValue, b, bValue);
}
