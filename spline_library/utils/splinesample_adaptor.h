#ifndef SPLINESAMPLE_ADAPTOR_H
#define SPLINESAMPLE_ADAPTOR_H

/***********************************************************************
 * Software License Agreement (BSD License)
 *
 * Copyright 2011 Jose Luis Blanco (joseluisblancoc@gmail.com).
 *   All rights reserved.
 * Copyright 2014 Elliott Mahler (join.together@gmail.com).
 *   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *************************************************************************/

#include "nanoflann.hpp"
#include <vector>
#include "spline_library/vector3d.h"

struct SplineSamples3D
{
    typedef double coord_t; //!< The type of each coordinate

    struct Point3D
    {
        double x,y,z;
        double t;

        Point3D(double cx, double cy, double cz, double ct)
            :x(cx),y(cy),z(cz),t(ct)
        {}
    };

    std::vector<Point3D> pts;
};


template <typename Derived>
struct SplineSampleAdaptor
{
    typedef typename Derived::coord_t coord_t;

    const Derived obj;

    /// The constructor that sets the data set source
    SplineSampleAdaptor(const Derived &obj_) : obj(obj_) {}

    /// CRTP helper method
    inline const Derived& derived() const { return obj; }

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const { return derived().pts.size(); }

    // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
    inline coord_t kdtree_distance(const coord_t *p1, const size_t idx_p2,size_t size) const
    {
        const coord_t dx=p1[0]-derived().pts[idx_p2].x;
        const coord_t dy=p1[1]-derived().pts[idx_p2].y;
        const coord_t dz=p1[2]-derived().pts[idx_p2].z;
        return dx*dx+dy*dy+dz*dz;
    }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate value, the
    //  "if/else's" are actually solved at compile time.
    inline coord_t kdtree_get_pt(const size_t idx, int dim) const
    {
        if (dim==0) return derived().pts[idx].x;
        else if (dim==1) return derived().pts[idx].y;
        else return derived().pts[idx].z;
    }

    // Optional bounding-box computation: return false to default to a standard bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
    //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX &bb) const { return false; }
};

class SplineSampleTree
{
    typedef SplineSampleAdaptor<SplineSamples3D> AdaptorType;
    typedef nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<double, AdaptorType>,
            AdaptorType, 3>
        TreeType;

public:
    SplineSampleTree(const SplineSamples3D &samples)
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

#endif // SPLINESAMPLE_ADAPTOR_H
