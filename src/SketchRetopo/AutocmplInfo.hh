#pragma once
#include <Eigen/Core>
#include <boost/optional.hpp>
#include <vector>
#include <kt84/geometry/PointNormal.hh>

struct AutocmplInfo {
    struct UVLine {
        Eigen::Vector2d uv0;
        Eigen::Vector2d uv1;
        boost::optional<Eigen::Vector2d> unnormalized_dir0;
        boost::optional<Eigen::Vector2d> unnormalized_dir1;
        
        UVLine() {}
        UVLine(const Eigen::Vector2d& uv0_,
               const Eigen::Vector2d& uv1_,
               const boost::optional<Eigen::Vector2d>& unnormalized_dir0_ = boost::none,
               const boost::optional<Eigen::Vector2d>& unnormalized_dir1_ = boost::none)
            : uv0(uv0_)
            , uv1(uv1_)
            , unnormalized_dir0(unnormalized_dir0_)
            , unnormalized_dir1(unnormalized_dir1_)
        {}
    };
    
    std::vector<UVLine> uv_lines;
    std::vector<kt84::PointNormal> corners;
    
    void add_uv_line(const Eigen::Vector2d& uv0,
                     const Eigen::Vector2d& uv1,
                     const boost::optional<Eigen::Vector2d>& unnormalized_dir0 = boost::none,
                     const boost::optional<Eigen::Vector2d>& unnormalized_dir1 = boost::none)
    {
        uv_lines.push_back(UVLine(uv0, uv1, unnormalized_dir0, unnormalized_dir1));
    }
    void set_corners(const kt84::PointNormal& pn0, const kt84::PointNormal& pn1, const kt84::PointNormal& pn2, const kt84::PointNormal& pn3) {
        corners.resize(4);
        corners[0] = pn0;
        corners[1] = pn1;
        corners[2] = pn2;
        corners[3] = pn3;
    }
    void set_corners(const kt84::PointNormal& pn0, const kt84::PointNormal& pn1, const kt84::PointNormal& pn2) {
        corners.resize(3);
        corners[0] = pn0;
        corners[1] = pn1;
        corners[2] = pn2;
    }
};
