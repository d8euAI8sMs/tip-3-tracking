#pragma once

#include <afxwin.h>

#include <vector>
#include <map>
#include <cstdint>

#include <util/common/math/complex.h>
#include <util/common/plot/plot.h>
#include <util/common/math/fft.h>
#include <util/common/math/vec.h>
#include <util/common/math/raster.h>
#include <util/common/geom/geom.h>

#include <opencv2/opencv.hpp>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif // !M_PI
#ifndef M_E
#define M_E 2.7182818284590452353602874713527
#endif // !M_E

namespace model
{

    /*****************************************************/
    /*                     params                        */
    /*****************************************************/

    struct parameters
    {
        size_t target_resolution_w;
        size_t target_resolution_h;
        size_t max_frames;
    };

    inline parameters make_default_parameters()
    {
        parameters p =
        {
            480, 360,
            500
        };
        return p;
    }

    /*****************************************************/
    /*                     data                          */
    /*****************************************************/

    struct video_data
    {
        static const size_t memory_threshold = (((1 << 10) << 10) << 10); // 1Gb
        std::vector < cv::Mat > frames;

        void to_cbitmap(CBitmap & bmp, int frame) const
        {
            auto & mat = frames[frame];
            std::vector < COLORREF > buf(mat.rows * mat.cols);
            for (size_t i = 0; i < mat.rows; ++i)
            for (size_t j = 0; j < mat.cols; ++j)
            {
                if (mat.channels() == 1)
                {
                    float v = mat.at < float > (i, j);
                    BYTE c = (BYTE) (v * 255);
                    buf[mat.cols * i + j] = RGB(c, c, c);
                }
                else
                {
                    auto v = mat.at < cv::Vec3f > (i, j);
                    BYTE c[3] = { (BYTE) (v[0] * 255), (BYTE) (v[1] * 255), (BYTE) (v[2] * 255) };
                    buf[mat.cols * i + j] = RGB(c[0], c[1], c[2]);
                }
            }
            bmp.DeleteObject();
            bmp.CreateBitmap(mat.cols, mat.rows, 1, sizeof(COLORREF) * 8, (LPCVOID) buf.data());
            bmp.SetBitmapDimension(mat.cols, mat.rows);
        }
        video_data & from_file(const std::string & path, size_t frame_limit = 0, cv::Size resolution = cv::Size())
        {
            cv::VideoCapture cap(path);
            frames.clear();
            size_t size = 0;
            cv::Mat frame, copy;
            while (cap.read(frame))
            {
                cv::cvtColor(frame, copy, cv::COLOR_BGR2GRAY);
                if (resolution.width != 0) cv::resize(copy, copy, resolution);
                copy.convertTo(copy, CV_32F, 1. / 255);
                size += sizeof(float) * copy.rows * copy.cols;
                if (size > memory_threshold) break;
                if (frame_limit != 0 && frames.size() >= frame_limit) break;
                frames.push_back(std::move(copy));
            }
            return *this;
        }
    };

    struct video_stream
    {
        video_data video;
        int frame = -1;
    };

    struct decorator_data
    {
        static const size_t decorator_markers = 1 << 0;
        static const size_t decorator_bbox    = 1 << 1;
        size_t decorator_mask = 0xff;
        std::vector < geom::point2d_t > markers;
        std::vector < plot::world_t > bounding_boxes;
        std::vector < std::pair < geom::point2d_t, geom::point2d_t > > opflow;
    };

    inline plot::drawable::ptr_t make_bmp_plot(const video_stream & s)
    {
        return plot::custom_drawable::create([&s] (CDC & dc, const plot::viewport & vp)
        {
            if (s.frame == -1) return;
            CBitmap b; s.video.to_cbitmap(b, s.frame);
            CDC memDC; memDC.CreateCompatibleDC(&dc);
            memDC.SelectObject(&b);
            dc.SetStretchBltMode(HALFTONE);
            auto wh = b.GetBitmapDimension();
            dc.StretchBlt(vp.screen.xmin, vp.screen.ymin,
                          vp.screen.width(), vp.screen.height(),
                          &memDC, 0, 0, wh.cx, wh.cy, SRCCOPY);
        });
    }

    inline plot::drawable::ptr_t make_decorator_plot(const decorator_data & d)
    {
        return plot::custom_drawable::create([&d] (CDC & dc, const plot::viewport & vp)
        {
            if (d.decorator_mask & decorator_data::decorator_markers)
            {
                auto brush = plot::palette::brush(RGB(100,255,100));
                for each (auto & marker in d.markers)
                {
                    geom::point2d_t m(marker.x * vp.screen.width(),
                                      marker.y * vp.screen.height());
                    CRect rect(m - geom::point2d_t(4, 4),
                               m + geom::point2d_t(4, 4));
                    dc.FillRect(rect, brush.get());
                }
            }
            if (d.decorator_mask & decorator_data::decorator_bbox)
            {
                auto pen = plot::palette::pen(RGB(255,100,100), 3);
                dc.SelectObject(pen.get());
                for each (auto & bbox in d.bounding_boxes)
                {
                    plot::screen_t b = {
                        (int) (bbox.xmin * vp.screen.width()),
                        (int) (bbox.xmax * vp.screen.width()),
                        (int) (bbox.ymin * vp.screen.height()),
                        (int) (bbox.ymax * vp.screen.height())
                    };
                    dc.MoveTo({ b.xmin, b.ymin });
                    dc.LineTo({ b.xmin, b.ymax });
                    dc.LineTo({ b.xmax, b.ymax });
                    dc.LineTo({ b.xmax, b.ymin });
                    dc.LineTo({ b.xmin, b.ymin });
                }
            }
            if (true)
            {
                auto pen = plot::palette::pen(RGB(0,100,255), 1);
                dc.SelectObject(pen.get());
                for each (auto & v in d.opflow)
                {
                    auto p1 = v.first;
                    auto p2 = p1 + v.second;
                    dc.MoveTo(CPoint(p1.x * vp.screen.width(), p1.y * vp.screen.height()));
                    dc.LineTo(CPoint(p2.x * vp.screen.width(), p2.y * vp.screen.height()));
                }
            }
        });
    }

    /*****************************************************/
    /*                     algo                          */
    /*****************************************************/

    class tracker
    {
    public:
        struct config
        {
            std::vector < geom::point < int > > trackpoints;
            float eigenvalue_threshold;
            float feature_threshold;
            size_t mip_count;
            bool blur;
        };
        static config make_default_config()
        {
            return
            {
                {},
                0.01,
                0.04,
                4,
                false
            };
        }
    private:
        const parameters & params;
    public:
        tracker(const parameters & params) : params(params)
        {
        }
    public:
        void track(
            std::vector < plot::world_t > & bb0,
            std::vector < std::pair < geom::point2d_t, geom::point2d_t > > & flow,
            const video_stream & s,
            const config & cfg) const
        {
            auto ref = s.video.frames[s.frame].clone();
            auto next = s.video.frames[s.frame + 1].clone();

            ref.convertTo(ref, CV_8U, 255);
            next.convertTo(next, CV_8U, 255);

            flow.clear();

            for (size_t i = 0; i < bb0.size(); ++i)
            {
                plot::screen_t bb = {
                    (int) std::ceil(bb0[i].xmin * ref.cols) - 5,
                    (int) std::ceil(bb0[i].xmax * ref.cols) + 5,
                    (int) std::ceil(bb0[i].ymin * ref.rows) - 5,
                    (int) std::ceil(bb0[i].ymax * ref.rows) + 5
                };

                if (bb.xmin < 0) bb.xmin = 0;
                if (bb.ymin < 0) bb.ymin = 0;
                if (bb.xmax >= ref.cols) bb.xmax = ref.cols - 1;
                if (bb.ymax >= ref.rows) bb.ymax = ref.rows - 1;

                if (bb.width() <= 5 || bb.height() <= 5) continue;
               
                std::vector<cv::Point2f> points[2];
                cv::Mat status, err;
                cv::goodFeaturesToTrack(
                    ref.colRange(bb.xmin, bb.xmax).rowRange(bb.ymin, bb.ymax),
                    points[0], 1000, cfg.feature_threshold, 1);

                if (points[0].size() < 5) continue;

                cv::calcOpticalFlowPyrLK(
                    ref.colRange(bb.xmin, bb.xmax).rowRange(bb.ymin, bb.ymax),
                    next.colRange(bb.xmin, bb.xmax).rowRange(bb.ymin, bb.ymax),
                    points[0], points[1], status, err);

                plot::world_t bbn = {
                    10000, 0,
                    10000, 0,
                };

                for (size_t p = 0; p < points[0].size(); ++p)
                {
                    if (status.at < std::uint8_t > (p, 0) == 0) continue;
                    auto pivot1 = points[0][p];
                    auto pivot2 = points[1][p];
                    pivot1.x = bb0[i].xmin + pivot1.x / ref.cols;
                    pivot1.y = bb0[i].ymin + pivot1.y / ref.rows;
                    pivot2.x = bb0[i].xmin + pivot2.x / ref.cols;
                    pivot2.y = bb0[i].ymin + pivot2.y / ref.rows;
                
                    bbn.xmin = std::fmin(bbn.xmin, pivot1.x);
                    bbn.xmax = std::fmax(bbn.xmax, pivot1.x);
                    bbn.ymin = std::fmin(bbn.ymin, pivot1.y);
                    bbn.ymax = std::fmax(bbn.ymax, pivot1.y);
                
                    bbn.xmin = std::fmin(bbn.xmin, pivot2.x);
                    bbn.xmax = std::fmax(bbn.xmax, pivot2.x);
                    bbn.ymin = std::fmin(bbn.ymin, pivot2.y);
                    bbn.ymax = std::fmax(bbn.ymax, pivot2.y);
                }
                
                bb0[i] = bbn;
                
                for (size_t p = 0; p < points[0].size(); ++p)
                {
                    auto pivot1 = points[0][p];
                    auto pivot2 = points[1][p];
                    pivot1.x = bb0[i].xmin + pivot1.x / ref.cols;
                    pivot1.y = bb0[i].ymin + pivot1.y / ref.rows;
                    pivot2.x = bb0[i].xmin + pivot2.x / ref.cols;
                    pivot2.y = bb0[i].ymin + pivot2.y / ref.rows;
                    if (status.at < std::uint8_t > (p, 0) == 0) continue;
                    flow.emplace_back(geom::point2d_t{ pivot1.x, pivot1.y },
                                      geom::point2d_t{ pivot2.x - pivot1.x, pivot2.y - pivot1.y });
                }
            }
        }
    };

    /*****************************************************/
    /*                     model                         */
    /*****************************************************/

    struct model_data
    {
        parameters params;
        video_stream source;
        decorator_data decorator;
    };
}