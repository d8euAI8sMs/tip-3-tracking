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
        void good_features_to_track(
            plot::world_t bb0,
            const cv::Mat & frame0,
            const cv::Mat & frame1,
            std::vector < geom::point2d_t > & pivots,
            const config & cfg) const
        {
            pivots.clear();

            plot::screen_t bb = {
                (int) std::ceil(bb0.xmin * frame0.cols) - 3,
                (int) std::ceil(bb0.xmax * frame0.cols) + 3,
                (int) std::ceil(bb0.ymin * frame0.rows) - 3,
                (int) std::ceil(bb0.ymax * frame0.rows) + 3
            };

            if (bb.width() <= 0 || bb.height() <= 0) return;

            std::vector < float > quality(bb.height() * bb.width());

            float maxquality = 0;

            for (size_t i = 0; i < bb.height(); ++i)
            for (size_t j = 0; j < bb.width(); ++j)
            {
                geom::point2d_t pivot = {
                    ((int) j + bb.xmin) / (float) frame0.cols,
                    ((int) i + bb.ymin) / (float) frame0.rows
                }, grad;
                quality[i * bb.width() + j] = track(
                    pivot, { 0, 0 }, frame0, frame1, grad, cfg);
                maxquality = std::fmax(maxquality, quality[i * bb.width() + j]);
            }

            for (size_t i = 0; i < bb.height(); ++i)
            for (size_t j = 0; j < bb.width(); ++j)
            {
                if (quality[i * bb.width() + j] < std::fmin(
                        cfg.eigenvalue_threshold, maxquality * cfg.feature_threshold))
                    continue;
                pivots.emplace_back(
                    ((int) j + bb.xmin) / (float) frame0.cols,
                    ((int) i + bb.ymin) / (float) frame0.rows
                );
            }
        }
        void track(
            std::vector < plot::world_t > & bb0,
            std::vector < std::pair < geom::point2d_t, geom::point2d_t > > & flow,
            const video_stream & s,
            const config & cfg) const
        {
            size_t rc = s.video.frames[0].rows;
            size_t cc = s.video.frames[0].cols;

            size_t maxmip = (size_t) std::fmax(std::floor(std::log2(std::fmin(rc, cc))) - 4, 0);

            size_t mip = std::fmin(maxmip, cfg.mip_count);

            std::vector < std::pair < cv::Mat, cv::Mat > > mips(mip);

            mips[0] = { s.video.frames[s.frame].clone(), s.video.frames[s.frame + 1].clone() };

            for (size_t i = 1; i < mip; ++i)
            {
                cv::Size mipsize(cc / (1 << i), rc / (1 << i));
                cv::resize(mips[0].first, mips[i].first, mipsize);
                cv::resize(mips[0].second, mips[i].second, mipsize);
            }

            if (cfg.blur)
            {
                for (size_t i = 0; i < mip; ++i)
                {
                    cv::medianBlur(mips[i].first, mips[i].first, 5);
                    cv::medianBlur(mips[i].second, mips[i].second, 5);
                    cv::GaussianBlur(mips[i].first, mips[i].first, cv::Size(5, 5), 0);
                    cv::GaussianBlur(mips[i].second, mips[i].second, cv::Size(5, 5), 0);
                }
            }

            flow.clear();

            for (size_t i = 0; i < bb0.size(); ++i)
            {
                std::vector < geom::point2d_t > pivots;

                good_features_to_track(bb0[i], mips.front().first, mips.front().second, pivots, cfg);

                std::vector < geom::point2d_t > grads(pivots.size());
                std::vector < geom::point2d_t > tmpgrads(pivots.size());
                std::vector < float > qual(pivots.size(), true);
                std::vector < bool > state(pivots.size(), true);

                for (size_t m = 0; m < mip; ++m)
                {
                    track(pivots, tmpgrads, grads, qual, state,
                          mips[mip - m - 1].first,
                          mips[mip - m - 1].second, cfg);
                }

                plot::world_t bbn = {
                    10000, 0,
                    10000, 0,
                };

                for (size_t p = 0; p < pivots.size(); ++p)
                {
                    if (!state[p]) continue;
                    auto pivot1 = pivots[p];
                    auto pivot2 = pivots[p] + grads[p];

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

                for (size_t p = 0; p < pivots.size(); ++p)
                {
                    if (!state[p]) continue;
                    flow.emplace_back(pivots[p], grads[p]);
                }
            }
        }
        void track(
            const std::vector < geom::point2d_t > & pc0,
            std::vector < geom::point2d_t > & tmpgrad,
            std::vector < geom::point2d_t > & grad,
            std::vector < float > & qual,
            std::vector < bool > & state,
            const cv::Mat & frame0,
            const cv::Mat & frame1,
            const config & cfg) const
        {
            geom::point2d_t newgrad;
            float maxqual = 0;
            for (size_t i = 0; i < pc0.size(); ++i)
            {
                if (!state[i]) continue;
                qual[i] = track(pc0[i], grad[i], frame0, frame1, newgrad, cfg);
                if (qual[i] < cfg.eigenvalue_threshold) state[i] = false;
                tmpgrad[i] = grad[i] + newgrad;
                maxqual = std::fmax(maxqual, qual[i]);
            }
            for (size_t i = 0; i < pc0.size(); ++i)
            {
                if (!state[i]) continue;
                if (qual[i] < maxqual * cfg.feature_threshold)
                {
                    state[i] = false;
                }
                else
                {
                    grad[i] = tmpgrad[i];
                }
            }
        }
        float track(
            const geom::point2d_t & pc0,
            const geom::point2d_t & grad,
            const cv::Mat & frame0,
            const cv::Mat & frame1,
            geom::point2d_t & out,
            const config & cfg) const
        {
            size_t wndsize = 13;

            float sx = wndsize / 3. / 2.;
            float sy = wndsize / 3. / 2.;

            cv::Mat wnd(wndsize, wndsize, CV_32F);

            for (size_t i = 0; i <= wndsize / 2; ++i)
            for (size_t j = 0; j <= wndsize / 2; ++j)
            {
                double d = std::exp(-1.0 * ((int)i * (int)i / sy / sy + (int)j * (int)j / sx / sx) / 2);
                wnd.at < float > (wndsize / 2 - (int)i, wndsize / 2 - (int)j) = (float) d;
                wnd.at < float > (wndsize / 2 + (int)i, wndsize / 2 + (int)j) = (float) d;
                wnd.at < float > (wndsize / 2 - (int)i, wndsize / 2 + (int)j) = (float) d;
                wnd.at < float > (wndsize / 2 + (int)i, wndsize / 2 - (int)j) = (float) d;
            }

            return track(pc0, grad, wnd, frame0, frame1, cfg, out);
        }
        float track(
            const geom::point2d_t & pc0,
            const geom::point2d_t & grad,
            const cv::Mat & wnd,
            const cv::Mat & frame0,
            const cv::Mat & frame1,
            const config & cfg,
            geom::point2d_t & out) const
        {
            cv::Mat ref, next;

            geom::point < int > pc = {
                (int) std::round(pc0.x * frame0.cols),
                (int) std::round(pc0.y * frame0.rows)
            };
            geom::point < int > pc2 = {
                (int) std::round((pc0.x + grad.x) * frame0.cols),
                (int) std::round((pc0.y + grad.y) * frame0.rows)
            };

            ref = frame0;
            next = frame1;

            cv::Mat a = cv::Mat::zeros(2, 2, CV_32FC1);
            cv::Mat b = cv::Mat::zeros(2, 1, CV_32FC1);
            
            for (size_t i = 0; i < wnd.rows; ++i)
            for (size_t j = 0; j < wnd.cols; ++j)
            {
                int r0 = pc.y - (int) i + wnd.rows / 2;
                int c0 = pc.x - (int) j + wnd.cols / 2;
                int r1 = pc2.y - (int) i + wnd.rows / 2;
                int c1 = pc2.x - (int) j + wnd.cols / 2;
                int ri[6] = { r0 - 1, r0, r0 + 1, r1 - 1, r1, r1 + 1 };
                int ci[6] = { c0 - 1, c0, c0 + 1, c1 - 1, c1, c1 + 1 };
                for (size_t k = 0; k < 6; ++k)
                {
                    if (ri[k] < 0) ri[k] = 0;
                    if (ri[k] >= ref.rows) ri[k] = ref.rows - 1;
                    if (ci[k] < 0) ci[k] = 0;
                    if (ci[k] >= ref.cols) ci[k] = ref.cols - 1;
                }
                // sobel operator to calculate more stable gradients
                float dy = (
                    2 * (ref.at < float > (ri[2], ci[1]) - ref.at < float > (ri[0], ci[1])) +
                    1 * (ref.at < float > (ri[2], ci[0]) - ref.at < float > (ri[0], ci[0])) +
                    1 * (ref.at < float > (ri[2], ci[2]) - ref.at < float > (ri[0], ci[2]))
                    ) / 4; /* possibly 8? */
                float dx = (
                    2 * (ref.at < float > (ri[1], ci[2]) - ref.at < float > (ri[1], ci[0])) +
                    1 * (ref.at < float > (ri[0], ci[2]) - ref.at < float > (ri[0], ci[0])) +
                    1 * (ref.at < float > (ri[2], ci[2]) - ref.at < float > (ri[2], ci[0]))
                    ) / 4; /* possibly 8? */
                float dt = next.at < float > (ri[3 + 1], ci[3 + 1]) - ref.at < float > (ri[1], ci[1]);
                float wi = wnd.at < float > (i, j);
                
                a.at < float > (0, 0) += wi * dx * dx;
                a.at < float > (1, 1) += wi * dy * dy;
                a.at < float > (1, 0) += wi * dx * dy;
                a.at < float > (0, 1) += wi * dx * dy;
                b.at < float > (0, 0) += - wi * dx * dt;
                b.at < float > (1, 0) += - wi * dy * dt;
            }

            std::vector < float > eigen;
            cv::eigen(a, eigen);

            cv::invert(a, a, cv::DECOMP_SVD);

            cv::Mat x = a * b;

            float dx = x.at < float > (0) / frame0.cols;
            float dy = x.at < float > (1) / frame0.rows;

            out = { dx, dy };

            return std::fmin(std::abs(eigen[0]), std::abs(eigen[1]));
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