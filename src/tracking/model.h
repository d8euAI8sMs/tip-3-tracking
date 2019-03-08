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
        static const size_t memory_threshold = (((1 << 10) << 10) << 10) >> 2; // 500Mb
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
            size_t pivot_neighborhood;
            size_t pivot_count;
            size_t mip_count;
            bool blur;
        };
        static config make_default_config()
        {
            return
            {
                {},
                1e-3,
                6,
                4,
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
            const video_stream & s,
            const config & cfg) const
        {
            size_t rc = s.video.frames[0].rows;
            size_t cc = s.video.frames[0].cols;

            size_t maxmip = (size_t) std::fmax(std::floor(std::log2(std::fmin(rc, cc))) - 4, 0);

            size_t mip = std::fmin(maxmip, cfg.mip_count);

            std::vector < std::pair < cv::Mat, cv::Mat > > mips(mip);

            mips[0] = { s.video.frames[s.frame], s.video.frames[s.frame + 1] };
            for (size_t i = 1; i < mip; ++i)
            {
                cv::Size mipsize(cc / (1 << i), rc / (1 << i));
                cv::resize(s.video.frames[s.frame], mips[i].first, mipsize);
                cv::resize(s.video.frames[s.frame + 1], mips[i].second, mipsize);
            }

            for (size_t i = 0; i < bb0.size(); ++i)
            {
                geom::point2d_t grad = { 0, 0 };
                for (size_t m = 0; m < mip; ++m)
                {
                    grad = grad + track(bb0[i], grad,
                                        mips[mip - m - 1].first,
                                        mips[mip - m - 1].second, cfg);
                }
                bb0[i] = { bb0[i].xmin + grad.x, bb0[i].xmax + grad.x,
                           bb0[i].ymin + grad.y, bb0[i].ymax + grad.y
                };
            }
        }
        geom::point2d_t track(
            const plot::world_t & bb0,
            const geom::point2d_t & grad,
            const cv::Mat & frame0,
            const cv::Mat & frame1,
            const config & cfg) const
        {
            plot::screen_t bb = {
                (int) std::ceil(bb0.xmin * frame0.cols),
                (int) std::ceil(bb0.xmax * frame0.cols),
                (int) std::ceil(bb0.ymin * frame0.rows),
                (int) std::ceil(bb0.ymax * frame0.rows)
            };
            geom::point < int > pc = {
                (int) std::round((bb0.xmax + bb0.xmin) / 2 * frame0.cols),
                (int) std::round((bb0.ymax + bb0.ymin) / 2 * frame0.rows)
            };
            geom::point < int > gc = {
                (int) std::round(grad.x * frame0.cols),
                (int) std::round(grad.y * frame0.rows)
            };

            std::vector < geom::point < int > > pivots(1 + cfg.pivot_count);
            pivots[0] = pc;
            for (size_t i = 1; i < cfg.pivot_count; ++i)
            {
                pivots[i] = {
                    pc.x + (rand() / (RAND_MAX + 1.0) - 0.5) * cfg.pivot_neighborhood,
                    pc.y + (rand() / (RAND_MAX + 1.0) - 0.5) * cfg.pivot_neighborhood
                };
            }

            float sx = bb.width() / 2.0 / 3.0;
            float sy = bb.height() / 2.0 / 3.0;

            size_t wndrows = (size_t) (3 * sy);
            size_t wndcols = (size_t) (3 * sx);
            if ((wndrows & 1) == 0) ++wndrows;
            if ((wndcols & 1) == 0) ++wndcols;
            if (wndrows < 3) wndrows = 3;
            if (wndcols < 3) wndcols = 3;

            cv::Mat wnd(wndrows, wndcols, CV_32F);

            for (size_t i = 0; i <= wndrows / 2; ++i)
            for (size_t j = 0; j <= wndcols / 2; ++j)
            {
                double d = std::exp(-1.0 * ((int)i * (int)i / sy / sy + (int)j * (int)j / sx / sx) / 2);
                wnd.at < float > (wndrows / 2 - (int)i, wndcols / 2 - (int)j) = (float) d;
                wnd.at < float > (wndrows / 2 + (int)i, wndcols / 2 + (int)j) = (float) d;
                wnd.at < float > (wndrows / 2 - (int)i, wndcols / 2 + (int)j) = (float) d;
                wnd.at < float > (wndrows / 2 + (int)i, wndcols / 2 - (int)j) = (float) d;
            }

            geom::point2d_t dd = { 0, 0 }, dt;
            size_t n = 0;

            for (size_t i = 0; i < pivots.size(); ++i)
            {
                if (!track(pivots[i], gc, wnd, frame0, frame1, cfg, dt)) continue;
                dd = dd + dt;
                ++n;
            }

            if (n > 0) dd = dd / (int) n;

            return dd;
        }
        bool track(
            plot::point < int > pc,
            plot::point < int > grad,
            const cv::Mat & wnd,
            const cv::Mat & frame0,
            const cv::Mat & frame1,
            const config & cfg,
            geom::point2d_t & out) const
        {
            cv::Mat ref, next;
            
            if (cfg.blur)
            {
                cv::GaussianBlur(frame0, ref, cv::Size(), wnd.cols / 3, wnd.rows / 3);
                cv::GaussianBlur(frame1, next, cv::Size(), wnd.cols / 3, wnd.rows / 3);
            }
            else
            {
                ref = frame0;
                next = frame1;
            }

            cv::Mat a = cv::Mat::zeros(2, 2, CV_32FC1);
            cv::Mat b = cv::Mat::zeros(2, 1, CV_32FC1);
            
            for (size_t i = 0; i < wnd.rows; ++i)
            for (size_t j = 0; j < wnd.cols; ++j)
            {
                int r0 = pc.y - (int) i + wnd.rows / 2;
                int c0 = pc.x - (int) j + wnd.cols / 2;
                int r1 = pc.y - (int) i + wnd.rows / 2 + grad.y;
                int c1 = pc.x - (int) j + wnd.cols / 2 + grad.x;
                int ri[6] = { r0 - 1, r0, r0 + 1, r1 - 1, r1, r1 + 1 };
                int ci[6] = { c0 - 1, c0, c0 + 1, c1 - 1, c1, c1 + 1 };
                for (size_t k = 0; k < 6; ++k)
                {
                    if (ri[k] < 0) ri[k] = 0;
                    if (ri[k] >= ref.rows) ri[k] = ref.rows - 1;
                    if (ci[k] < 0) ci[k] = 0;
                    if (ci[k] >= ref.cols) ci[k] = ref.cols - 1;
                }
                float dy = (ref.at < float > (ri[2], ci[1]) - ref.at < float > (ri[0], ci[1])) / 2;
                float dx = (ref.at < float > (ri[1], ci[2]) - ref.at < float > (ri[1], ci[0])) / 2;
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
            
            if (eigen[0] < cfg.eigenvalue_threshold ||
                eigen[1] < cfg.eigenvalue_threshold) return false;

            cv::invert(a, a, cv::DECOMP_SVD);

            cv::Mat x = a * b;

            float dx = x.at < float > (0) / frame0.cols;
            float dy = x.at < float > (1) / frame0.rows;

            out = { dx, dy };

            return true;
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