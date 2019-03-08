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
    };

    inline parameters make_default_parameters()
    {
        parameters p =
        {
            480, 360
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

        video_data & to_cbitmap(CBitmap & bmp, int frame)
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
            return *this;
        }
        video_data & from_file(const std::string & path, cv::Size resolution = cv::Size())
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
        std::vector < geom::point < int > > markers;
        std::vector < plot::screen_t > bounding_boxes;
    };

    inline plot::drawable::ptr_t make_bmp_plot(video_stream & s)
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

    inline plot::drawable::ptr_t make_decorator_plot(decorator_data & d)
    {
        return plot::custom_drawable::create([&d] (CDC & dc, const plot::viewport & vp)
        {
            if (d.decorator_mask & decorator_data::decorator_markers)
            {
                auto brush = plot::palette::brush(RGB(100,100,255));
                for each (auto & marker in d.markers)
                {
                    CRect rect(marker - geom::point < int > (3, 3),
                               marker + geom::point < int > (3, 3));
                    dc.FillRect(rect, brush.get());
                }
            }
            if (d.decorator_mask & decorator_data::decorator_bbox)
            {
                auto pen = plot::palette::pen(RGB(255,100,100), 3);
                dc.SelectObject(pen.get());
                for each (auto & bbox in d.bounding_boxes)
                {
                    CRect rect(bbox.xmin, bbox.ymin, bbox.xmax, bbox.ymax);
                    dc.Rectangle(rect);
                }
            }
        });
    }

    /*****************************************************/
    /*                     algo                          */
    /*****************************************************/

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