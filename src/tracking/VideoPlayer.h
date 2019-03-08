#pragma once

#include <util/common/gui/PlotControl.h>
#include "model.h"

// CVideoPlayer

class CVideoPlayer : public CPlotControl
{
	DECLARE_DYNAMIC(CVideoPlayer)

public:
	CVideoPlayer();
	virtual ~CVideoPlayer();

protected:
	DECLARE_MESSAGE_MAP()
    const model::model_data * m_data;
public:
    std::function < void (const geom::point2d_t &, bool lbutton) > callback;
    void Init(const model::model_data & data)
    {
        m_data = &data;
        plot_layer.with(model::make_bmp_plot(m_data->source));
        plot_layer.with(model::make_decorator_plot(m_data->decorator));
    }
    afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
    afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
};


