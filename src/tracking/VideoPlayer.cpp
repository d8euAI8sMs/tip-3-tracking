// VideoPlayer.cpp : implementation file
//

#include "stdafx.h"
#include "tracking.h"
#include "VideoPlayer.h"


// CVideoPlayer

IMPLEMENT_DYNAMIC(CVideoPlayer, CPlotControl)

CVideoPlayer::CVideoPlayer()
{
    triple_buffered = true;
}

CVideoPlayer::~CVideoPlayer()
{
}


BEGIN_MESSAGE_MAP(CVideoPlayer, CPlotControl)
    ON_WM_LBUTTONDOWN()
    ON_WM_RBUTTONDOWN()
END_MESSAGE_MAP()



// CVideoPlayer message handlers

void CVideoPlayer::OnLButtonDown(UINT nFlags, CPoint point)
{
    if (callback)
    {
        CRect r; GetClientRect(&r);
        callback({ (float) point.x / r.Width(), (float) point.y / r.Height() }, true);
    }

    CPlotControl::OnLButtonDown(nFlags, point);
}

void CVideoPlayer::OnRButtonDown(UINT nFlags, CPoint point)
{
    if (callback)
    {
        CRect r; GetClientRect(&r);
        callback({ (float) point.x / r.Width(), (float) point.y / r.Height() }, false);
    }

    CPlotControl::OnRButtonDown(nFlags, point);
}
