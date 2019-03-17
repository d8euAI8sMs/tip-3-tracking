// trackingDlg.h : header file
//

#pragma once

#include <util/common/gui/SimulationDialog.h>
#include <util/common/gui/PlotControl.h>
#include "afxcmn.h"

#include "model.h"
#include "VideoPlayer.h"

// CTrackingDlg dialog
class CTrackingDlg : public CSimulationDialog
{
// Construction
public:
    CTrackingDlg(CWnd* pParent = NULL);    // standard constructor

// Dialog Data
    enum { IDD = IDD_TRACKING_DIALOG };

    protected:
    virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
    HICON m_hIcon;

    // Generated message map functions
    virtual BOOL OnInitDialog();
    afx_msg void OnPaint();
    afx_msg HCURSOR OnQueryDragIcon();
    DECLARE_MESSAGE_MAP()
public:
    model::model_data m_data;
    model::tracker::config m_cfg;
    CVideoPlayer m_videoCtrl;
    afx_msg void OnBnClickedButton1();
    afx_msg void OnBnClickedButton2();
    CSliderCtrl m_frameSlider;
    afx_msg void OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
    afx_msg void OnBnClickedButton3();
    void OnSimulation() override;
    BOOL m_bBlur;
};
